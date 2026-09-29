#!/usr/bin/env python3
# NOTE: 本脚本源自项目内部流水线, 含两台计算机的绝对路径。
# 复现时请全局替换:
#   /data/wcf/MD_trend_analysis  -> 本地仓库根 (analysis 数据目录的父目录)
#   /data/md_extend              -> 远程机工作根
#   /home/rcdb/anaconda3/envs/AmberTools23 -> 远程 AMBERHOME (AmberTools23)
#   /data/Tools/Amber22          -> 本地 AMBERHOME (Amber22)
#   rcdb@100.76.253.60 / rcdb@100.97.18.74:10021 -> 计算机地址(需自备 ssh key)
"""
批量趋势收敛分析 (配合 13_batch_mmgbsa_all.py 的 500ns 重算):

每任务输出:
  1. rmsd_bb.dat       — 延伸全程骨架 RMSD (stride 10, 0.2ns/点, 参考=延伸首帧)
  2. dG_segK.dat       — 每段末 20ns (100帧) MM/GBSA → dG 随累计采样时间演化
  3. trend_summary.json — RMSD 100ns 窗口统计 + 收敛判定 + dG 演化

协议与 12 脚本一致: igb=5, saltcon=0.154, interval=20.
用法:  python3 14_batch_trend_analysis.py [--local] [--jobs 8] [--status]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


def get_paths(is_local: bool) -> dict:
    if is_local:
        return {
            "runs_root": Path("/data/wcf/MD_trend_analysis/md_extend/runs_nvt"),
            "out_root": Path("/data/wcf/MD_trend_analysis/extend_analysis"),
            "cpptraj": "/data/Tools/Amber22/bin/cpptraj",
            "mmpbsa": "/data/Tools/Amber22/bin/MMPBSA.py",
            "py3": sys.executable,
            "amberhome_env": "/data/Tools/Amber22",
        }
    return {
        "runs_root": Path("/data/md_extend/runs"),
        "out_root": Path("/data/md_extend/extend_analysis"),
        "cpptraj": "/data/Tools/Amber22/bin/cpptraj",
        "mmpbsa": "/home/rcdb/anaconda3/envs/AmberTools23/bin/MMPBSA.py",
        "py3": "/home/rcdb/anaconda3/envs/AmberTools23/bin/python",
        "amberhome_env": "/home/rcdb/anaconda3/envs/AmberTools23",
    }


IS_LOCAL = "local" in " ".join(sys.argv)
PATHS = get_paths(IS_LOCAL)
RUNS_ROOT, OUT_ROOT = PATHS["runs_root"], PATHS["out_root"]
CPPTRAJ, MMPBSA, PY3 = PATHS["cpptraj"], PATHS["mmpbsa"], PATHS["py3"]
AMBERHOME_ENV = PATHS["amberhome_env"]

STRIP_MASK = ":WAT,Na+,Cl-,K+,Cs+"
MASKS = {
    "SH3_VNAR": {"receptor": ":1-105", "ligand": ":106-211"},
    "HCG_VNAR": {"receptor": ":1-281", "ligand": ":282-414"},
}


def last_seg_done(rundir: Path, max_seg: int = 6) -> int:
    for k in range(max_seg, 0, -1):
        out, rst = rundir / f"md_seg{k}.out", rundir / f"md_seg{k}.rst"
        if out.exists() and rst.exists():
            tail = out.read_text(errors="ignore")[-2000:]
            if "STOP" in tail or "Final" in tail:
                return k
    return 0


def run_cmd(cmd, logfile, cwd=None):
    with open(logfile, "w") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd)
    return p.returncode


def calc_rmsd(task_dir: Path, system: str, nseg: int) -> str:
    outdir = OUT_ROOT / system / task_dir.name
    outdir.mkdir(parents=True, exist_ok=True)
    dat = outdir / "rmsd_bb.dat"
    if dat.exists():
        # 缓存校验: 点数须覆盖全部段 (0.1ns/点, 每段1000点), 否则重算
        n_data = sum(1 for l in open(dat) if not l.startswith("#"))
        if n_data >= nseg * 1000:
            return "skip"
        dat.unlink()
    cin = outdir / "_trend_rmsd.in"
    lines = [f"parm {task_dir/'system.prmtop'}"]
    for k in range(1, nseg + 1):
        lines.append(f"trajin {task_dir/f'md_seg{k}.nc'} 1 last 10")
    lines += [
        f"strip {STRIP_MASK}",
        "autoimage",
        "rms bb first out rmsd_bb.dat",
        "run",
        "quit",
    ]
    cin.write_text("\n".join(lines) + "\n")
    rc = run_cmd([CPPTRAJ, "-i", str(cin)], outdir / "_trend_rmsd.log", cwd=outdir)
    return "ok" if rc == 0 else "cpptraj_err"


def calc_dg_seg(task_dir: Path, system: str, k: int) -> float | None:
    outdir = OUT_ROOT / system / task_dir.name
    outdir.mkdir(parents=True, exist_ok=True)
    masks = MASKS[system]
    nowat = outdir / f"md_last20ns_seg{k}.nc"
    final = outdir / f"dG_seg{k}.dat"

    # MM/GBSA 已算过则直接解析
    if final.exists():
        m = re.search(r"DELTA TOTAL\s+(-?[\d.]+)\s+(-?[\d.]+)",
                      final.read_text(errors="ignore"))
        if m:
            return float(m.group(1))

    if not nowat.exists():
        cin = outdir / f"_prep_seg{k}.in"
        cin.write_text(
            f"parm {task_dir/'system.prmtop'}\n"
            f"trajin {task_dir/f'md_seg{k}.nc'} 8001 last\n"
            f"strip {STRIP_MASK}\n"
            "autoimage\n"
            f"trajout {nowat}\nrun\nquit\n"
        )
        rc = run_cmd([CPPTRAJ, "-i", str(cin)], outdir / f"_prep_seg{k}.log",
                     cwd=outdir)
        if rc != 0:
            return None

    # 复用 13 脚本的拓扑
    prot = outdir / "protein_only.prmtop"
    if not prot.exists():
        return None

    cin = outdir / f"_mmpbsa_seg{k}.in"
    cin.write_text(
        "&general\n"
        "  startframe=1, endframe=2000, interval=20,\n"
        "  verbose=1, entropy=0,\n"
        f"  receptor_mask='{masks['receptor']}', ligand_mask='{masks['ligand']}',\n"
        "/\n&gb\n  igb=5, saltcon=0.154,\n/\n"
    )
    cmd = " ".join([
        MMPBSA, "-O", "-i", str(cin), "-o", str(final),
        "-cp", str(prot),
        "-rp", str(outdir / "receptor_mmpbsa.prmtop"),
        "-lp", str(outdir / "ligand_mmpbsa.prmtop"),
        "-y", str(nowat),
    ])
    with open(outdir / f"_mmpbsa_seg{k}.log", "w") as lf:
        subprocess.run(["bash", "-c", f"export AMBERHOME={AMBERHOME_ENV} && {cmd}"],
                       stdout=lf, stderr=subprocess.STDOUT, cwd=outdir)
    if not final.exists():
        return None
    m = re.search(r"DELTA TOTAL\s+(-?[\d.]+)\s+(-?[\d.]+)",
                  final.read_text(errors="ignore"))
    return float(m.group(1)) if m else None


def window_stats(v: np.ndarray, n_win: int, win_ns: float = 100.0) -> list:
    n_pt = int(win_ns / 0.1)  # 每 100ns 窗的点数 (0.1ns/点: ntwx=5000×dt=0.002ps=10ps/帧, stride 10)
    stats = []
    for w in range(n_win):
        seg = v[w * n_pt:(w + 1) * n_pt]
        if len(seg) < 10:
            stats.append(None)
            continue
        stats.append({
            "window": f"{w*100:.0f}-{(w+1)*100:.0f}ns",
            "mean": float(seg.mean()), "std": float(seg.std()), "n": len(seg)})
    return stats


def analyze_task(task_dir: Path, system: str) -> dict:
    name = task_dir.name
    outdir = OUT_ROOT / system / name
    outdir.mkdir(parents=True, exist_ok=True)
    nseg = last_seg_done(task_dir)
    rec = {"system": system, "name": name, "nseg": nseg,
           "total_ext_ns": nseg * 100}

    # 1) RMSD
    st = calc_rmsd(task_dir, system, nseg)
    dat = outdir / "rmsd_bb.dat"
    if dat.exists():
        arr = np.loadtxt(dat, comments="#")
        v = arr[:, 1]
        rec["rmsd"] = {
            "mean_last100ns": float(v[-500:].mean()),
            "max": float(v.max()),
            "windows": window_stats(v, nseg),
            "converged": None,
        }
        wins = rec["rmsd"]["windows"]
        if None not in wins and len(wins) >= 2:
            # 旧判据: 末两窗漂移 + 末窗斜率 (5% 或 0.1Å 绝对阈值)
            m1, m2 = wins[-2]["mean"], wins[-1]["mean"]
            drift = abs(m2 - m1)
            slope = (v[-500:] - v[-1000:-500]).mean() if len(v) >= 1000 else 0.0
            thr = max(0.05 * abs(m1), 0.1)
            rec["rmsd"]["converged"] = bool(drift < thr and abs(slope) < thr)
            rec["rmsd"]["drift_last2win"] = float(drift)
        # 新判据 (2026-09-06 采纳): 末1/3轨迹 std<2Å 且 |线性斜率|<1Å/100ns
        n3 = max(len(v) // 3, 10)
        tv, vv = arr[-n3:, 0], v[-n3:]
        lin_slope = np.polyfit(tv * 0.1, vv, 1)[0] * 100.0  # Å/100ns (0.1ns/点)
        rec["rmsd"]["tail_std"] = float(vv.std())
        rec["rmsd"]["tail_slope_per100ns"] = float(lin_slope)
        rec["rmsd"]["converged_v2"] = bool(
            vv.std() < 2.0 and abs(lin_slope) < 1.0)

    # 2) dG 演化 (每段末20ns)
    dgs, dgs_std = [], []
    for k in range(1, nseg + 1):
        dg = calc_dg_seg(task_dir, system, k)
        dgs.append(dg)
        final = outdir / f"dG_seg{k}.dat"
        m = re.search(r"DELTA TOTAL\s+(-?[\d.]+)\s+(-?[\d.]+)",
                      final.read_text(errors="ignore")) if final.exists() else None
        dgs_std.append(float(m.group(2)) if m else None)
    rec["dG_evolution"] = {"seg_end_ns": [100 + k * 100 for k in range(1, nseg + 1)],
                           "dG": dgs, "dG_std": dgs_std}
    if nseg >= 2 and dgs[-1] is not None and dgs[-2] is not None:
        drift = abs(dgs[-1] - dgs[-2])
        # 阈值 10% (2026-09-16 由 5% 上调): 体系 ΔG 量级 -100~-120 kcal/mol 时,
        # 5% 仅 5-6 kcal/mol, 与 MM/GBSA 段间 SD(10-12) 同量级 -> 会把方法噪声误判为未收敛。
        thr = max(0.10 * abs(dgs[-2]), 1.0)
        rec["dG_evolution"]["drift_last2seg"] = drift
        rec["dG_evolution"]["converged"] = bool(drift < thr)
        rec["dG_evolution"]["converged_thr_frac"] = 0.10
        rec["dG_evolution"]["converged_thr"] = float(thr)

    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--local", action="store_true")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    tasks = []
    for system in MASKS:
        root = RUNS_ROOT / system
        if not root.is_dir():
            continue
        for d in sorted(root.iterdir()):
            if d.is_dir() and (d / "system.prmtop").exists():
                seg = last_seg_done(d)
                if seg > 0:
                    tasks.append((d, system))

    print(f"任务: {len(tasks)}")
    if args.status:
        for d, s in tasks:
            print(f"  {s}/{d.name}: seg={last_seg_done(d)}")
        return

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(analyze_task, d, s): (d, s) for d, s in tasks}
        for f in as_completed(futs):
            r = f.result()
            results.append(r)
            rm = r.get("rmsd", {})
            dg = r.get("dG_evolution", {}).get("dG", [])
            print(f"[{r['name']}] nseg={r['nseg']} "
                  f"rmsd_conv={rm.get('converged')} "
                  f"last100ns_mean={rm.get('mean_last100ns')} "
                  f"dG_evo={['%.1f' % x if x else 'NA' for x in dg]}",
                  flush=True)

    outj = OUT_ROOT / "trend_summary.json"
    outj.write_text(json.dumps(results, indent=1))
    print(f"\n汇总 -> {outj}")


if __name__ == "__main__":
    main()
