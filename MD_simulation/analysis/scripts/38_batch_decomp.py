#!/usr/bin/env python3
# NOTE: 本脚本源自项目内部流水线, 含两台计算机的绝对路径。
# 复现时请全局替换:
#   /data/wcf/MD_trend_analysis  -> 本地仓库根 (analysis 数据目录的父目录)
#   /data/md_extend              -> 远程机工作根
#   /home/rcdb/anaconda3/envs/AmberTools23 -> 远程 AMBERHOME (AmberTools23)
#   /data/Tools/Amber22          -> 本地 AMBERHOME (Amber22)
#   rcdb@100.76.253.60 / rcdb@100.97.18.74:10021 -> 计算机地址(需自备 ssh key)
"""Step 38: 补跑 per-residue decomposition (MM/GBSA + &decomp)。

口径与 500 ns 统一重算 (12 脚本 FINAL_RESULTS_MMGBSA_BINDING_500ns.dat) 完全一致:
  轨迹 = 各体系 md_last20ns_nowat.nc (seg4 末 20ns, 2000 帧)
  采样 = startframe=1, endframe=2000, interval=20 → 100 帧
  GB   = igb=5, saltcon=0.154
  新增 = &decomp idecomp=1, print_res="within 4"
拓扑复用 12/14 脚本产物: protein_only.prmtop / receptor_mmpbsa.prmtop / ligand_mmpbsa.prmtop
输出: extend_analysis/{SYSTEM}/{name}/DECOMP_MMGBSA_500ns.dat (FINAL_DECOMP 格式)

用法:
  python3 38_batch_decomp.py --local  [--jobs 4] [--only NAME] [--status]
  python3 38_batch_decomp.py --remote [--jobs 4] [--only NAME] [--status]
  (远程作业在 100.76.253.60 上执行, 本脚本通过 ssh 提交; 远程机需已放好本脚本)
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

IS_LOCAL = "--local" in sys.argv or "--remote" not in sys.argv
if IS_LOCAL:
    OUT_ROOT = Path("/data/wcf/MD_trend_analysis/extend_analysis")
    MMPBSA = "/data/Tools/Amber22/bin/MMPBSA.py"
    AMBERHOME = "/data/Tools/Amber22"
else:
    OUT_ROOT = Path("/data/md_extend/extend_analysis")
    MMPBSA = "/home/rcdb/anaconda3/envs/AmberTools23/bin/MMPBSA.py"
    AMBERHOME = "/home/rcdb/anaconda3/envs/AmberTools23"

MASKS = {
    "SH3_VNAR": {"receptor": ":1-105", "ligand": ":106-211"},
    "HCG_VNAR": {"receptor": ":1-281", "ligand": ":282-392"},
}


def jobs_from_finals(only: str | None) -> list[tuple[str, Path]]:
    """扫描已完成 500ns 统一重算的体系 (两种 FINAL 文件名都认), 且具备全部前置."""
    out = []
    for system in MASKS:
        for d in sorted((OUT_ROOT / system).iterdir()):
            if not d.is_dir():
                continue
            if only and only not in d.name:
                continue
            finals = [d / "FINAL_RESULTS_MMGBSA_BINDING_500ns.dat",
                      d / "FINAL_MMGBSA_500ns.dat"]
            if not any(f.exists() for f in finals):
                continue
            pre = ["md_last20ns_nowat.nc", "protein_only.prmtop",
                   "receptor_mmpbsa.prmtop", "ligand_mmpbsa.prmtop"]
            missing = [f for f in pre if not (d / f).exists()]
            if missing:
                print(f"SKIP {system}/{d.name}: missing {missing}")
                continue
            out.append((system, d))
    return out


def hist_masks(d: Path) -> tuple[str, str] | None:
    """从该体系历史 MMPBSA 输入提取掩码 (对该体系拓扑必然有效), 按优先级。"""
    for f in ("mmpbsa_gb_extend.in", "_mmpbsa_seg4.in", "_mmpbsa_seg1.in", "mmpbsa_gb.in"):
        p = d / f
        if not p.exists():
            continue
        t = p.read_text(errors="ignore")
        mr = re.search(r"receptor_mask='([^']+)'", t)
        ml = re.search(r"ligand_mask='([^']+)'", t)
        if mr and ml:
            return mr.group(1), ml.group(1)
    return None


def run_one(system: str, d: Path) -> tuple[str, str]:
    name = d.name
    dec = d / "DECOMP_MMGBSA_500ns.dat"
    if dec.exists():
        return name, "skip"
    masks = hist_masks(d) or (MASKS[system]["receptor"], MASKS[system]["ligand"])
    cin = d / "_mmpbsa_decomp.in"
    cin.write_text(
        "&general\n"
        "  startframe=1, endframe=2000, interval=20,\n"
        "  verbose=1, entropy=0,\n"
        f"  receptor_mask='{masks[0]}', ligand_mask='{masks[1]}',\n"
        "/\n"
        "&gb\n"
        "  igb=5, saltcon=0.154,\n"
        "/\n"
        "&decomp\n"
        # 注意: 本地 Amber22 的 MMPBSA.py 不支持 "within 4" 语法(SelectionError),
        # 故用 "all" 全残基分解, 汇总时按界面距离(<4A, 从 _MMPBSA_complex.pdb 算)筛选。
        '  idecomp=1, print_res="all",\n'
        "/\n"
    )
    cmd = " ".join([
        MMPBSA, "-O", "-i", str(cin),
        "-o", str(d / "DECOMP_SUMMARY_500ns.dat"),
        "-do", str(dec),
        "-cp", str(d / "protein_only.prmtop"),
        "-rp", str(d / "receptor_mmpbsa.prmtop"),
        "-lp", str(d / "ligand_mmpbsa.prmtop"),
        "-y", str(d / "md_last20ns_nowat.nc"),
    ])
    bash_cmd = f"source {AMBERHOME}/amber.sh && {cmd}"
    with open(d / "_decomp.log", "w") as lf:
        p = subprocess.run(["bash", "-c", bash_cmd], stdout=lf,
                           stderr=subprocess.STDOUT, cwd=d)
    if dec.exists() and dec.stat().st_size > 1000:
        return name, "ok"
    return name, f"rc={p.returncode}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", action="store_true")
    ap.add_argument("--remote", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    jobs = jobs_from_finals(a.only)
    if a.status:
        for system, d in jobs:
            dec = d / "DECOMP_MMGBSA_500ns.dat"
            print(f"{'DONE' if dec.exists() else 'todo'}  {system}/{d.name}")
        print(f"total {len(jobs)}")
        return

    t0 = time.time()
    done = fail = 0
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_one, s, d): (s, d) for s, d in jobs}
        for f in as_completed(futs):
            s, d = futs[f]
            try:
                name, st = f.result()
            except Exception as e:
                name, st = d.name, f"EXC {e}"
            print(f"[{time.time()-t0:7.1f}s] {s}/{name}: {st}", flush=True)
            if st == "ok":
                done += 1
            elif st == "skip":
                pass
            else:
                fail += 1
    print(f"== done={done} fail={fail} skip={len(jobs)-done-fail} "
          f"({time.time()-t0:.0f}s) ==")


if __name__ == "__main__":
    main()
