#!/usr/bin/env python3
# NOTE: 本脚本源自项目内部流水线, 含两台计算机的绝对路径。
# 复现时请全局替换:
#   /data/wcf/MD_trend_analysis  -> 本地仓库根 (analysis 数据目录的父目录)
#   /data/md_extend              -> 远程机工作根
#   /home/rcdb/anaconda3/envs/AmberTools23 -> 远程 AMBERHOME (AmberTools23)
#   /data/Tools/Amber22          -> 本地 AMBERHOME (Amber22)
#   rcdb@100.76.253.60 / rcdb@100.97.18.74:10021 -> 计算机地址(需自备 ssh key)
"""Step 40: PB (inp=1, 经典 SASA 非极性) 全量重算 — 替换 inp=2 的失效结果。

根因 (2026-09-20 诊断): 旧 mmpbsa_pb.in 用 inp=2 (Tan-Luo 非极性分解),
EDISPER 色散回填中位 +132.5 kcal/mol 淹没全部有利项 → 42/48 体系 ΔG_PB 为正。
修复: inp=1 + cavity_surften=0.0072 (PPI 文献标准 SASA 配置), 本地 10 帧验证
S85I_G91D +57.0 → -44.9。PB 计算量约为 GB 的 10-100 倍 (有限差分网格求解)。

口径与 500ns 统一重算一致: md_last20ns_nowat.nc (seg4 末 20ns 2000 帧),
interval=20 → 100 帧。掩码自适应 (hist_masks, HCG 含 tag :282-414)。
输出: {outdir}/FINAL_MMPBSA_PB_INP1.dat (+ PB_INP1_PER_FRAME.csv)

用法: python3 40_batch_pb_inp1.py --local|--remote1|--remote2 [--jobs N]
      [--only SUBSTR] [--status] [--exclude SUBSTR]
remote1 = 100.76.253.60 /data/md_extend/extend_analysis (AmberTools23)
remote2 = 100.97.18.74 /data/pb_inp1 (at23 conda env, 输入由 remote1 rsync 过去)
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MODE = next((m for m in ("--local", "--remote1", "--remote2") if m in sys.argv), "--local")
if MODE == "--local":
    OUT_ROOT = Path("/data/wcf/MD_trend_analysis/extend_analysis")
    MMPBSA = "/data/Tools/Amber22/bin/MMPBSA.py"
    AMBERHOME = "/data/Tools/Amber22"
elif MODE == "--remote1":
    OUT_ROOT = Path("/data/md_extend/extend_analysis")
    MMPBSA = "/home/rcdb/anaconda3/envs/AmberTools23/bin/MMPBSA.py"
    AMBERHOME = "/home/rcdb/anaconda3/envs/AmberTools23"
else:
    OUT_ROOT = Path("/data/pb_inp1/extend_analysis")
    MMPBSA = "/home/rcdb/anaconda3/envs/at23/bin/MMPBSA.py"
    AMBERHOME = "/home/rcdb/anaconda3/envs/at23"

MASKS = {
    "SH3_VNAR": {"receptor": ":1-105", "ligand": ":106-211"},
    "HCG_VNAR": {"receptor": ":1-281", "ligand": ":282-414"},   # HCG ligand 含 tag!
}
OUT_DAT = "FINAL_MMPBSA_PB_INP1.dat"
OUT_EO = "PB_INP1_PER_FRAME.csv"


def hist_masks(d: Path) -> tuple[str, str] | None:
    """从该体系历史 MMPBSA 输入提取掩码 (对该体系拓扑必然有效)。"""
    for f in ("mmpbsa_pb.in", "mmpbsa_gb_extend.in", "_mmpbsa_seg4.in", "_mmpbsa_seg1.in", "mmpbsa_gb.in"):
        p = d / f
        if not p.exists():
            continue
        t = p.read_text(errors="ignore")
        mr = re.search(r"receptor_mask='([^']+)'", t)
        ml = re.search(r"ligand_mask='([^']+)'", t)
        if mr and ml:
            return mr.group(1), ml.group(1)
    return None


def jobs_all(only: str | None, exclude: str | None) -> list[tuple[str, Path]]:
    out = []
    for system in MASKS:
        if not (OUT_ROOT / system).is_dir():
            continue
        for d in sorted((OUT_ROOT / system).iterdir()):
            if not d.is_dir():
                continue
            finals = [d / "FINAL_RESULTS_MMGBSA_BINDING_500ns.dat", d / "FINAL_MMGBSA_500ns.dat"]
            if not any(f.exists() for f in finals):
                continue
            if only and only not in d.name:
                continue
            if exclude and exclude in d.name:
                continue
            pre = ["md_last20ns_nowat.nc", "protein_only.prmtop",
                   "receptor_mmpbsa.prmtop", "ligand_mmpbsa.prmtop"]
            if any(not (d / f).exists() for f in pre):
                continue
            out.append((system, d))
    return out


def run_one(system: str, d: Path) -> tuple[str, str]:
    name = d.name
    if (d / OUT_DAT).exists():
        return name, "skip"
    masks = hist_masks(d) or (MASKS[system]["receptor"], MASKS[system]["ligand"])
    cin = d / "_pb_inp1.in"
    cin.write_text(
        "&general\n"
        "  startframe=1, endframe=2000, interval=20,\n"
        "  verbose=1, entropy=0,\n"
        f"  receptor_mask='{masks[0]}', ligand_mask='{masks[1]}',\n"
        "/\n"
        "&pb\n"
        "  ipb=2, inp=1, indi=1.0, exdi=80.0, istrng=0.154,\n"
        "  radiopt=0, cavity_surften=0.0072, cavity_offset=0.0,\n"
        "/\n"
    )
    cmd = " ".join([
        MMPBSA, "-O", "-i", str(cin),
        "-o", str(d / OUT_DAT), "-eo", str(d / OUT_EO),
        "-cp", str(d / "protein_only.prmtop"),
        "-rp", str(d / "receptor_mmpbsa.prmtop"),
        "-lp", str(d / "ligand_mmpbsa.prmtop"),
        "-y", str(d / "md_last20ns_nowat.nc"),
    ])
    bash_cmd = f"export AMBERHOME={AMBERHOME} && {cmd}"
    with open(d / "_pb_inp1.log", "w") as lf:
        p = subprocess.run(["bash", "-c", bash_cmd], stdout=lf,
                           stderr=subprocess.STDOUT, cwd=d)
    if (d / OUT_DAT).exists() and (d / OUT_DAT).stat().st_size > 500:
        return name, "ok"
    return name, f"rc={p.returncode}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", action="store_true")
    ap.add_argument("--remote1", action="store_true")
    ap.add_argument("--remote2", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--exclude", type=str, default=None)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    jobs = jobs_all(a.only, a.exclude)
    if a.status:
        n_done = 0
        for system, d in jobs:
            ok = (d / OUT_DAT).exists()
            n_done += ok
            print(f"{'DONE' if ok else 'todo'}  {system}/{d.name}")
        print(f"total {len(jobs)} done {n_done}")
        return

    t0 = time.time()
    done = fail = skip = 0
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_one, s, d): (s, d) for s, d in jobs}
        for f in as_completed(futs):
            s, d = futs[f]
            try:
                name, st = f.result()
            except Exception as e:
                name, st = d.name, f"EXC {e}"
            print(f"[{time.time()-t0:8.1f}s] {s}/{name}: {st}", flush=True)
            done += st == "ok"; fail += st not in ("ok", "skip"); skip += st == "skip"
    print(f"== done={done} fail={fail} skip={skip} ({time.time()-t0:.0f}s) ==")


if __name__ == "__main__":
    main()
