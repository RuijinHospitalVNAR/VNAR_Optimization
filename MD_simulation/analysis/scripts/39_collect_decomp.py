#!/usr/bin/env python3
# NOTE: 复现时把 EA 根路径 /data/wcf/MD_trend_analysis 替换为你的数据根
# (含 extend_analysis/{SYSTEM}/{name}/)。
"""Collect per-residue decomposition (MM/PBSA, inp=1) -> CSV summary.

输入: extend_analysis/{SYSTEM}/{name}/DECOMP_MMGBSA_PB_INP1.dat (41 脚本产物)
输出: extend_analysis/decomp_pb_all_residues.csv / decomp_pb_hotspots.csv
DELTAs 段 R/L 行各用 receptor/ligand 本地编号; 换算复合物编号:
SH3 ligand+105 (>=184 尾); HCG ligand+281 (>392 tag)。
"""
from __future__ import annotations
import csv
import re
from pathlib import Path

EA = Path("/data/wcf/MD_trend_analysis/extend_analysis")
HOT_THR = 1.0

ROW_RE = re.compile(
    r"^([A-Z]{3})\s*(\d+),([RL]) [A-Z]{3}\s*(\d+),"
    + r"".join(r"(-?[\d.]+(?:[eE][-+]?\d+)?)," for _ in range(17))
    + r"(-?[\d.]+(?:[eE][-+]?\d)?)$")


def parse_decomp(path: Path) -> list[dict]:
    rows = []
    in_delta = False
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("DELTAS:"):
            in_delta = True
            continue
        if not in_delta:
            continue
        m = ROW_RE.match(line.strip())
        if not m:
            if line.strip() and not line.startswith(",") and rows:
                break  # Sidechain/Backbone 段
            continue
        g = m.groups()
        rows.append({
            "resname": g[0], "resid": int(g[3]), "part": g[2],
            "vdw": float(g[7]), "eel": float(g[10]),
            "egb": float(g[13]),   # PB decomp 中此列 = 极性溶剂化 (EPB)
            "enpol": float(g[16]),
            "total": float(g[19]),
            "sd": float(g[20]), "sem": float(g[21]),
        })
    return rows


def to_complex(system: str, r: dict) -> tuple[int, str]:
    if system == "SH3_VNAR":
        if r["part"] == "R":
            return r["resid"], "VNAR"
        cid = 105 + r["resid"]
        return cid, ("antigen" if cid <= 183 else "antigen-tail")
    if r["part"] == "R":
        return r["resid"], "antigen"
    cid = 281 + r["resid"]
    return cid, ("VNAR" if cid <= 392 else "tag")


def main():
    all_rows, hot_rows = [], []
    nsys = 0
    for system in ("SH3_VNAR", "HCG_VNAR"):
        for d in sorted(EA.joinpath(system).iterdir()):
            dec = d / "DECOMP_MMGBSA_PB_INP1.dat"
            if not dec.exists():
                continue
            rows = parse_decomp(dec)
            if not rows:
                print(f"WARN empty: {system}/{d.name}")
                continue
            nsys += 1
            for r in rows:
                cid, chain = to_complex(system, r)
                rec = {"system": system, "name": d.name, "chain": chain, "resid": cid,
                       **{k: v for k, v in r.items() if k not in ("resid", "part")}}
                all_rows.append(rec)
                if abs(r["total"]) >= HOT_THR and chain != "tag":
                    hot_rows.append(rec)
    fields = ["system", "name", "chain", "resid", "resname",
              "vdw", "eel", "egb", "enpol", "total", "sd", "sem"]
    with open(EA / "decomp_pb_all_residues.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        w.writerows(all_rows)
    hot_rows.sort(key=lambda r: (r["system"], r["name"], r["total"]))
    with open(EA / "decomp_pb_hotspots.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        w.writerows(hot_rows)
    print(f"systems={nsys}  all={len(all_rows)}  hot(|dG|>={HOT_THR})={len(hot_rows)}")


if __name__ == "__main__":
    main()
