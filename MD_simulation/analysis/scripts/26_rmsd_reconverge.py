#!/usr/bin/env python3
# NOTE: 本脚本源自项目内部流水线, 含两台计算机的绝对路径。
# 复现时请全局替换:
#   /data/wcf/MD_trend_analysis  -> 本地仓库根 (analysis 数据目录的父目录)
#   /data/md_extend              -> 远程机工作根
#   /home/rcdb/anaconda3/envs/AmberTools23 -> 远程 AMBERHOME (AmberTools23)
#   /data/Tools/Amber22          -> 本地 AMBERHOME (Amber22)
#   rcdb@100.76.253.60 / rcdb@100.97.18.74:10021 -> 计算机地址(需自备 ssh key)
"""重算 RMSD 收敛统计，基于 unified_rmsd.dat（完整 0~500ns，PBC autoimage 修复后）。

背景：14 脚本的 rmsd 字段基于 rmsd_bb.dat（旧，只覆盖延伸段、参考延伸起点、9/4），
而 15 脚本的 RMSD 曲线读 unified_rmsd.dat（16 脚本 9/10 重算，覆盖 0~500ns、参考原始第 1 帧）。
二者口径不一致 → 图上曲线与收敛标签对不上。本脚本只重算 rmsd 字段，保留
dG_evolution / nseg / total_ext_ns 等其余字段不变。

收敛判据（与 14 脚本 208-209 行一致）：
  converged_v2 = 末 100ns 时间窗 std < 2.0 Å 且 |线性斜率| < 1.0 Å/100ns
  converged(旧) = 末两窗 mean 漂移 < thr 且 末窗斜率 < thr，thr = max(0.05*m1, 0.1)
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path("/data/wcf/MD_trend_analysis/extend_analysis")
SUMMARY_FILES = [ROOT / "trend_summary.json",
                 ROOT / "trend_summary_remote.json"]


def _stats(p: Path) -> dict | None:
    """对单个 unified_rmsd*.dat 计算全部收敛统计"""
    if not p.exists():
        return None
    a = np.loadtxt(p, comments=("#", "@"))
    t, v = a[:, 0], a[:, 1]  # t: ns, v: Å
    # 以数据实际终点为准(而非 100+nseg*100): 跨服务器迁移体系实际长度非整百
    # (S85T 609 / S86R_Y92D 677 / S85I_G91D_H95Y 769 ns), 用整百假设会使"末100ns"窗错位
    t_end = float(t[-1])
    nseg = max(int(round((t_end - 100.0) / 100.0)), 1)

    # 完整 0~t_end 分窗，每 100ns 一窗（原始 100ns + nseg 段），末窗含到轨迹尾
    wins = []
    for w in range(nseg + 1):
        lo = w * 100
        hi = (w + 1) * 100
        m = (t >= lo) if w == nseg else ((t >= lo) & (t < hi))
        if m.sum() < 10:
            wins.append(None)
        else:
            wins.append({"window": f"{lo}-{round(t_end) if w == nseg else hi}ns",
                         "mean": float(v[m].mean()),
                         "std": float(v[m].std()),
                         "n": int(m.sum())})

    m_last = t >= (t_end - 100)
    mean_last100 = float(v[m_last].mean()) if m_last.sum() else float("nan")

    # 旧判据：末两窗 mean 漂移 + 末窗线性斜率
    if None not in wins and len(wins) >= 2:
        m1, m2 = wins[-2]["mean"], wins[-1]["mean"]
        drift = abs(m2 - m1)
        mw = (t >= (t_end - 100)) & (t <= t_end)
        slope_win = float(np.polyfit(t[mw], v[mw], 1)[0] * 100.0) if mw.sum() >= 2 else 0.0
        thr = max(0.05 * abs(m1), 0.1)
        converged = bool(drift < thr and abs(slope_win) < thr)
        drift_last2win = float(drift)
    else:
        converged, drift_last2win = None, None

    # 新判据：末 100ns 时间窗 std < 2.0 Å 且 |斜率| < 1.0 Å/100ns
    m_tail = t >= (t_end - 100)
    tv, vv = t[m_tail], v[m_tail]
    if len(vv) >= 10:
        lin_slope = float(np.polyfit(tv, vv, 1)[0] * 100.0)
        tail_std = float(vv.std())
        conv2 = bool(tail_std < 2.0 and abs(lin_slope) < 1.0)
    else:
        lin_slope = tail_std = 0.0
        conv2 = False

    # 末 50ns 窗（辅助闸: 识别"晚期跃迁落在末100ns窗内")
    m50 = t >= (t_end - 50)
    if m50.sum() >= 10:
        s50 = float(v[m50].std())
        l50 = float(np.polyfit(t[m50], v[m50], 1)[0] * 100.0)
        conv50 = bool(s50 < 2.0 and abs(l50) < 1.0)
    else:
        s50 = l50 = 0.0
        conv50 = False

    return {
        "mean_last100ns": round(mean_last100, 6),
        "max": float(v.max()),
        "windows": wins,
        "converged": converged,
        "drift_last2win": drift_last2win,
        "tail_std": round(tail_std, 6),
        "tail_slope_per100ns": round(lin_slope, 6),
        "converged_v2": conv2,
        "tail50_std": round(s50, 6),
        "tail50_slope_per100ns": round(l50, 6),
        "converged_50ns": conv50,
    }


def recompute_rmsd(system: str, name: str, nseg: int) -> tuple[dict, dict] | None:
    """收敛判据用**核心区**(排除固有无序段/纯化标签); 另附全复合物(监控用)。

    核心区掩码: SH3 :1-183 (去抗原脯氨酸尾184-211) / HCG :1-392 (去 VNAR tag 393-414)。
    见 handoff 十七节: 全复合物掩码把抗原无序段的摆动误报为"复合物未收敛"。
    """
    d = ROOT / system / name
    core = _stats(d / "unified_rmsd_core.dat")
    whole = _stats(d / "unified_rmsd.dat")
    if core is None and whole is None:
        return None
    # 回退: 无核心区文件时用全复合物(并标注), 避免静默丢数据
    if core is None:
        core = dict(whole); core["fallback_whole"] = True
    if whole is not None:
        whole = {k: whole[k] for k in
                 ("mean_last100ns", "tail_std", "tail_slope_per100ns", "converged_v2",
                  "max")}
    return core, whole


def main():
    for f in SUMMARY_FILES:
        if not f.exists():
            continue
        rs = json.load(open(f))
        n_upd = n_miss = n_fb = 0
        for r in rs:
            got = recompute_rmsd(r["system"], r["name"], r["nseg"])
            if got is None:
                n_miss += 1
                continue
            core, whole = got
            if core.pop("fallback_whole", False):
                n_fb += 1
            r["rmsd"] = core                     # 判据字段 = 核心区
            r["rmsd_whole"] = whole              # 监控字段 = 全复合物
            n_upd += 1
        json.dump(rs, open(f, "w"), ensure_ascii=False, indent=2)
        print(f"{f.name}: 更新 {n_upd}, 缺文件 {n_miss}, 回退全复合物 {n_fb}")


if __name__ == "__main__":
    main()
