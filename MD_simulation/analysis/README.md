# MD 分析流水线（SH3_VNAR / HCG_VNAR，48 体系 × 500–700 ns）

本目录包含 500 ns 延伸 MD 的完整分析与 MM/PBSA 结合自由能计算脚本（论文 Part 3 全部计算脚本），
以及汇总结果数据集。执行环境: AMBER22 (本地) / AmberTools23 (远程)；三机并行（本地 + 两台远程，ssh key 自备）。

## 脚本（按流水线顺序）

| 脚本 | 作用 |
|---|---|
| `14_batch_trend_analysis.py` | 每段末 20 ns MM/GBSA 趋势 + dG 演化收敛判定（`--local/--remote`） |
| `16_unified_rmsd.py` | 统一参考（AF3 起始姿态）重算 RMSD，含核心区掩码（SH3 `:1-183` / HCG `:1-392`） |
| `26_rmsd_reconverge.py` | 核心 RMSD 末 100 ns 窗收敛判定（std<2 Å 且 slope<1 Å/100ns） |
| `40_batch_pb_inp1.py` | **MM/PBSA 主计算**（PB `ipb=2/inp=1/indi=1.0/exdi=80.0/istrng=0.154/cavity_surften=0.0072`；末 20 ns 抽 100 帧；三机分流） |
| `41_batch_pb_decomp.py` | MM/PBSA + per-residue decomposition（`idecomp=1, print_res="all"`） |
| `39_collect_decomp.py` | decomp 明细汇总 → 热点残基表（本地编号→复合物编号换算内置） |

⚠️ 各脚本头部有 NOTE 说明需要替换的绝对路径/机器地址（源自内部流水线原样入库，保证与论文结果可对齐）。

**关于 14/38 号脚本的 GB 参数**：`14_batch_trend_analysis.py`（dG 段间演化收敛判定）与
`38_batch_decomp.py`（GB 版 decomp，已被 41 号 PB 版取代）在流水线历史上使用 MM/GBSA（igb=5）——
收敛判定只关心 ΔG 的**段间漂移**而非绝对值，GB 计算快 10–100 倍，适合趋势监测；
**论文报告的结合能数值全部来自 MM/PBSA（40/41 号脚本，inp=1）**。两者用途不同，不冲突。

## 汇总数据集（data/）

| 文件 | 内容 |
|---|---|
| `pb_inp1_vs_pb_inp2.csv` | 48 体系 MM/PBSA ΔG（inp=1 修复后 vs inp=2 旧值对照），含 SD/SEM |
| `decomp_pb_hotspots.csv` | 逐残基热点表（\|ΔG_res\| ≥ 1.0 kcal/mol，1,229 行；tag 已剔除，链已标注） |
| `decomp_pb_all_residues.csv` | 全残基分解明细（12,639 行） |

## 关键结论口径

- **非极性模型必须用 `inp=1`**：`inp=2`（Tan–Luo）的 EDISPER 项在本类 PPI 上约 +130 kcal/mol，
  曾导致 42/48 体系 ΔG 为正；修复后 48/48 全负（SH3 中位 −65.0，HCG 中位 −97.4）。
- 收敛体系才参与排序：结构判据（核心 RMSD 末 100 ns）+ 能量判据（段间 ΔG 漂移 < max(10%·\|ΔG\|, 1.0)）。
