# VNAR Optimization — 计算流程、软件参数与图表脚本

VNAR（鲨鱼单域抗体）亲和力成熟课题的计算全流程仓库：序列设计 → AF3 结构预测 → 聚类/静态分析 → 分子动力学模拟 → 结合自由能，以及论文全部图表的生成脚本。

## 目录结构

```
VNAR_Optimization/
├── README.md               ← 本文档：全流程软件参数与执行命令
├── AF3_input/              ← AF3 输入示例与流水线配置（含种子方案、阈值表）
├── MD_simulation/          ← AMBER MD 完整协议（mdin 文件 + 执行脚本）
├── 脚本/                   ← 所有 Python 绘图/分析脚本
├── 绘图结果/               ← 生成的 SVG/PNG/PDF 图表文件
└── Bamboo_MSA_18seq_Fisher.csv ← Fisher 分析依赖数据
```

## 计算流水线总览

```
IgGM 亲和力成熟（序列设计）
  → AntiBMPNN / SaProt 辅助打分
  → AlphaFold3 批量结构预测（突变体库）
  → 2STEP 两阶段聚类（Part 1）
  → PyRosetta 界面能量分析（Part 2）
  → AMBER MD 500 ns 动态模拟（Part 3）
  → MM/GBSA 结合自由能与残基分解
```

---

# 一、IgGM — 抗体序列与结构协同设计（亲和力成熟）

工具：[IgGM](https://github.com/tencent-ailab/IgGM)（Tencent AI Lab，ESM-650M + antibody design trunk + IGSO3）。用于对 2D4D2 VNAR–SH3 复合体做亲和力成熟采样。

## 实际运行命令（每 GPU 一个进程，4 GPU 并行）

```bash
cd /data/Tools/IgGM-master

export CUDA_VISIBLE_DEVICES="$GPU_ID"
python design.py \
    --fasta 2D4D2_maturation/2d4d2_sh3_design_corrected.fasta \
    --antigen 2D4D2_maturation/2d4d2_sh3.pdb \
    --output outputs/2D4D2_maturation_0116/gpu_X \
    --num_samples 10000 \
    --chunk_size 32 \
    --steps 10 \
    --temperature 1.0 \
    --max_antigen_size 400 \
    --run_task affinity_maturation \
    --fasta_origin 2D4D2_maturation/2d4d2_sh3_origin_corrected.fasta
```

## 关键参数

| 参数 | 值 | 说明 |
|------|-----|------|
| `--run_task` | `affinity_maturation` | 亲和力成熟模式（`--fasta_origin` 提供起始序列） |
| `--num_samples` | 10,000 / GPU | 总样本 40,000（4 × GPU 4-7） |
| `--steps` | 10 | 蒙特卡洛采样步数 |
| `--temperature` | 1.0 | 采样温度 |
| `--chunk_size` | 32 | 批处理样本数 |
| `--max_antigen_size` | 400 | 抗原残基数上限 |
| `--fasta` | 设计 FASTA（CDR3 区用 `X` 占位） | 抗体(VNAR)+抗原两段 FASTA，A 链 VNAR / B 链 SH3 |

输入 FASTA（A 链 = VNAR 2D4D2，尾部 `XXXXXXXXXXXXXXXX` 为 CDR3 设计区；B 链 = SH3 结构域序列）。

---

# 二、AntiBMPNN / SaProt — 辅助序列设计打分

两个模型作为**辅助先验**参与突变体优先级排序（非最终判据）：

| 工具 | 版本/权重 | 用途 |
|------|-----------|------|
| **AntiBMPNN** | `AntiBMPNN-main`（antibody-specific ProteinMPNN 变体，权重于 `/data/Tools/AntiBMPNN/antibmpnn_model_weights/`） | 抗体特异性 GNN，结构条件下的残基突变概率/ΔΔG 辅助打分 |
| **SaProt** | `westlake-repl/SaProt_650M_AF2`（HuggingFace） | 结构感知 token（AA + FoldSeek 3Di）的 1280 维残基嵌入，辅助稳定性打分 |

AntiBMPNN 基于 [ProteinMPNN 运行框架](https://github.com/dauparas/ProteinMPNN)（`Running_AntiBMPNN_run.py`，参数体系同 `protein_mpnn_run.py`：`--path_to_model_weights`、`--model_name`、`--temperature`、`--seed` 等）。

> **已知局限**（论文 Discussion）：AntiBMPNN 与 SaProt 均为残基级打分器，不编码跨链结合能项——对 E96F（SH3 体系）与 Q86S（HCG 体系）这两个实验致害突变给出假阳性推荐。案例的 Rosetta/MD 结构分析见 `MD_simulation/README.md` 与论文 Discussion 部分。

---

# 三、AlphaFold3 — 复合物结构批量预测

工具：AlphaFold 3（`alphafold3` v3.0.1，官方 `run_alphafold.py`）。输入格式、种子方案与筛选阈值详见 [`AF3_input/README.md`](AF3_input/README.md)。

## 实际运行命令

**初始复合物建模（含数据管线，MSA/templates 在线生成）：**

```bash
CUDA_VISIBLE_DEVICES=$GPU_ID python /data/Tools/AF3/alphafold3/run_alphafold.py \
    --json_path=input/hcg_vnar_trimer.json \
    --output_dir=output/hcg_vnar_trimer \
    --max_template_date=3000-12-01 \
    --run_data_pipeline=True \
    --run_inference=true
```

**突变体库批量预测（复用 WT MSA/templates，仅推理）：**

```bash
CUDA_VISIBLE_DEVICES=$GPU_ID python /data/Tools/AF3/alphafold3/run_alphafold.py \
    --json_path=$JOB_FILE \
    --output_dir="$OUTPUT_BASE/$JOB_NAME" \
    --run_data_pipeline=false \
    --run_inference=true
```

## 关键参数

| 参数 | 值 | 说明 |
|------|-----|------|
| `modelSeeds` | HCG: `[26226, 116740, 288390, 670488, 777573]`；SH3: `[87231, 49455, 37084, 89841, 63891]` | 每突变体 5 独立种子；体系内种子恒定，仅序列变化 |
| `--run_data_pipeline` | True（初始）/ false（突变体库） | 突变体输入替换 VNAR 链序列与 MSA 首序列（query），保留链 A/B 的 MSA |
| `--max_template_date` | 3000-12-01 | 禁用模板日期过滤 |
| 突变体库规模 | HCG trimer 170 变体；SH3 组合库 313 变体 | 见 `AF3_input/README.md` |
| `dialect` | `alphafold3` | 标准 AF3 JSON `sequences` 数组 |

MSA 替换脚本：`batch_fasta_to_json.py`（FASTA → AF3 JSON，同步更新 `unpairedMsa`/`pairedMsa` 的 query 首序列）。

**置信度过滤阈值**（`AF3_input/full_pipeline_example.yaml`）：

| 指标 | 阈值 |
|------|------|
| pLDDT | ≥ 0.7 |
| ipSAE | ≥ 0.6 |
| clashes | ≤ 5 |
| pDockQ | ≥ 0.2 |
| ipTM | ≥ 0.6 |

---

# 四、2STEP 两阶段聚类 — Part 1 结构分型

工具：`/data/Tools/IgGM-master/2STEP/`（自研脚本）。

## 第一步：粗聚类（结合模式分型）

```bash
python AF3_Cluster_Corse_v1.py
```

- 依据抗原接触集的 **Jaccard 距离**聚类（接触距离阈值 `contact_cutoff = 5.0 Å`）
- 配置：`CHAIN_CONFIG`（抗体链/抗原链、接触阈值）
- 输出：聚类结果 `.pkl`/`.csv`、可视化图、`coarse_clusters/` 结构文件夹

## 第二步：精细聚类（结构相似性）

```bash
python AF3_Cluster_fine_v1.py   # 配置见 config_fine_clustering.txt
```

- 基于 **Foldseek** 结构比对 + **US-align** TM-score 对粗聚类内结构再分型
- 关键配置：`COARSE_RESULTS_FILE`、`PDB_DIR`、`COARSE_CLUSTERS_DIR`、`N_JOBS=4`

---

# 五、PyRosetta — Part 2 界面能量静态分析

主脚本：`scripts/part2/part2_run_pyrosetta_static_relax_interface.py`（入口 `scripts/run_pyrosetta_static.py`）。

## 实际运行命令

```bash
# CSV 模式：候选结构 + 界面定义（ligand/receptor 列），先 FastRelax 再算界面能
python scripts/run_pyrosetta_static.py \
    --csv_path candidates.csv \
    --output_dir rosetta_static_out \
    --relax true --fixbb true --fixed_chain A \
    --batch_idx 1

# 目录模式：直接对聚类结构做静态界面分析
python scripts/run_pyrosetta_static.py \
    --pdb_dir fine_clusters/cluster_0 \
    --output_dir rosetta_static_out \
    --binder_chain B --target_chain A \
    --batch_idx 0
```

## 关键参数

| 参数 | 值 | 说明 |
|------|-----|------|
| `--relax` | `true` | FastRelax 松弛后评估（管线配置 `relax=true`） |
| `--fixbb` / `--fixed_chain` | `true` / 受体链 | 固定指定链骨架（PPIFlow 风格 repack-only） |
| `--dump_top_n` | 30 | 松弛后输出 top 30 结构（`*_relaxed.pdb`） |
| 打分器 | `InterfaceAnalyzerMover` | 输出 `interface_dG`、`interface_delta_sasa`、`complexed_sasa`（`rosetta_static_<batch>.csv`） |

Part 2 逻辑对齐 Germinal（`score_interface`）与 PPIFlow（`relax_complex.py`）的 InterfaceAnalyzerMover + FastRelax 流程。过滤阈值同 AF3 一节；通过后按 `interface_score` 排序送入 Part 3 MD。

---

# 六、AMBER MD — Part 3 动态模拟与 MM/GBSA

500 ns 分子动力学模拟（SH3/HCG 全突变体）与结合自由能计算的**完整协议、实际运行命令、全部 mdin 输入文件**见 [`MD_simulation/README.md`](MD_simulation/README.md)。要点：

- tleap 构系：ff14SB + TIP3P，八面体盒子 8 Å buffer，Na⁺/Cl⁻ 中和
- 预平衡：两阶段最小化（200 kcal·mol⁻¹·Å⁻² 约束 → 无约束）→ NVT 升温（2.0 约束）→ NPT 加压（2.0 约束）→ NPT 无约束
- 生产：100 ns NPT + 4×100 ns NVT 延伸（`06_run_md_extend.py` GPU 调度器）
- 分析：cpptraj 去水/RMSD → MM/GBSA（igb=5, saltcon=0.154）→ per-residue decomposition

---

---

# 七、绘图脚本 — Fisher 统计分析与论文图表

本目录还包含 VNAR 论文所有图表的生成脚本及输出结果。

### 一、环形图 (Donut Chart) — 数据集构成展示

| 脚本 | 功能 |
|------|------|
| `generate_donut_svg.py` | 生成 **物种分布环形图** 和 **可溶性标签环形图**（2张独立SVG），简单版。 |
| `generate_final_svgs.py` | **整合版**，一次性生成 3 张图：① 物种分布 Donut（Nature 配色，含拉丁学名）、② 可溶性标签 Donut、③ 3 集 Venn 图。SVG 文字保持可编辑（`svg.fonttype=none`）。 |

**输出文件：**
- `Fig_Chart_Species_Donut.svg` — 物种分布（Nurse/Bamboo/Wobbegong/Dogfish，N=257）
- `Fig_Chart_Labels_Donut.svg` — 可溶性标注（SOLUBLE 27 / INSOLUBLE 10 / Unlabeled 220）

### 二、韦恩图 (Venn Diagram) — Fisher 分析交叉验证

| 脚本 | 功能 |
|------|------|
| `generate_venn_svg.py` | 生成 **2 集 Venn 图**：All-species Fisher (9 位点) vs Species-controlled Fisher (6 位点)，标注 ROBUST 重叠位点（43G/83I/94T）。 |
| `calc_venn3.py` | **3 集 Venn 数据运算**：计算 L1a (All-species)、L1b (Species-controlled)、L2 (Bamboo-only) 三层 Fisher 分析的集合交集，输出各区数值。依赖 `MSA_Position_Mapping_ALL_to_Bamboo.csv` 位置映射表。 |
| `compute_venn_overlap.py` | **3 集重叠详细计算 + 位置映射**：将 L1a/L1b 的 ALL_MSA 编号映射到 Bamboo MSA 编号，计算与 L2 的交集，输出 ROBUST 位点详情及 Bamboo 特异性位点。 |
| `prepare_chart_data.py` | **图表数据准备**：归一化统计物种分布与标签类别，输出 `chart_data.json` 供图表脚本使用。 |

**输出文件：**
- `Fig_Chart_Venn_Fisher.svg` — 3 集 Venn（L1a ∩ L1b ∩ L2 三重叠分析）

### 三、序列 Logo (Sequence Logo) — 保守性可视化

| 脚本 | 功能 |
|------|------|
| `generate_sequence_logo_corrected.py` | **全物种 MSA Sequence Logo**（260 序列，205 列），基于 MSA 保守残基分析修正区域边界。9 区域 Partition（FR1-CDR1-FR2-HV2-FR3a-HV4-FR3b-CDR3-FR4），氨基酸按理化性质着色。 |
| `generate_sequence_logo_bamboo.py` | **Bamboo shark MSA Sequence Logo**（18 序列 SOL+INSOL，149 列），标记 14 个 Fisher 显著位点（橙色 ▼）和保守 Cys（金色 ★）。依赖 `Bamboo_MSA_18seq_Fisher.csv`。 |
| `script_sequence_logo_bamboo.py` | Bamboo Logo 另一版本（与上功能基本相同）。 |

**输出文件：**
- `VNAR_SequenceLogo_205col_Corrected.svg/.pdf` — 全物种 Logo（260 seqs）
- `VNAR_SequenceLogo_Bamboo_149col.svg/.pdf` — Bamboo Logo（18 seqs）
- `VNAR_SequenceLogo_Bamboo_22seq_FullMSA.svg/.pdf` — Bamboo Logo 全 MSA 版（22 seqs，含 ENGINEERING）

### 四、Fisher 统计分析与可视化

| 脚本 | 功能 |
|------|------|
| `bamboo_18seq_fisher.py` | **Bamboo shark Fisher 精确检验**：对 18 条（8 SOLUBLE + 10 INSOLUBLE）Bamboo MSA 进行逐列 Fisher 检验（单侧 greater），识别可溶性关联高频氨基酸位点，输出 `Bamboo_MSA_18seq_Fisher.csv`。并验证 3 个 ROBUST 位点在 Bamboo MSA 中的对应位置。 |
| `verify_bamboo_regions.py` | **Bamboo MSA 区域边界验证**：逐列分析 Bamboo MSA 氨基酸分布、熵值与保守性，验证 9 区域分区边界（特别是 Cys21/Trp35/Cys82/Gly112 四个保守锚点）。输出区域保守性统计表。 |

### 五、Fisher 层次分析图与综合图表

以下脚本位于 `figures/` 及 `high_freq_analysis/` 子目录中：

| 脚本 | 功能 |
|------|------|
| `chart1_manhattan_plot.py` | **曼哈顿图 (Fig2)**：Fisher -log10(p) 按位置排列，9 区域着色，标注显著性阈值线。 |
| `chart2_species_conformity.py` | **物种一致性分析图 (Fig3)**：堆叠条形图 + 分组条形图展示显著位点在各物种中的氨基酸保留率。 |
| `chart3_bamboo_shark_forest.py` | **森林图 (Fig5)**：Bamboo Fisher 显著位点的 OR 及 95%CI。 |
| `generate_visualization_plots.py` | **综合 6 图**：Volcano + Manhattan + Frequency + Multiple Testing + Forest + Table。 |
| `create_fisher_hierarchy_diagram.py` | **Fisher 层次分析逻辑推理图**：展示三层 Fisher 分析（L1a→L1b→L2）的递进逻辑关系。 |

**对应输出文件（Fig1–Fig6 系列，位于 `绘图结果/`）：**
- `Fig1_Volcano_Plot` — 火山图
- `Fig2_Manhattan_Plot` — 曼哈顿图
- `Fig3_Frequency_Comparison` — 频率对比图
- `Fig4_Multiple_Testing` — 多重检验校正图
- `Fig5_Forest_Plot` — 森林图
- `Fig6_Significant_Positions_Table` — 显著位点汇总表
- `Fig_Fisher_Hierarchy_Diagram.svg` — Fisher 层次分析逻辑图

---

## 依赖数据文件

| 文件 | 来源 | 用途 |
|------|------|------|
| `Bamboo_MSA_18seq_Fisher.csv` | `bamboo_18seq_fisher.py` 生成 | Logo 显著位点标记、统计分析 |
| `MSA_Position_Mapping_ALL_to_Bamboo.csv` | `high_freq_analysis/` | ALL_MSA → Bamboo MSA 列号映射 |
| `chart_data.json` | `prepare_chart_data.py` 生成 | Donut/Venn 图表数据 |

---

## 技术栈

- **绘图引擎**：Matplotlib（主）+ `matplotlib_venn`（Venn 图）
- **数据**：NumPy, Pandas, SciPy (`fisher_exact`)
- **图像格式**：SVG（矢量，可编辑）+ PDF（印刷）+ PNG（预览）
- **风格**：Nature 期刊配色（Navy/Slate/Ochre/Terracotta），Arial 字体，`svg.fonttype=none` 保证 Illustrator 可编辑

## 重新生成图表

```bash
# 生成 Donut/Venn 图（3 张）
python 脚本/generate_final_svgs.py

# 生成全物种 Sequence Logo
python 脚本/generate_sequence_logo_corrected.py

# 生成 Bamboo Sequence Logo（18 条）
python 脚本/generate_sequence_logo_bamboo.py

# 运行 Bamboo Fisher 检验（生成 CSV）
python 脚本/bamboo_18seq_fisher.py
```
