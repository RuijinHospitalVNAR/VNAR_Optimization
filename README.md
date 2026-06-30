# VNAR Production Protocol — 绘图脚本说明

本目录包含 VNAR 论文所有图表的生成脚本及输出结果。

## 目录结构

```
绘图脚本/
├── README.md
├── 脚本/           ← 所有 Python 绘图/分析脚本
├── 绘图结果/        ← 生成的 SVG/PNG/PDF 图表文件
└── Bamboo_MSA_18seq_Fisher.csv  ← 依赖数据文件
```

---

## 脚本说明

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
