# VNAR Optimization — Reproducible computational workflow for VNAR affinity maturation

This repository hosts the complete computational pipeline used for the affinity maturation of VNAR (shark single-domain antibody) candidates reported in the accompanying manuscript. It provides the software parameters and run commands for every stage — sequence design → AlphaFold3 (AF3) structure prediction → clustering/static analysis → molecular dynamics (MD) simulation → binding free-energy calculation — together with the scripts that generated all figures in the paper.

A parallel Chinese version of this document is available in [`README.md`](README.md).

## Repository layout

```
VNAR_Optimization/
├── README.md               ← Chinese summary of the full pipeline and parameters
├── README_EN.md            ← This document (English)
├── AF3_input/              ← Representative AF3 input JSON files and pipeline configuration (seed scheme, thresholds)
├── MD_simulation/          ← Complete AMBER MD protocol (mdin input files + run scripts)
├── 脚本/                   ← All Python plotting/analysis scripts
├── 绘图结果/               ← Generated figure files (SVG/PNG/PDF)
└── Bamboo_MSA_18seq_Fisher.csv ← Input data used by the Fisher analysis
```

## Computational pipeline at a glance

```
IgGM affinity maturation (sequence design)
  → AntiBMPNN / SaProt auxiliary scoring
  → AlphaFold3 batch structure prediction (mutant libraries)
  → 2STEP two-stage clustering (Part 1)
  → PyRosetta interface-energy analysis (Part 2)
  → AMBER MD 500 ns dynamic simulation (Part 3)
  → MM/GBSA binding free energy and per-residue decomposition
```

---

# 1. IgGM — antibody sequence/structure co-design (affinity maturation)

Tool: [IgGM](https://github.com/tencent-ailab/IgGM) (Tencent AI Lab; ESM-650M + antibody design trunk + IGSO3). It was used to sample affinity-matured sequences for the 2D4D2 VNAR–SH3 complex.

## Run command (one process per GPU; 4 GPUs in parallel)

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

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--run_task` | `affinity_maturation` | Affinity-maturation mode (`--fasta_origin` supplies the starting sequence) |
| `--num_samples` | 10,000 / GPU | Total 40,000 samples (4 × GPU 4–7) |
| `--steps` | 10 | Number of Monte Carlo sampling steps |
| `--temperature` | 1.0 | Sampling temperature |
| `--chunk_size` | 32 | Batch size of samples |
| `--max_antigen_size` | 400 | Upper limit of antigen residue count |
| `--fasta` | Design FASTA (CDR3 region occupied by `X`) | Two-chain FASTA (VNAR + antigen); chain A = VNAR, chain B = SH3 |

Input FASTA: chain A = VNAR 2D4D2 with the trailing `XXXXXXXXXXXXXXXX` marking the CDR3 design region; chain B = the SH3 domain sequence.

---

# 2. AntiBMPNN / SaProt — auxiliary sequence-design scoring

Both models were used as **auxiliary priors** for ranking candidate mutations (not as the final criterion):

| Tool | Version / weights | Purpose |
|------|-------------------|---------|
| **AntiBMPNN** | `AntiBMPNN-main` (antibody-specific ProteinMPNN variant; weights in `/data/Tools/AntiBMPNN/antibmpnn_model_weights/`) | Antibody-specific graph neural network; residue mutation probability / ΔΔG auxiliar score conditioned on structure |
| **SaProt** | `westlake-repl/SaProt_650M_AF2` (Hugging Face) | Structure-aware token (AA + FoldSeek 3Di) 1280-dimensional residue embedding; auxiliary stability score |

AntiBMPNN follows the [ProteinMPNN framework](https://github.com/dauparas/ProteinMPNN) (run via `Running_AntiBMPNN_run.py`; parameters mirror `protein_mpnn_run.py`: `--path_to_model_weights`, `--model_name`, `--temperature`, `--seed`, etc.).

> **Known limitation** (paper Discussion): AntiBMPNN and SaProt are residue-level scorers that do not encode cross-chain binding energetics — they produced false-positive recommendations for the experimentally deleterious mutations E96F (SH3 system) and Q86S (HCG system). Rosetta/MD structural analyses of these cases are described in `MD_simulation/README.md` and the paper Discussion.

---

# 3. AlphaFold3 — batch complex-structure prediction

Tool: AlphaFold 3 (`alphafold3` v3.0.1, official `run_alphafold.py`). Input format, seed scheme, and screening thresholds are detailed in [`AF3_input/README.md`](AF3_input/README.md).

## Run commands

**Initial complex modelling (with data pipeline; MSA/templates generated online):**

```bash
CUDA_VISIBLE_DEVICES=$GPU_ID python /data/Tools/AF3/alphafold3/run_alphafold.py \
    --json_path=input/hcg_vnar_trimer.json \
    --output_dir=output/hcg_vnar_trimer \
    --max_template_date=3000-12-01 \
    --run_data_pipeline=True \
    --run_inference=true
```

**Batch prediction of mutant libraries (reuse WT MSA/templates, inference only):**

```bash
CUDA_VISIBLE_DEVICES=$GPU_ID python /data/Tools/AF3/alphafold3/run_alphafold.py \
    --json_path=$JOB_FILE \
    --output_dir="$OUTPUT_BASE/$JOB_NAME" \
    --run_data_pipeline=false \
    --run_inference=true
```

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `modelSeeds` | HCG: `[26226, 116740, 288390, 670488, 777573]`; SH3: `[87231, 49455, 37084, 89841, 63891]` | 5 independent seeds per variant; seeds are fixed within a system and only the sequence varies |
| `--run_data_pipeline` | True (initial) / false (mutant libraries) | For mutants, the VNAR chain sequence and the first (query) MSA sequence are replaced while the MSA of chains A/B is retained |
| `--max_template_date` | 3000-12-01 | Disables template date filtering |
| Mutant-library sizes | HCG trimer 170 variants; SH3 combination library 313 variants | See `AF3_input/README.md` |
| `dialect` | `alphafold3` | Standard AF3 JSON `sequences` array |

MSA-replacement script: `batch_fasta_to_json.py` (FASTA → AF3 JSON, updating the query first sequence of `unpairedMsa`/`pairedMsa`).

**Confidence-filter thresholds** (`AF3_input/full_pipeline_example.yaml`):

| Metric | Threshold |
|--------|-----------|
| pLDDT | ≥ 0.7 |
| ipSAE | ≥ 0.6 |
| clashes | ≤ 5 |
| pDockQ | ≥ 0.2 |
| ipTM | ≥ 0.6 |

---

# 4. 2STEP two-stage clustering — Part 1 structure grouping

Tool: `/data/Tools/IgGM-master/2STEP/` (in-house scripts).

## Step 1: coarse clustering (binding-mode grouping)

```bash
python AF3_Cluster_Corse_v1.py
```

- Clusters by **Jaccard distance** of antigen contact sets (contact-distance cutoff `contact_cutoff = 5.0 Å`)
- Configuration: `CHAIN_CONFIG` (antibody chain / antigen chain, contact threshold)
- Outputs: clustering results `.pkl`/`.csv`, visualization plots, `coarse_clusters/` structure directories

## Step 2: fine clustering (structural similarity)

```bash
python AF3_Cluster_fine_v1.py   # configuration in config_fine_clustering.txt
```

- Re-groups structures within each coarse cluster using **Foldseek** structural alignment + **US-align** TM-score
- Key settings: `COARSE_RESULTS_FILE`, `PDB_DIR`, `COARSE_CLUSTERS_DIR`, `N_JOBS=4`

---

# 5. PyRosetta — Part 2 static interface-energy analysis

Main script: `scripts/part2/part2_run_pyrosetta_static_relax_interface.py` (entry point `scripts/run_pyrosetta_static.py`).

## Run commands

```bash
# CSV mode: candidate structures + interface definition (ligand/receptor columns); FastRelax then interface energy
python scripts/run_pyrosetta_static.py \
    --csv_path candidates.csv \
    --output_dir rosetta_static_out \
    --relax true --fixbb true --fixed_chain A \
    --batch_idx 1

# Directory mode: static interface analysis directly on clustered structures
python scripts/run_pyrosetta_static.py \
    --pdb_dir fine_clusters/cluster_0 \
    --output_dir rosetta_static_out \
    --binder_chain B --target_chain A \
    --batch_idx 0
```

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--relax` | `true` | Evaluate after FastRelax relaxation (pipeline config `relax=true`) |
| `--fixbb` / `--fixed_chain` | `true` / receptor chain | Fix the backbone of the specified chain (PPIFlow-style repack-only) |
| `--dump_top_n` | 30 | Write top-30 relaxed structures (`*_relaxed.pdb`) |
| Scorer | `InterfaceAnalyzerMover` | Outputs `interface_dG`, `interface_delta_sasa`, `complexed_sasa` (`rosetta_static_<batch>.csv`) |

Part 2 mirrors the Germinal (`score_interface`) and PPIFlow (`relax_complex.py`) `InterfaceAnalyzerMover` + FastRelax workflows. Filter thresholds are the same as in the AF3 section; surviving variants are sorted by `interface_score` and advanced to Part 3 (MD).

---

# 6. AMBER MD — Part 3 dynamic simulation and MM/GBSA

The **complete protocol, run commands, and all `mdin` input files** for the 500 ns all-atom MD simulations (all SH3/HCG variants) and the binding free-energy calculations are provided in [`MD_simulation/README.md`](MD_simulation/README.md). Key points:

- tleap construction: ff14SB + TIP3P, octahedral box with an 8 Å buffer, neutralized with Na⁺/Cl⁻
- Pre-equilibration: two-stage minimization (200 kcal·mol⁻¹·Å⁻² restraint → unrestrained) → NVT heating (2.0 restraint) → NPT pressurization (2.0 restraint) → unrestrained NPT
- Production: 100 ns NPT + 4 × 100 ns NVT extensions (`06_run_md_extend.py` GPU scheduler)
- Analysis: cpptraj water removal/RMSD → MM/GBSA (`igb=5`, `saltcon=0.154`) → per-residue decomposition

---

# 7. Plotting scripts — Fisher statistical analysis and manuscript figures

All figure-generation scripts and their outputs for the VNAR paper are also included.

### a) Donut charts — dataset composition

| Script | Function |
|--------|----------|
| `generate_donut_svg.py` | Generates the **species-distribution** and **solubility-label** donut charts as 2 separate SVGs (simple version). |
| `generate_final_svgs.py` | **Integrated version** producing 3 figures in one run: ① species-distribution Donut (Nature palette, including Latin names), ② solubility-label Donut, ③ 3-set Venn diagram. SVG text is kept editable (`svg.fonttype=none`). |

**Outputs:** `Fig_Chart_Species_Donut.svg` (Nurse/Bamboo/Wobbegong/Dogfish, N=257), `Fig_Chart_Labels_Donut.svg` (SOLUBLE 27 / INSOLUBLE 10 / Unlabeled 220).

### b) Venn diagrams — Fisher cross-validation

| Script | Function |
|--------|----------|
| `generate_venn_svg.py` | 2-set Venn: All-species Fisher (9 positions) vs Species-controlled Fisher (6 positions), marking the ROBUST overlap positions (43G/83I/94T). |
| `calc_venn3.py` | 3-set Venn data: intersection of L1a (All-species), L1b (Species-controlled), and L2 (Bamboo-only) Fisher analyses. Depends on the `MSA_Position_Mapping_ALL_to_Bamboo.csv` mapping. |
| `compute_venn_overlap.py` | Detailed 3-set overlap + position mapping: maps ALL_MSA numbers of L1a/L1b to Bamboo MSA numbers, computes intersections with L2, and outputs ROBUST and Bamboo-specific positions. |
| `prepare_chart_data.py` | Data preparation: normalized species/solubility counts, outputs `chart_data.json`. |

**Output:** `Fig_Chart_Venn_Fisher.svg` (3-set Venn, L1a ∩ L1b ∩ L2).

### c) Sequence logos — conservation visualization

| Script | Function |
|--------|----------|
| `generate_sequence_logo_corrected.py` | All-species MSA sequence logo (260 sequences, 205 columns) with corrected region boundaries; 9 regions (FR1–CDR1–FR2–HV2–FR3a–HV4–FR3b–CDR3–FR4), amino acids colored by physicochemical property. |
| `generate_sequence_logo_bamboo.py` | Bamboo-shark MSA logo (18 sequences SOL+INSOL, 149 columns) marking 14 Fisher-significant positions (orange ▼) and conserved Cys (gold ★). Depends on `Bamboo_MSA_18seq_Fisher.csv`. |
| `script_sequence_logo_bamboo.py` | Variant of the Bamboo logo script (functionally similar). |

**Outputs:** `VNAR_SequenceLogo_205col_Corrected.svg/.pdf` (260 seqs), `VNAR_SequenceLogo_Bamboo_149col.svg/.pdf` (18 seqs), `VNAR_SequenceLogo_Bamboo_22seq_FullMSA.svg/.pdf` (22 seqs, including ENGINEERING).

### d) Fisher statistical analysis and visualization

| Script | Function |
|--------|----------|
| `bamboo_18seq_fisher.py` | Bamboo-shark Fisher exact test: column-wise Fisher test (one-sided, greater) over 18 Bamboo MSA sequences (8 SOLUBLE + 10 INSOLUBLE) to identify solubility-associated high-frequency positions; outputs `Bamboo_MSA_18seq_Fisher.csv`. Also verifies the 3 ROBUST positions in the Bamboo MSA. |
| `verify_bamboo_regions.py` | Bamboo MSA region-boundary validation: per-column amino-acid distribution, entropy, and conservation; verifies the 9-region partition boundaries (especially the four conserved anchors Cys21/Trp35/Cys82/Gly112). Outputs a region-conservation table. |

### e) Fisher hierarchical analyses and composite figures

Scripts under `figures/` and `high_freq_analysis/`:

| Script | Function |
|--------|----------|
| `chart1_manhattan_plot.py` | **Manhattan plot (Fig2)**: Fisher −log10(p) by position, colored by the 9 regions, with significance-threshold lines. |
| `chart2_species_conformity.py` | **Species-conformity analysis (Fig3)**: stacked + grouped bar charts showing residue retention of significant positions across species. |
| `chart3_bamboo_shark_forest.py` | **Forest plot (Fig5)**: OR and 95% CI of Bamboo Fisher-significant positions. |
| `generate_visualization_plots.py` | **Combined 6-figure generator**: Volcano + Manhattan + Frequency + Multiple Testing + Forest + Table. |
| `create_fisher_hierarchy_diagram.py` | **Fisher hierarchy logic diagram**: the L1a→L1b→L2 three-layer Fisher logic. |

**Corresponding outputs (Fig1–Fig6 series, in `绘图结果/`):** `Fig1_Volcano_Plot`, `Fig2_Manhattan_Plot`, `Fig3_Frequency_Comparison`, `Fig4_Multiple_Testing`, `Fig5_Forest_Plot`, `Fig6_Significant_Positions_Table`, `Fig_Fisher_Hierarchy_Diagram.svg`.

---

## Dependent data files

| File | Source | Purpose |
|------|--------|---------|
| `Bamboo_MSA_18seq_Fisher.csv` | generated by `bamboo_18seq_fisher.py` | Logo significant-position marking; statistical analysis |
| `MSA_Position_Mapping_ALL_to_Bamboo.csv` | `high_freq_analysis/` | ALL_MSA → Bamboo MSA column mapping |
| `chart_data.json` | generated by `prepare_chart_data.py` | Donut/Venn chart data |

## Technology stack

- **Plotting engine**: Matplotlib (primary) + `matplotlib_venn` (Venn diagrams)
- **Data**: NumPy, Pandas, SciPy (`fisher_exact`)
- **Image formats**: SVG (vector, editable) + PDF (print) + PNG (preview)
- **Style**: Nature journal palette (Navy/Slate/Ochre/Terracotta), Arial font, `svg.fonttype=none` for Illustrator editability

## Regenerate figures

```bash
# Donut/Venn figures (3 images)
python 脚本/generate_final_svgs.py

# All-species sequence logo
python 脚本/generate_sequence_logo_corrected.py

# Bamboo sequence logo (18 sequences)
python 脚本/generate_sequence_logo_bamboo.py

# Run Bamboo Fisher test (generates CSV)
python 脚本/bamboo_18seq_fisher.py
```

---

## Citation

If you use this pipeline or its scripts in your work, please cite the accompanying manuscript.

## License

To be determined by the authors; contact the repository owner for reuse permissions.
