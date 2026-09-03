# AF3 Inputs — AlphaFold3 Prediction (Representative Examples)

This directory holds **representative AlphaFold3 (AF3, `alphafold3` v3.0.1) input JSON
files** used for structure prediction in the manuscript. The full prediction set (HCG
Round-2 trimer complex, 170 variants; SH3 combinatorial library, 313 variants) is large
because every input embeds per-chain MSAs and templates. To keep this repository lean, we
deposit one representative input per case study (MSAs/templates stripped; sequences +
`modelSeeds` kept) plus the shared pipeline configuration, and document the fixed seed
scheme below.

## Files

| File | Represents |
|---|---|
| `HCG_WT_original.example.json` | HCG VNAR (Case Study 1) trimer complex, wild-type parent — 5 seeds, 3 chains (α/β hCG antigen + VNAR with C-terminal His/FLAG tags). |
| `SH3_C87S_C93Y.example.json` | SH3 domain (Case Study 2) library variant C87S_C93Y — 5 seeds, 2 chains. |
| `full_pipeline_example.yaml` | Pipeline configuration (AF3 filtering → PyRosetta → AMBER thresholds). |

## AF3 seed scheme

Focused Round-2 mutation panels (HCG 170 / SH3 313 variants) were each run with **five
independent random seeds** per variant:

- HCG trimer inputs: `modelSeeds = [26226, 116740, 288390, 670488, 777573]`
- SH3 library inputs: `modelSeeds = [87231, 49455, 37084, 89841, 63891]`

Seed identity is held constant per case study; only the variant CDR3/domain sequence
differs across inputs. The initial complex (docking-landscape) modelling used an expanded
random-seed scan (input specification available from the authors on request).

All inputs use `"dialect": "alphafold3"`, standard AF3 `sequences` array (`protein`
chains with `id`, `sequence`, `modifications`), and each `modelSeeds` entry yields one
sampled model. The `full_pipeline_example.yaml` thresholds are:

| Parameter | Threshold |
|---|---|
| `plddt_threshold` | 0.7 |
| `ipsae_threshold` | 0.6 |
| `clashes_threshold` | 5 |
| `pdockq_threshold` | 0.2 |
| `iptm_threshold` | 0.6 |

PyRosetta (Part 2): `relax = true`, `dump_top_n = 30`.
AMBER MD (Part 3): `production_ns = 100`, `npt_ns = 1.0`, `tmp = 310.0`,
`forcefield = amber14sb_parmbsc1`.