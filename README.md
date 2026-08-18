# Topology-aware thermodynamics for DNA probe design under fixed stringency

**Manuscript:** Van der Snickt et al., *Topology-aware thermodynamics for DNA probe design under fixed stringency.

**Repository:** (https://github.com/OncoRNALab/Topology-aware-thermodynamics-for-DNA-probe-design-)

---

## Overview

This repository contains the code, pre-computed data, and supplementary tables needed to reproduce the figures and analyses reported in the manuscript. The central contribution is a **retained-box framework** that bridges topology-only scoring (S_ECI) and scalar thermodynamics (nearest-neighbor ΔG) into a single partition function, **Z_box,NN**, which includes S_ECI while capturing spatial information not included in scalar ΔG.

Three scoring layers are implemented in the standalone `zbox_nn_engine` package:

| Score | What it captures |
|-------|-----------------|
| **S_ECI** | Topology only — sum of squared retained island lengths |
| **Scalar NN ΔG** | Thermodynamics only — full-duplex nearest-neighbor free energy |
| **Z_box,NN** | Both — Boltzmann-weighted NN free energies over all retained boxes ≥ k |

---

## Repository structure

```
ECI_02/
├── README.md                              ← This file
├── supplementary_tables/                  ← All 16 CSV files (S1–S16)
├── R/                                     ← Figure generation scripts
│   ├── Figure 2 and 3.R                   ← Figure 2 (6-panel empirical proof) and Figure 3 (4-panel portrait)
├── python/                                ← Pipeline & baseline models
│   ├── generate_supplementary_tables.py   ← Regenerates S12–S16 from raw data
│   └── dg_model.py                        ← Scalar NN baseline (reference)
├── zbox_nn_engine/                        ← Z_box/N_box engine (pip-installable)
│   ├── README.md                          ← Engine documentation
│   ├── __init__.py
│   ├── nn_parameters.py                   ← SantaLucia 2004 unified NN parameters
│   ├── island_detection.py                ← Retained island identification
│   ├── scores.py                          ← S_ECI, N_box, Z_box,NN computation
│   └── setup.py
└── data/                                  ← Pre-computed data for Figure 2
    └── fig_data/                          ← 11 CSV files for Figure 2 panels
```

---

## Quick start

### Generate figures (R)

```r
# Install required packages (one-time)
install.packages(c("ggplot2", "ggprism", "patchwork", "svglite", "dplyr", "readr"))

# Figure 2 and 3 — reads from supplementary_tables/
setwd("/path/to/Supplementary tables/")
source("R/Figure2.R")
fig2 <- make_figure2()   # produces Fig. 2.pdf, Fig. 2.png, Fig. 2.svg
source("R/Figure3.R")
fig3 <- make_figure3()   # produces Fig. 3.pdf, Fig. 3.png, Fig. 3.svg


### Regenerate computed tables from raw data (Python)

```bash
# Install the engine
pip install -e ./zbox_nn_engine

# Run the pipeline (requires raw data files — see Raw Data below)
python python/generate_supplementary_tables.py \
    --raw-data-dir /path/to/raw_data/ \
    --output-dir supplementary_tables/
```

This regenerates S12–S16. Tables S1–S11 are literature-derived and ship as static CSVs.

---

## Supplementary table index

| Table | Description |
|-------|-------------|
| S1 | Scoring definitions (formal/retained islands, S_ECI, N_box, Z_box,NN, ΔG_NN) |
| S2 | Evidence-layer map (6 external datasets and their role in the manuscript) |
| S3 | Seringhaus et al. [7] reconstruction (ACT1/HBG2, centered vs staggered) |
| S4 | Deng et al. [8] even-vs-random mismatch aggregate (0–7 mismatches) |
| S5 | Deng et al. [8] MPDNN model improvement (R² 0.637 → 0.880) |
| S6 | Rennie et al. [9] summary (Lmax explains 43% of variance) |
| S7 | Mechanistic and thermodynamic bridge (Naiser, Hadiwikarta) |
| S8 | Stringency calibration (Deng 50/45/42 °C formamide conditions) |
| S9 | Affymetrix fixed-mismatch rows (PM, clustered, distributed probes) |
| S10 | Affymetrix correlation summary (2–4 MM and strict 3-MM subsets) |
| S11 | HPV edge-case audit (7 CP cases with discrepancy flags) |
| S12 | Single-model R² and partial correlations at primary conditions (KRAS 68 °C/1000 mM; BRAF 60 °C/1000 mM) |
| S13 | Per-condition R² and partial correlations across the stringency sweep (16 temperature/salt combinations) |
| S14 | Top-quantile probe identification rates with 95% bootstrap CIs (top 5%, 10%, 25%) |
| S15 | Per-target R² for the 10-target generalization panel |
| S16 | Per-probe scatter data for Figure 3 (1211 probes, 27 columns) |

---

## Z_box,NN calculation engine

The `zbox_nn_engine/` directory is a standalone, pip-installable Python package implementing the retained-box scoring framework. It can be used independently of this repository for future probe-design applications.

```bash
pip install -e ./zbox_nn_engine
```

```python
from zbox_nn_engine import compute_all_scores

scores = compute_all_scores(
    probe_seq="ctcttgcctacgccaca",
    mt_target_seq="ctcttgcctacgccaca",   # mutant target (perfect match)
    wt_target_seq="ctcttgcctacgccacc",   # wild-type target (1 mismatch)
    temperature_C=68,
    salt_mM=1000,
    threshold=2,
)
print(scores["Zbox_ratio"])
```

See [`zbox_nn_engine/README.md`](zbox_nn_engine/README.md) for full API documentation, scoring formulas, and usage examples.

---

## Raw data

Raw hybridization data for the 17-mer KRAS and BRAF probe panels originate from Van der Snickt et al. [11] and are available from that source (https://github.com/OncoRNALab/Probe-enrichment-platform). The following files are required to regenerate S12–S16 via the Python pipeline:

| File | Description |
|------|-------------|
| `Data KRAS c.34 G>T 17mer.txt` | KRAS c.34 G>T hybridization data  |
| `Data BRAF c.1799 T>A 17mer.txt` | BRAF c.1799 T>A hybridization data |
| `Simulated_KRAS_17mer.txt` | Scalar NN binding-free-energy predictions (KRAS) |
| `Simulated_BRAF_17mer.txt` | Scalar NN binding-free-energy predictions (BRAF) |
| `Data 10 prospective targets.txt` | 10-target generalization panel |

External datasets reanalyzed in S1–S11 (Seringhaus, Deng, Rennie, Naiser, Affymetrix, HPV typing) were obtained from the published sources cited in the manuscript reference list.

---
