# Selective Sweep Simulation

This repository contains all scripts and configurations used to simulate genealogies under selection, infer tree-based statistics, and analyze signals of natural selection using tree imbalance and pairwise distance metrics.  
It integrates both **SLiM-based forward simulations** and **msprime-based coalescent simulations**, along with R-based post-analysis pipelines.

---

## 📁 Repository Structure

### `SliM/`
Contains the **forward-time simulation pipeline** and replication of the paper  
**“Robust detection of natural selection using a probabilistic model of tree imbalance.”**

- **`PSlim.py`** — main Python driver for running SLiM simulations and recapitating trees with msprime.  
- **`recap.py`** — post-simulation recapitation script (extends SLiM output to coalescent-compatible `.trees` files).  
- **`Slim.txt`** — core SLiM script defining population demography, mutation, and selection parameters.  
- **`Simulation---Constant_Directional_Selection.ipynb`** *(optional)* — Jupyter workflow used to control SLiM simulations, trigger recapitation, and compute genome-scan statistics.  
- **`env_sim312.yml` / `env_sweep312.yml`** — reproducible environment specifications (see below).

This module reproduces **Dilber & Terhorst’s** genome-scan framework and computes tree-based statistics (e.g., `bsfs`, `TajD`, `btree`, `iColless`) to visualize local deviations under selection.

---

### `Simulation/`
Contains **msprime-based coalescent simulations** for controlled experiments on tree-based distance metrics.

- **`run.py`** — simulates multiple evolutionary scenarios (varying selection strength and allele frequency) and extracts genealogies at multiple genomic positions.  
- **`genomic_position.py`** — performs fine-scale chromosome simulations and exports genealogies with explicit positional metadata.

These simulations form the basis for downstream **phylogenetic distance analyses** and the development of the **F-matrix–based distance metric**.

---

### `Analysis/`
Contains **R-based analysis and visualization scripts** for post-simulation inference, using the [`phylodyn`](https://cran.r-project.org/package=phylodyn) framework and custom tools.

- **`phylodyn_rework.Rmd`** — main analysis pipeline that computes pairwise tree distances, performs MDS embeddings, and visualizes distributions under selection vs. neutrality.  
- **`chromosome_analysis.Rmd`** — examines individual chromosomes, plotting total distance per locus and detecting local sweep signatures.  
- **`giant_MDS.Rmd`** — builds large-scale MDS plots using the full all-vs-all distance matrix.

Outputs include `.json` distance data, MDS coordinates, and summary plots used for main figures.

---

### `Results/`
Contains serialized simulation outputs and analysis-ready JSON files.  
*(Note: these large files are excluded from GitHub but can be regenerated from the above scripts.)*

---

## ⚙️ Environment Setup

The SliM project uses **two reproducible conda environments** to maintain compatibility between modern SLiM/tskit workflows and legacy BIM/statistics pipelines.

| Environment | Purpose | Key Packages |
|--------------|----------|---------------|
| `sim312` | SLiM → recap → msprime (modern tree sequences) | `tskit >= 0.6`, `pyslim >= 1.0` |
| `sweep312` | BIM estimation (legacy tskit ABI) | `tskit == 0.5.8`, `BIM` |

To create them from scratch:
```bash
conda env create -f SliM/env_sim312.yml
conda env create -f SliM/env_sweep312.yml