# Selective Sweep Simulation

This repository contains all scripts and configurations used to simulate genealogies under selection, infer tree-based statistics, and analyze signals of natural selection using tree imbalance and pairwise distance metrics.  
It includes both **Slim-based forward simulations** and **msprime-based coalescent simulations**, along with R-based post-analysis pipelines.

---

## 📁 Repository Structure

### `SliM/`
Contains the **forward-time simulation pipeline** and replication of the paper  
[**“Robust detection of natural selection using a probabilistic model of tree imbalance.”**](https://academic.oup.com/genetics/article/220/3/iyac009/6511494?login=false)

- **`PSlim.py`** — main Python driver for running SLiM simulations and recapitating trees with msprime.  
- **`recap.py`** — post-simulation recapitation script (extends SLiM output to coalescent-compatible `.trees` files).  
- **`Slim.txt`** — core SLiM script defining population demography, mutation, and selection parameters.  
- **`Simulation---Constant_Directional_Selection.ipynb`** — Jupyter workflow used to control SLiM simulations and trigger recapitation.  
- **`Calculation---Constant_Directional_Selection.ipynb`** — Jupyter workflow used to calculate summary statsitics and produce visuals.
- **`binary_to_newick.ipynb`** — Jupyter workflow used to convert `.trees` files to newick format for analysis in R. 
- **`env_sim312.yml` / `env_sweep312.yml`** — reproducible environment specifications (see below).

This module reproduces [**Dilber & Terhorst’s** genome-scan framework](https://github.com/jthlab/bim-paper/tree/main) and computes tree-based statistics (e.g., `bsfs`, `TajD`, `btree`, `iColless`) to visualize local deviations under selection.

---

### `Simulation/`
Contains **msprime-based backwards-in-time coalescent simulations** for controlled experiments on tree-based distance metrics.

- **`run.py`** — simulates multiple evolutionary scenarios (varying selection strength and allele frequency) and extracts genealogies at multiple genomic positions.  
- **`genomic_position.py`** — performs fine-scale chromosome simulations and exports genealogies with explicit positional metadata.

These simulations were used for the initial analysis of the F-matrix and distance metric behaviour under selection.

---

### `Analysis/`
Contains **R-based analysis and visualization scripts** for post-simulation inference, using the [`phylodyn`](https://cran.r-project.org/package=phylodyn) framework and custom tools. Most files can be found in `phylodyn` folder.

#### Part 1. Analysis of Pairwise Distance Distribution Under Different Scenarios
- **`phylodyn_rework.Rmd`** — main analysis pipeline that computes pairwise tree distances, performs MDS embeddings, and visualizes distributions under selection vs. neutrality.  
- **`chromosome_analysis.Rmd`** — examines individual chromosomes, plotting total distance per locus and detecting local sweep signatures.  
- **`giant_MDS.Rmd`** — builds large-scale MDS plots using the full all-vs-all distance matrix.

Outputs include `.json` distance data, MDS coordinates, and summary plots used for main figures.

#### Part 2. Leverage Pairwise Distances to Detect Natural Selection Along Chromosomes
- **`load_and_store_data.Rmd`** — Load .tree files simulated using Slim, calculate the F/FW-matrices for each genealogy, the weighted average F/FW-matrices for each genome, and estimated empirical mean F/FW-matrix under neutrality. Saves the results to `data` folder (not included in this repo).
- **`chromosome_wide_Fmat_analysis`** - Compute the distance of averaged F-matrices to the empirical mean F-matrix, compares the distribution, calculates ROC_AUC scores, etc. The results are stored in `roc` folder (not included in this repo) to be loaded in `Calculation---Constant_Directional_Selection.ipynb`.
- **`chromosome_wide_FWmat_analysis`** - Compute the distance of averaged FW-matrices to the empirical mean FW-matrix, compares the distribution, calculates ROC_AUC scores, etc. The results are stored in `roc` folder (not included in this repo) to be loaded in `Calculation---Constant_Directional_Selection.ipynb`.

---

### `Results/`
Contains serialized simulation outputs and analysis-ready JSON files.  
*(Note: these large files are excluded from GitHub but can be regenerated from the above scripts.)*

---

## ⚙️ Environment Setup

The SliM project (simulation using Slim) uses **two reproducible conda environments** to maintain compatibility between modern SLiM/tskit workflows and legacy BIM/statistics pipelines.

| Environment | Purpose | Key Packages |
|--------------|----------|---------------|
| `sim312` | SLiM → recap → msprime (modern tree sequences) | `tskit >= 0.6`, `pyslim >= 1.0` |
| `sweep312` | BIM estimation (legacy tskit ABI) | `tskit == 0.5.8`, `BIM` |

To create them from scratch:
```bash
conda env create -f SliM/env_sim312.yml
conda env create -f SliM/env_sweep312.yml