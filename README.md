# Selective Sweep Simulation

This repository contains all scripts and configurations used to simulate genealogies under selection, infer tree-based statistics, and analyze signals of natural selection using tree imbalance and pairwise distance metrics.  
It includes both **SLiM-based forward simulations** and **msprime-based coalescent simulations**, along with R-based post-analysis pipelines.

---

## 🔄 Workflow

The analysis pipeline follows this sequential workflow:

1. **Simulate Trees** (`Simulation/Slim(primary)/Simulation.ipynb`)
   - Run SLiM forward-time simulations to generate tree sequences
   - Uses the `sim312` conda environment
   - Outputs binary `.trees` files stored in `trees/` directory

2. **Convert to Newick Format** (`Simulation/Slim(primary)/binary_to_newick.ipynb`)
   - Convert binary tree sequence files to Newick format (`.tree` files)
   - Generate CSV files with genomic interval breakpoints
   - Should be run **once** after all simulations are complete
   - Outputs to `trees_newick/` directory

3. **Load Trees and Calculate F-Matrices** (`Analysis/phylodyn/load_and_store_data.Rmd`)
   - Load Newick tree files in R
   - Convert genealogies into F-matrices and weighted F-matrices (FW-matrices)
   - Calculate weighted average F/FW-matrices for each genome
   - Estimate empirical mean F/FW-matrix under neutrality
   - Store results in `data/` folder

4. **Calculate Distance Metrics** (`Analysis/phylodyn/chromosome_wide_Fmat_analysis.Rmd` and `Analysis/phylodyn/chromosome_wide_FWmat_analysis.Rmd`)
   - Compute distance of averaged F-matrices to empirical mean F-matrix
   - Compute distance of averaged FW-matrices to empirical mean FW-matrix
   - Compare distributions and calculate ROC_AUC scores
   - Results stored in `roc/` folder

5. **Calculate Summary Statistics and Compare Results** (`Analysis/Calculation.ipynb`)
   - Use BIM (Software for β-Imbalance) package to calculate additional summary statistics
   - Compare results from F-matrix analysis with BIM-based statistics
   - Uses the `sweep312` conda environment
   - Must be run after completing simulations and R scripts for F-matrix calculations

---

## 📁 Repository Structure

### `Simulation/Slim(primary)/`
Contains the **forward-time simulation pipeline** and replication of the paper  
[**"Robust detection of natural selection using a probabilistic model of tree imbalance."**](https://academic.oup.com/genetics/article/220/3/iyac009/6511494?login=false)

- **`Simulation.ipynb`** — Main Jupyter notebook to run SLiM simulations (use `sim312` conda environment)
- **`binary_to_newick.ipynb`** — Converts binary `.trees` files to Newick format for R analysis
- **`PSlim.py`** — Main Python driver for running SLiM simulations and recapitating trees with msprime
- **`recap.py`** — Post-simulation recapitation script (extends SLiM output to coalescent-compatible `.trees` files)
- **`Slim.txt`** — Core SLiM script defining population demography, mutation, and selection parameters
- **`env_sim312.yml`** — Conda environment specification for simulations

This module reproduces [**Dilber & Terhorst's** genome-scan framework](https://github.com/jthlab/bim-paper/tree/main) and computes tree-based statistics (e.g., `bsfs`, `TajD`, `btree`, `iColless`) to visualize local deviations under selection.

---

### `Simulation/msprime(secondary)/`
Contains **msprime-based backwards-in-time coalescent simulations** for controlled experiments on tree-based distance metrics.

- **`run.py`** — Simulates multiple evolutionary scenarios (varying selection strength and allele frequency) and extracts genealogies at multiple genomic positions

These simulations were used for the initial analysis of the F-matrix and distance metric behaviour under selection.

---

### `Analysis/`
Contains **R-based analysis and visualization scripts** for post-simulation inference, using the [`phylodyn`](https://cran.r-project.org/package=phylodyn) framework and custom tools.

#### Key Analysis Scripts
- **`Calculation.ipynb`** — Main analysis notebook using BIM package to calculate summary statistics and compare results (use `sweep312` conda environment)
- **`PSlim2.py`** — Python utilities for analysis workflow

#### `Analysis/phylodyn/`
Contains R Markdown scripts for F-matrix analysis:

- **`load_and_store_data.Rmd`** — Load Newick tree files, calculate F/FW-matrices for each genealogy, compute weighted averages, and estimate empirical mean matrices. Saves results to `data/` folder
- **`load_and_store_data_utils.R`** — Utility functions for data loading and processing
- **`chromosome_wide_Fmat_analysis.Rmd`** — Compute distance of averaged F-matrices to empirical mean, compare distributions, calculate ROC_AUC scores. Results stored in `roc/` folder
- **`chromosome_wide_FWmat_analysis.Rmd`** — Compute distance of averaged FW-matrices to empirical mean, compare distributions, calculate ROC_AUC scores. Results stored in `roc/` folder

---

### `Results/`
Contains serialized simulation outputs and analysis-ready JSON files.  
*(Note: these large files are excluded from GitHub but can be regenerated from the above scripts.)*

---

## ⚙️ Environment Setup

The project uses **two reproducible conda environments** to maintain compatibility between modern SLiM/tskit workflows and legacy BIM/statistics pipelines.

| Environment | Purpose | Key Packages |
|--------------|----------|---------------|
| `sim312` | SLiM → recap → msprime (modern tree sequences) | `tskit >= 0.6.0`, `pyslim >= 1.0.6` |
| `sweep312` | BIM estimation and analysis (legacy tskit ABI) | `tskit == 0.5.8`, `BIM` |

To create them from scratch:
```bash
conda env create -f Simulation/Slim(primary)/env_sim312.yml
conda env create -f Analysis/env_sweep312.yml