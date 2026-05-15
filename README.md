# Selective Sweep Simulation

This repository contains all scripts and configurations used to simulate genealogies under selection, infer tree-based statistics, and analyze signals of natural selection using tree imbalance and pairwise distance metrics.  
It includes both **SLiM-based forward simulations** and **msprime-based coalescent simulations**, along with R-based post-analysis pipelines.

---

## ⚙️ Environment Setup

The project uses **three reproducible conda environments** to maintain compatibility between modern SLiM/tskit workflows, legacy BIM/statistics pipelines, and the 1000 Genomes real-data analysis.


| Environment   | Purpose                                          | Key Packages                                                  |
| ------------- | ------------------------------------------------ | ------------------------------------------------------------- |
| `sim312`      | SLiM → recap → msprime (modern tree sequences)   | `tskit >= 0.6.0`, `pyslim >= 1.0.6`                           |
| `sweep312`    | BIM estimation and analysis (legacy tskit ABI)   | `tskit == 0.5.8`, `BIM`                                       |
| `genomes1000` | 1000 Genomes tree sequence analysis and plotting | `tskit`, `tszip`, `msprime`, `pyslim`, `pandas`, `matplotlib` |


To create them from scratch:

```bash
conda env create -f "Simulation/Slim(primary)/env_sim312.yml"
conda activate sim312
conda install ipykernel
python -m ipykernel install --user --name sim312 --display-name "sim312"
```

```bash
conda env create -f Analysis/env_sweep312.yml
conda activate sweep312
conda install ipykernel
python -m ipykernel install --user --name sweep312 --display-name "sweep312"
```

```bash
conda env create -f 1000Genomes/genomes1000.yml
conda activate genomes1000
conda install ipykernel
python -m ipykernel install --user --name genomes1000 --display-name "genomes1000"
```

---

## 🔄 Workflow (Simulation based analysis)

The analysis pipeline for the *"Detecting selection at genome-wide level"* and *"Localizing selection"* sections in the paper follows this sequential workflow:

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

## 📁 Repository Structure (Simulation based analysis)

### `Simulation/Slim(primary)/`

Contains the **forward-time simulation pipeline** and replication of the paper  
["Robust detection of natural selection using a probabilistic model of tree imbalance."](https://academic.oup.com/genetics/article/220/3/iyac009/6511494?login=false)

- `**Simulation.ipynb`** — Main Jupyter notebook to run SLiM simulations (use `sim312` conda environment)
- `**binary_to_newick.ipynb`** — Converts binary `.trees` files to Newick format for R analysis
- `**PSlim.py`** — Main Python driver for running SLiM simulations and recapitating trees with msprime
- `**recap.py`** — Post-simulation recapitation script (extends SLiM output to coalescent-compatible `.trees` files)
- `**Slim.txt`** — Core SLiM script defining population demography, mutation, and selection parameters
- `**env_sim312.yml*`* — Conda environment specification for simulations

This module reproduces **[Dilber & Terhorst's** genome-scan framework](https://github.com/jthlab/bim-paper/tree/main) and computes tree-based statistics (e.g., `bsfs`, `TajD`, `btree`, `iColless`) to visualize local deviations under selection.

---

### `Simulation/msprime(secondary)/`

Contains **msprime-based backwards-in-time coalescent simulations** for controlled experiments on tree-based distance metrics.

- `**run.py`** — Simulates multiple evolutionary scenarios (varying selection strength and allele frequency) and extracts genealogies at multiple genomic positions

These simulations were used for the initial analysis of the F-matrix and distance metric behaviour under selection.

---

### `Analysis/`

Contains **R-based analysis and visualization scripts** for post-simulation inference, using the `[phylodyn](https://cran.r-project.org/package=phylodyn)` framework and custom tools.

#### Key Analysis Scripts

- `**Calculation.ipynb`** — Main analysis notebook using BIM package to calculate summary statistics and compare results (use `sweep312` conda environment)
- `**PSlim2.py`** — Python utilities for analysis workflow

#### `Analysis/phylodyn/`

Contains R Markdown scripts for F-matrix analysis:

- `**load_and_store_data.Rmd`** — Load Newick tree files, calculate F/FW-matrices for each genealogy, compute weighted averages, and estimate empirical mean matrices. Saves results to `data/` folder
- `**load_and_store_data_utils.R`** — Utility functions for data loading and processing
- `**chromosome_wide_Fmat_analysis.Rmd`** — Compute distance of averaged F-matrices to empirical mean, compare distributions, calculate ROC_AUC scores. Results stored in `roc/` folder
- `**chromosome_wide_FWmat_analysis.Rmd`** — Compute distance of averaged FW-matrices to empirical mean, compare distributions, calculate ROC_AUC scores. Results stored in `roc/` folder

---

## 🔄 Workflow (Real data analysis)

The analysis pipeline for the *"Real Data Analysis"* section in the paper follows this sequential workflow:

1. Create and activate the `genomes1000` environment:
2. Download the inferred tree sequence data from [Kelleher et al.](https://zenodo.org/records/3051855) and place the chr2 file at `1000Genomes/1kg_chr2.trees.tsz`.
3. Run `1000Genomes/load_1kg_tsz.ipynb` from inside the `1000Genomes/` directory. This generates the `.trees`, `.tree`, and `_breaks.csv` inputs used by the R analyses.
4. Run the R Markdown notebooks that create the result CSVs consumed by `1000GenomesAnalysis.ipynb`. For the main chr2 figures, run:
  - `chr2_50samples_d1f.Rmd`
  - `chr2_50samples_d1fw.Rmd`
  - `chr2_2inds_allpops_d1f.Rmd`
  - `chr2_52samples_vs_allpops2inds_d1f.Rmd`
  - `chr2_52samples_signed_d1f.Rmd`
5. Open `1000Genomes/1000GenomesAnalysis.ipynb` with the `genomes1000` kernel and run the notebook. The notebook expects the generated CSVs in `1000Genomes/results/`.
6. For the ARGweaver sections at the end of the notebook, first run `raw_data/argweaver_pipeline.ipynb` or provide equivalent ARGweaver outputs under `1000Genomes/raw_data/data/runs/`, then run `chr2_20samples_d1f_argweaver.Rmd`. The empirical-neutral ARGweaver plots also require the neutral SLiM conversion outputs described in the later cells of `1000GenomesAnalysis.ipynb`.

---

## 📁 Repository Structure (Real data analysis)

### `1000Genomes/`

Contains the empirical 1000 Genomes workflow used to create chromosome 2 tree-based statistics and plots. Before running this folder, download the inferred 1000 Genomes tree sequence data from Zenodo:

[https://zenodo.org/records/3051855](https://zenodo.org/records/3051855)

Place the chromosome tree sequence archive where `load_1kg_tsz.ipynb` expects it, for example:

```text
1000Genomes/1kg_chr2.trees.tsz
```

The notebook can load the file directly with `tskit.load(...)`, or fall back to explicit `tszip.decompress(...)` when needed.

#### Main 1000 Genomes Scripts

- `**genomes1000.yml**` - Conda environment for the 1000 Genomes workflow, including `tskit`, `tszip`, `msprime`, `pyslim`, `pandas`, `matplotlib`, and Jupyter.
- `**load_1kg_tsz.ipynb**` - Loads the Zenodo-inferred `.trees.tsz` chromosome tree sequence, samples population-specific haplotypes, restricts to the chr2 LCT analysis interval, simplifies each subset, writes `.trees` files, converts local trees to one-Newick-per-line `.tree` files, and writes breakpoint CSVs under `1000Genomes/results/`. It also creates the `chr2_2inds_allpops` combined all-population subset.
- `**chr2_50samples_d1f.Rmd**` - Computes per-tree D1F distances for chr2 population tree files against a Kingman reference. It also builds population-specific weighted F references, grouped CEU/JPT references, CEU/JPT `F_ST^Tree` summaries, and binned between-population D1F outputs.
- `**chr2_50samples_d1fw.Rmd**` - Computes D1FW distances using FW matrices and a shared chromosome-wide weighted FW reference across the selected populations.
- `**chr2_2inds_allpops_d1f.Rmd**` - Computes per-tree D1F for the `chr2_2inds_allpops` subset created by `load_1kg_tsz.ipynb`.
- `**chr2_52samples_vs_allpops2inds_d1f.Rmd**` - Builds/caches F matrices for population-specific trees and the all-population two-individual subset, computes binned D1F(population vs allpops), and produces local `F_ST^Tree` CSVs for CEU comparisons.
- `**chr2_52samples_signed_d1f.Rmd**` - Computes signed D1F for selected 52-tip population tree sets and the all-population subset using Kingman, balanced, and unbalanced reference matrices.
- `**1000GenomesAnalysis.ipynb**` - Main plotting and analysis notebook. It reads the CSV outputs from the R Markdown scripts, creates D1F and D1FW location plots, computes changepoint segment p-values, summarizes `F_ST^Tree`, and includes the ARGweaver CEU/JPT comparison sections.

#### ARGweaver / Raw-Data Helpers

- `**raw_data/argweaver_pipeline.ipynb**` - End-to-end ARGweaver workflow from raw 1000 Genomes Phase 3 VCFs to posterior Newick tree files. It selects samples, subsets a VCF with `bcftools`, converts phased genotypes to `.sites`, runs `arg-sample`, and writes `.tree` plus breakpoint outputs.
- `**raw_data/vcf2sites.py**` - Command-line helper that converts a phased, population-subsetted VCF into ARGweaver `.sites` format for a requested chromosome interval.
- `**raw_data/gz2trees.ipynb**` - Notebook helper that converts ARGweaver `.smc.gz` posterior samples into cleaned Newick `.tree` files and breakpoint tables.
- `**raw_data/gz2trees.py**` - Minimal script for converting ARGweaver `.smc.gz` files to `.trees` files through `argweaver.smc.smc2tskit`; useful as a conversion experiment or lightweight helper.
- `**chr2_20samples_d1f_argweaver.Rmd**` - Computes post-burn-in per-tree D1F for ARGweaver CEU/JPT outputs, first against a theoretical Kingman reference and then against empirical neutral references from matched neutral SLiM/ARGweaver-style runs.

