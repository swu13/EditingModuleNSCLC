# NSCLC  Analysis Pipeline

This project integrates RNA editing events, gene expression profiles, and clinical data from TCGA lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) to perform RNA editing network construction, module identification, immune infiltration analysis, and predictive modeling.

---

## Project Structure

```
/data_d/ZJJin/NSCLC/
├── preprocess_data.sh          # Data download, merging, and annotation
├── figure1.R                   # Figure 1: RNA editing events predict EIF2AK2 expression
├── figure2.R                   # Figure 2: Editing network, module identification, survival analysis, lymphatic metastasis
├── figure3.R                   # Figure 3: Gene network and immune correlation for a target module
├── figure4.R                   # Figure 4: Editing correlations, PCA, mediation analysis, immune subtypes
├── figure5_build_features.R            # Build four feature sets for modeling
├── train_model.py              # Module with cross-validation
├── Figure5_plot.R              # Final plots: radar chart, ridge plot, DeLong test
├── Expression/                 # Expression, clinical, and survival data (preprocessed)
├── Editing/                    # RNA editing event data (preprocessed)
├── Results/                    # Intermediate results (annotations, network edges, modules)
├── Figure2_output/
├── Figure3_output/
├── Figure4_output/
└── Figure5/
    ├── features/               # Feature set CSVs
    ├── cv_results/             # Cross-validation results
    └── plots/                  # Final figures
```

---

## Scripts

### 1. `preprocess_data.sh` — Data Preprocessing

Downloads and processes raw data from TCGA, including TPM expression matrices, clinical info, survival data, and Gencode annotation files (probemap and GFF3). Merges LUAD and LUSC expression profiles, converts Ensembl IDs to gene symbols, and retains tumor samples only (`.01A`). Filters RNA editing events with a valid sample ratio > 0.2 (missing values filled with 0) and annotates editing sites to genes using the GFF3 file.

**Output**
- `Expression/Exp_TPM_data.csv`, `combined_clinical_01A.tsv`, `TCGA-lung.survival.tsv`
- `Editing/Filled_lung_0.2.txt`
- `Results/editing_events_unique_gene_annotation.txt`

**Dependencies**: Bash, wget, gunzip; R: `dplyr`, `rtracklayer`, `GenomicRanges`

```bash
bash preprocess_data.sh /data_d/ZJJin/NSCLC
```

---

### 2. `figure1.R` — Editing Events Predict EIF2AK2 Expression

Selects editing events associated with EIF2AK2 and builds univariate and multivariate linear regression models with 10-fold cross-validation repeated 100 times. Outputs scatter plots, coefficient plots, MSE comparison boxplots, and predicted vs. actual scatter plots.

**Input**: `Expression/Exp_tpm_01A_data.csv`, `Editing/Filled_lung_0.2.txt`, `Results/editing_events_unique_gene_annotation.txt`

**Output**: `Figure1_output/`

**Dependencies**: `dplyr`, `tidyr`, `ggplot2`, `caret`, `ggpubr`, `gridExtra`, `doParallel`, `forcats`

```bash
Rscript figure1.R /data_d/ZJJin/NSCLC
```

---

### 3. `figure2.R` — Editing Network, Module Identification, and Survival Analysis

Computes Pearson correlations between editing events to build a co-editing network, applies PageRank for node filtering, and uses Walktrap clustering to identify modules (size 30–100). Performs PCA per module and extracts PC1 for Kaplan-Meier survival analysis (high vs. low PC1 groups) and Cox regression. Also tests the association between module PC1 and lymph node metastasis (N0 vs. N+).

**Input**: `Editing/Filled_lung_0.2.txt`, `Results/editing_events_unique_gene_annotation.txt`, `Expression/Exp_TPM_data.csv`, `Expression/TCGA-lung.survival.tsv`, `Expression/combined_clinical_01A.tsv`

**Output**: `Figure2_output/` (module PC matrix, KM curves, Cox results); `Results/filtered_edges.csv`, `Results/merged_df_filtered.rds`

**Dependencies**: `igraph`, `survival`, `survminer`, `ggplot2`

```bash
Rscript figure2.R /data_d/ZJJin/NSCLC
```

---

### 4. `figure3.R` — Gene Network and Immune Correlation for Target Module

Extracts editing events from a specified module (default: module 14), annotates them to genes, and constructs a gene–gene co-editing network. Analyzes correlations between module PC1 and immune cell infiltration scores (nTreg, iTreg, Tex). Reproduces Cox regression results and KM curves from `figure2.R` for the target module.

**Input**: `Editing/Filled_lung_0.2.txt`, `Results/editing_events_unique_gene_annotation.txt`, `Results/merged_df_filtered.rds`, `Results/filtered_edges.csv`, `Figure2_output/module_vectors.csv`, `Results/ImmuCellAI.csv`

**Output**: `Figure3_output/` (gene network PDF, immune scatter plots, t-test results, Cox and KM outputs)

**Dependencies**: `igraph`, `ggplot2`, `ggpubr`

```bash
Rscript figure3.R /data_d/ZJJin/NSCLC          # default module 14
Rscript figure3.R /data_d/ZJJin/NSCLC 5        # specify module 5
```

---

### 5. `figure4.R` — Editing Correlations, Mediation Analysis, and Immune Subtypes

Computes Pearson correlations between editing events in EIF2AK2 and PRKAR2A. Extracts PC1 from editing event PCA and integrates gene expression (EIF2AK2, PRKAR2A, DDX3X). Performs t-tests stratified by DDX3X expression level and mediation analysis (DDX3X → editing PC1 → target gene expression). Tests module PC1 differences across immune subtypes (C2 vs. others) using Wilcoxon/t-test and violin plots.

**Input**: `Editing/Filled_lung_0.2.txt`, `Expression/Exp_TPM_data.csv`, `Expression/mmc2.xlsx` (optional)

**Output**: `Figure4_output/` (correlation tables, mediation RDS objects, violin plots)

**Dependencies**: `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, `corrplot`, `mediation`, `readxl`, `viridis`

```bash
Rscript figure4.R /data_d/ZJJin/NSCLC
```

---

### 6. `build_features.R` — Feature Set Construction

Selects features significantly associated with lymph node metastasis from module PC1 values, expression PCA components (top 30 PCs), and clinical variables. Builds four feature sets for downstream modeling:

| Group | Features |
|-------|----------|
| Group 1 | Clinical only (gender, age, subtype, T stage) |
| Group 2 | Clinical + RNA editing module PC1 |
| Group 3 | Clinical + expression PCA components |
| Group 4 | Clinical + RNA editing + expression PCA |

**Input**: `Figure2_output/module_vectors.csv`, `Expression/Exp_TPM_data.csv`, `Expression/combined_clinical_01A.tsv`

**Output**: `Figure5/features/Group1_Clinical.csv`, `Group2_Clinical_Edit.csv`, `Group3_Clinical_Exp.csv`, `Group4_All_Features.csv`

**Dependencies**: `dplyr`, `tidyr`

```bash
Rscript build_features.R /data_d/ZJJin/NSCLC
```

---

### 7. `train_model.py` — AutoGluon Modeling

Runs stratified 5-fold cross-validation (3 train / 1 validation / 1 test, 20 combinations) on all four feature sets using AutoGluon. Saves per-fold metrics, test set predicted probabilities, feature importances, and a performance summary per model.

**Input**: `Figure5/features/Group*.csv`

**Output** (under `Figure5/cv_results/<group>/`)
- `all_cv_results.csv` — per-fold metrics for each model
- `model_performance_summary.csv` — mean ± SD across folds
- `all_test_predictions.csv` — test set prediction probabilities
- `feature_importance/` — feature importance CSVs

**Dependencies**: Python 3.8+; `autogluon`, `pandas`, `numpy`, `scikit-learn`

```bash
python train_model.py --base_dir /data_d/ZJJin/NSCLC
```

---

### 8. `plot_results.R` — Summary Figures

Reads performance summaries across all four feature groups and generates a radar chart comparing test AUC across models. Plots a ridge plot showing AUC gain attributable to RNA editing features (NeuralNetFastAI model). Performs DeLong tests comparing clinical vs. clinical+RNA editing, and clinical+expression vs. all features.

**Input**: `Figure5/cv_results/*/model_performance_summary.csv`, `all_cv_results.csv`, `all_test_predictions.csv`

**Output**: `Figure5/plots/Figure5_radar.png`, `Figure5_ridgeplot.pdf`, `DeLong_test_results.csv`

**Dependencies**: `tidyverse`, `ggridges`, `ggsci`, `pROC`, `scales`

```bash
Rscript plot_results.R /data_d/ZJJin/NSCLC
```

---

## Recommended Execution Order

```
1. preprocess_data.sh    →  Generate all base data
2. figure2.R             →  Network construction and module identification  
3. figure1.R             →  (Optional) EIF2AK2 prediction analysis           
4. figure4.R             →  Mediation and immune subtype analysis            
5. figure3.R             →  Target module analysis        (requires figure2.R output)
6. build_features.R      →  Prepare modeling feature sets (requires figure2.R output)
7. train_model.py        →  Module (Python environment required)
8. plot_results.R        →  Generate final Figure 5 plots
```

