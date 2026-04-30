# Single-Cell RNA-seq Analysis Pipeline

A Jupyter Notebook implementing a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline using [Scanpy](https://scanpy.readthedocs.io/), from raw 10x Genomics data to cell clustering and marker gene identification.

## Dataset

Mouse brain nuclei isolated and sequenced using the 10x Genomics Chromium platform, comparing three protocols:
- Chromium Nuclei Isolation
- SaltyEZ
- Complex Tissue DP

Input file: `M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5`

## Pipeline Overview

### Module 0 — Configuration & Data Extraction
- Sets up Scanpy global parameters (verbosity, figure DPI, color map)
- Defines working directories for data, figures, and cache
- Loads the filtered feature-barcode matrix from the 10x `.h5` file

### Module 1 — Quality Control & Filtering
- Identifies mitochondrial genes (`mt-` prefix)
- Computes QC metrics: genes per cell, total counts, and % mitochondrial reads
- Filters out low-quality cells and likely doublets:
  - Minimum 500 genes / cell
  - Maximum 6,500 genes / cell
  - Minimum 500 total counts / cell
  - Less than 5% mitochondrial content
- Filters genes expressed in fewer than 3 cells
- Normalizes counts to 10,000 reads per cell and applies log1p transformation
- Selects highly variable genes to reduce noise

### Module 2 — Dimensionality Reduction & Clustering
- Regresses out confounding factors (`total_counts`, `pct_counts_mt`)
- Scales the data (max value = 10)
- Runs PCA (40 components) and plots elbow curves to determine optimal PCs
- Builds a k-nearest neighbor graph (k=15, top 9 PCs)
- Computes UMAP embedding
- Performs Leiden clustering (resolution = 0.3)

### Module 4 — Marker Gene Identification
- Uses Wilcoxon rank-sum test to find differentially expressed genes across Leiden clusters
- Extracts and prints the top 10 marker genes per cluster
- Visualizes the expression of a gene of interest (`Egfr`) on the UMAP and as a dot plot

## Output Files

| File | Description |
|---|---|
| `figures/_qc_metrics.png` | Violin plots of QC metrics before filtering |
| `figures/_filtered.png` | Violin plots of QC metrics after filtering |
| `figures/_elbow_plot.png` | PCA variance ratio (log scale) |
| `figures/_elbow_plot_linear.png` | PCA variance ratio (linear scale) |
| `figures/_clusters.png` | UMAP colored by Leiden cluster |
| `figures/_rank_genes.png` | Top marker genes per cluster |
| `figures/_Egfr_expression.png` | UMAP colored by Egfr expression |
| `figures/_Egfr_preference.png` | Dot plot of Egfr expression per cluster |

## Requirements

```
scanpy
pandas
numpy
matplotlib
anndata
```

Install dependencies:

```bash
pip install scanpy pandas numpy matplotlib
```

## Usage

1. Place the `.h5` input file inside the `data/` directory.
2. Open the notebook and run all cells sequentially.
3. Figures will be saved automatically to the `figures/` directory.

## Project Structure

```
.
├── data/
│   └── M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5
├── figures/
├── notebook.ipynb
└── README.md
```
