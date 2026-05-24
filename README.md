# AstroSynapse
### Mapping Neuroligin Isoform Expression in Astrocytic Synaptic Modulation

A computational single-cell RNA-seq pipeline for profiling Neuroligin (Nlgn1/2/3) isoform expression across astrocyte subtypes and building co-expression networks to study astrocytic contributions to tripartite synapse signaling.

---

## Overview

Astrocytes are increasingly recognized as active participants in synaptic transmission through the tripartite synapse. This pipeline characterizes how Neuroligin isoforms — synaptic adhesion molecules typically studied in neurons — are differentially expressed across astrocyte subtypes and co-expressed with other synaptic genes.

**Pipeline steps:**

| Script | Description | Output |
|--------|-------------|--------|
| `01_preprocess.py` | QC filtering, normalization, HVG selection, PCA | `data/01_preprocessed.h5ad` |
| `02_cluster_annotate.py` | Clustering, cell-type annotation, astrocyte extraction | `data/02_astrocytes.h5ad` |
| `03_astrocyte_subcluster.py` | Astrocyte sub-clustering and subtype annotation | `data/03_astrocytes.h5ad` |
| `04_neuroligin_expression.py` | Nlgn1/2/3 expression analysis, DE testing | `figures/04_*.png`, `results/04_*.csv` |
| `05_coexpression.py` | Pearson correlation network, GO enrichment | `figures/05_*.png`, `results/05_*.csv` |

---

## Requirements

```bash
pip install -r requirements.txt
```

**Key dependencies:** `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `networkx`, `statsmodels`, `gprofiler-official`, `leidenalg`

Python 3.8+ recommended.

---

## Setup

**1. Download a dataset**

See `00_download_instructions.txt` for recommended datasets. The primary recommended dataset is:

> Bayraktar et al. (2020) — *Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ transcriptomic map* — GSE152371

**2. Edit `config.py`**

```python
DATA_PATH       = "data/your_file.h5ad"   # path to your downloaded .h5ad file
SPECIES         = "mouse"                  # "mouse" or "human"
ASTROCYTES_ONLY = True                     # True if file contains only astrocytes
```

---

## Running the Pipeline

Run scripts in order from inside the `neuroligin_project/` directory:

```bash
python 01_preprocess.py
python 02_cluster_annotate.py    # skip if ASTROCYTES_ONLY = True
python 03_astrocyte_subcluster.py
python 04_neuroligin_expression.py
python 05_coexpression.py
```

> **Note:** Scripts 02 and 03 require manual annotation steps. After each runs, inspect the dot plot and marker CSV in `figures/` and `results/`, then fill in the cluster annotation dictionary before re-running.

---

## Outputs

### Figures
| File | Description |
|------|-------------|
| `01a_qc_before_filtering.png` | QC metrics histogram |
| `01b_pca_variance_ratio.png` | PCA elbow plot |
| `03a_umap_astrocyte_clusters.png` | Astrocyte sub-cluster UMAP |
| `03c_umap_astrocyte_subtypes.png` | Annotated astrocyte subtype UMAP |
| `04a_feature_plots_nlgn.png` | UMAP coloured by Nlgn1/2/3 expression |
| `04b_violin_nlgn_by_subtype.png` | Nlgn expression per subtype |
| `04c_heatmap_nlgn_by_subtype.png` | Mean expression heatmap |
| `05a_coexpr_dendrogram.png` | Gene co-expression dendrogram |
| `05b_coexpr_network.png` | Nlgn3 co-expression network graph |
| `05c_GO_enrichment.png` | GO/KEGG enrichment bar chart |

### Results (CSV)
- `04_mean_nlgn_expression.csv` — mean Nlgn expression per astrocyte subtype
- `04_nlgn_differential_expression.csv` — DE results with FDR correction
- `05_correlation_matrix.csv` — full gene-gene Pearson correlation matrix
- `05_gene_modules.csv` — hierarchical co-expression module assignments
- `05_Nlgn3_module_genes.csv` — genes in the Nlgn3 co-expression module
- `05_top_coexpr_Nlgn*.csv` — top 20 co-expressed genes per Nlgn isoform

---

## Project Structure

```
neuroligin_project/
├── config.py                      # All settings — edit before running
├── utils.py                       # Shared helper functions
├── requirements.txt               # Python dependencies
├── 00_download_instructions.txt   # How to get the data
├── 01_preprocess.py
├── 02_cluster_annotate.py
├── 03_astrocyte_subcluster.py
├── 04_neuroligin_expression.py
├── 05_coexpression.py
├── data/                          # .h5ad files (not tracked by git)
├── figures/                       # Output figures (not tracked by git)
└── results/                       # Output CSVs (not tracked by git)
```

---

## Authors

Sathvik Loke & Pranith Valleri
