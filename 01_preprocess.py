# 01_preprocess.py
# ─────────────────────────────────────────────────────────────────
# Step 1: Quality control, normalisation, HVG selection, and PCA.
#
# Input  : config.DATA_PATH  (raw .h5ad file)
# Output : data/01_preprocessed.h5ad
#
# Run this script first regardless of whether ASTROCYTES_ONLY is
# True or False.
# ─────────────────────────────────────────────────────────────────

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend — safe on any machine
import matplotlib.pyplot as plt
import scanpy as sc

# Make sure config and utils are importable when running from any directory,
# and that relative paths in config.py resolve correctly regardless of cwd.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
os.chdir(_SCRIPT_DIR)
import config
import utils

# ── Setup ─────────────────────────────────────────────────────────
sc.settings.verbosity = 2
sc.settings.figdir = config.FIGURES_DIR
utils.ensure_dirs(config.DATA_DIR, config.FIGURES_DIR, config.RESULTS_DIR)

utils.print_section("STEP 1 — Quality Control & Preprocessing")

# ── 1. Load data ──────────────────────────────────────────────────
print(f"\nLoading: {config.DATA_PATH}")
if not os.path.exists(config.DATA_PATH):
    raise FileNotFoundError(
        f"Data file not found: {config.DATA_PATH}\n"
        "Please download the dataset and update DATA_PATH in config.py.\n"
        "See 00_download_instructions.txt for details."
    )

adata = sc.read_h5ad(config.DATA_PATH)
print(f"  Loaded: {adata.n_obs:,} cells  x  {adata.n_vars:,} genes")

# Ensure string indices and unique gene names
adata.var_names = adata.var_names.astype(str)
adata.obs_names = adata.obs_names.astype(str)
adata.var_names_make_unique()

# ── 2. Annotate mitochondrial genes ───────────────────────────────
# Mouse gene names start with lowercase 'mt-'; human with 'MT-'
mt_prefix = "mt-" if config.SPECIES == "mouse" else "MT-"
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
print(f"  Mitochondrial genes found: {adata.var['mt'].sum()}")

# ── 3. Compute QC metrics ─────────────────────────────────────────
# Adds to adata.obs:
#   n_genes_by_counts, total_counts, pct_counts_mt
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

# ── 4. QC plots (BEFORE filtering) ────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].hist(adata.obs["n_genes_by_counts"], bins=100,
             color="steelblue", edgecolor="none")
axes[0].axvline(config.MIN_GENES_PER_CELL, color="red",  linestyle="--",
                label=f"min = {config.MIN_GENES_PER_CELL}")
axes[0].axvline(config.MAX_GENES_PER_CELL, color="orange", linestyle="--",
                label=f"max = {config.MAX_GENES_PER_CELL}")
axes[0].set_xlabel("Genes per cell")
axes[0].set_ylabel("Number of cells")
axes[0].set_title("Genes per cell")
axes[0].legend(fontsize=8)

axes[1].hist(adata.obs["total_counts"], bins=100,
             color="steelblue", edgecolor="none")
axes[1].set_xlabel("Total counts (UMIs)")
axes[1].set_title("Counts per cell")

axes[2].hist(adata.obs["pct_counts_mt"], bins=100,
             color="steelblue", edgecolor="none")
axes[2].axvline(config.MAX_MITO_PERCENT, color="red", linestyle="--",
                label=f"max = {config.MAX_MITO_PERCENT} %")
axes[2].set_xlabel("% Mitochondrial counts")
axes[2].set_title("Mitochondrial %")
axes[2].legend(fontsize=8)

plt.suptitle("QC metrics — before filtering", fontsize=12, y=1.02)
plt.tight_layout()
qc_fig_path = os.path.join(config.FIGURES_DIR, "01a_qc_before_filtering.png")
utils.save_fig(fig, qc_fig_path)

# ── 5. Filter cells and genes ─────────────────────────────────────
n_cells_before = adata.n_obs
n_genes_before = adata.n_vars

sc.pp.filter_cells(adata, min_genes=config.MIN_GENES_PER_CELL)
sc.pp.filter_cells(adata, max_genes=config.MAX_GENES_PER_CELL)

# Filter on mitochondrial percentage (boolean indexing, then copy to avoid
# storing a view — views can cause issues in subsequent in-place operations)
adata = adata[adata.obs["pct_counts_mt"] < config.MAX_MITO_PERCENT].copy()

sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS_PER_GENE)

print(f"  Cells  : {n_cells_before:,}  →  {adata.n_obs:,}  "
      f"(removed {n_cells_before - adata.n_obs:,})")
print(f"  Genes  : {n_genes_before:,}  →  {adata.n_vars:,}  "
      f"(removed {n_genes_before - adata.n_vars:,})")

# ── 6. Store raw counts ───────────────────────────────────────────
# Preserve integer counts in a layer before any normalisation.
# Useful if you later want to use count-based models (e.g. scVI).
adata.layers["counts"] = adata.X.copy()

# ── 7. Normalise and log-transform ────────────────────────────────
# Normalise each cell to 10,000 total counts (CPM-like), then log1p.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save the log-normalised values BEFORE scaling.
# This layer is used downstream for co-expression analysis (not scaled data).
adata.layers["log_norm"] = adata.X.copy()
print("  Normalised (CPM) and log1p-transformed.")

# ── 8. Identify highly variable genes (HVGs) ──────────────────────
# HVGs capture the most informative genes for dimensionality reduction.
# subset=False marks genes but keeps all in the object.
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=config.N_TOP_HVG,
    subset=False,
    flavor="seurat",
)
n_hvg = adata.var["highly_variable"].sum()
print(f"  Highly variable genes identified: {n_hvg:,}")

# ── 9. Scale for PCA ──────────────────────────────────────────────
# Scaling centres each gene to mean=0, std=1. This overwrites adata.X,
# but the original log-normalised values are safely stored in adata.layers['log_norm'].
sc.pp.scale(adata, max_value=10)

# ── 10. PCA ───────────────────────────────────────────────────────
# use_highly_variable=True (default when the column exists) restricts
# PCA to the HVGs, which reduces noise from lowly-expressed genes.
sc.tl.pca(adata, n_comps=config.N_PCS_COMPUTE, use_highly_variable=True)
print(f"  PCA computed ({config.N_PCS_COMPUTE} components).")

# Elbow plot — helps decide how many PCs to use for the neighbour graph
fig2, ax2 = plt.subplots(figsize=(7, 4))
variance_ratio = adata.uns["pca"]["variance_ratio"]
ax2.plot(range(1, len(variance_ratio) + 1), variance_ratio, "o-", markersize=4)
ax2.axvline(config.N_PCS_USE, color="red", linestyle="--",
            label=f"N_PCS_USE = {config.N_PCS_USE}")
ax2.set_xlabel("Principal component")
ax2.set_ylabel("Variance explained (ratio)")
ax2.set_title("PCA — variance ratio (elbow plot)")
ax2.legend()
plt.tight_layout()
utils.save_fig(fig2, os.path.join(config.FIGURES_DIR, "01b_pca_variance_ratio.png"))

# ── 11. Save ──────────────────────────────────────────────────────
out_path = os.path.join(config.DATA_DIR, "01_preprocessed.h5ad")
adata.write_h5ad(out_path)

utils.print_section("Step 1 complete")
print(f"  Output : {out_path}")
utils.summarise_adata(adata)
print()
if config.ASTROCYTES_ONLY:
    print("  ASTROCYTES_ONLY = True")
    print("  Next step → run  03_astrocyte_subcluster.py")
else:
    print("  ASTROCYTES_ONLY = False")
    print("  Next step → run  02_cluster_annotate.py")
