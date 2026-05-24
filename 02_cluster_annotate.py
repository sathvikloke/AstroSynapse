# 02_cluster_annotate.py
# ─────────────────────────────────────────────────────────────────
# Step 2 (only needed when ASTROCYTES_ONLY = False):
#   Build neighbour graph → UMAP → Leiden clustering →
#   cell-type annotation → extract astrocytes.
#
# Input  : data/01_preprocessed.h5ad
# Output : data/02_annotated.h5ad       (all cell types, labelled)
#          data/02_astrocytes.h5ad      (astrocytes only, for step 3)
# ─────────────────────────────────────────────────────────────────

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
os.chdir(_SCRIPT_DIR)
import config
import utils

# Guard: skip if not needed
if config.ASTROCYTES_ONLY:
    print("ASTROCYTES_ONLY = True — this script is not needed.")
    print("Proceed directly to 03_astrocyte_subcluster.py.")
    sys.exit(0)

sc.settings.verbosity = 2
sc.settings.figdir = config.FIGURES_DIR
utils.ensure_dirs(config.DATA_DIR, config.FIGURES_DIR, config.RESULTS_DIR)

utils.print_section("STEP 2 — Clustering & Cell-Type Annotation")

# ── 1. Load preprocessed data ─────────────────────────────────────
in_path = os.path.join(config.DATA_DIR, "01_preprocessed.h5ad")
print(f"\nLoading: {in_path}")
adata = sc.read_h5ad(in_path)
utils.summarise_adata(adata)

# ── 2. Build neighbour graph ──────────────────────────────────────
# Uses the PCA embedding (already computed in step 1).
print(f"\nBuilding k-NN graph  (k={config.N_NEIGHBORS}, n_pcs={config.N_PCS_USE}) ...")
sc.pp.neighbors(adata, n_neighbors=config.N_NEIGHBORS, n_pcs=config.N_PCS_USE)

# ── 3. UMAP ───────────────────────────────────────────────────────
print("Computing UMAP ...")
sc.tl.umap(adata)

# ── 4. Leiden clustering ──────────────────────────────────────────
print(f"Leiden clustering (resolution = {config.LEIDEN_RES_FULL}) ...")
sc.tl.leiden(adata, resolution=config.LEIDEN_RES_FULL, key_added="leiden")
n_clusters = adata.obs["leiden"].nunique()
print(f"  Clusters found: {n_clusters}")

# ── 5. Plot UMAP coloured by Leiden cluster ───────────────────────
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color="leiden", ax=ax, show=False, title="Leiden clusters (all cells)")
plt.tight_layout()
utils.save_fig(fig, os.path.join(config.FIGURES_DIR, "02a_umap_leiden.png"))

# ── 6. Marker gene dot plot ────────────────────────────────────────
# Shows expression of known cell-type markers across clusters to guide annotation.
markers = config.CELL_TYPE_MARKERS[config.SPECIES]
# Flatten to a list and keep only genes present in the dataset
all_marker_genes = [g for genes in markers.values() for g in genes]
present_markers, _ = utils.find_genes(adata.var_names, all_marker_genes,
                                       label="cell-type markers")

if present_markers:
    # Build ordered dict: cell_type → [genes present in dataset]
    markers_present = {}
    for cell_type, gene_list in markers.items():
        found, _ = utils.find_genes(adata.var_names, gene_list, verbose=False)
        if found:
            markers_present[cell_type] = found

    # BUG FIX: sc.pl.dotplot manages its own multi-axes layout (dots + colorbar +
    # size legend).  Passing an external ax= from plt.subplots() causes a layout
    # conflict and the figure may not save correctly.  Use return_fig=True instead
    # to get the DotPlot object and save it directly.
    dp = sc.pl.dotplot(
        adata,
        var_names=markers_present,
        groupby="leiden",
        show=False,
        return_fig=True,
        title="Cell-type markers by Leiden cluster",
    )
    dotplot_path = os.path.join(config.FIGURES_DIR, "02b_dotplot_markers.png")
    # Retrieve the underlying matplotlib Figure via the mainplot axes so that
    # utils.save_fig() (which calls savefig with bbox_inches="tight") works
    # correctly even across different scanpy versions.
    dp_fig = dp.get_axes()["mainplot_ax"].get_figure()
    utils.save_fig(dp_fig, dotplot_path)
    plt.close("all")
    print(f"  Saved: {dotplot_path}")
    print("  Use this dot plot to fill in CLUSTER_TO_CELLTYPE below.")
else:
    print("  Warning: no marker genes found in dataset. Check SPECIES in config.py.")

# ── 7. Compute per-cluster marker genes ───────────────────────────
# This finds the top differentially expressed genes for each cluster
# and saves them to a CSV so you can look them up.
print("\nComputing per-cluster marker genes (Wilcoxon test) ...")
# BUG FIX: default uses adata.X which is scaled after PCA prep.
# layer="log_norm" ensures the test runs on biologically meaningful values.
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                        n_genes=20, layer="log_norm")

# Save top markers to CSV
marker_df_rows = []
for cluster in adata.obs["leiden"].cat.categories:
    genes  = adata.uns["rank_genes_groups"]["names"][cluster]
    scores = adata.uns["rank_genes_groups"]["scores"][cluster]
    pvals  = adata.uns["rank_genes_groups"]["pvals_adj"][cluster]
    for gene, score, pval in zip(genes, scores, pvals):
        marker_df_rows.append({
            "cluster": cluster,
            "gene":    gene,
            "score":   float(score),
            "pval_adj": float(pval),
        })

marker_df = pd.DataFrame(marker_df_rows)
marker_csv = os.path.join(config.RESULTS_DIR, "02_cluster_markers.csv")
marker_df.to_csv(marker_csv, index=False)
print(f"  Top marker genes saved to: {marker_csv}")

# ── 8. MANUAL ANNOTATION — edit this dict ─────────────────────────
# Look at the dot plot (02b_dotplot_markers.png) and the marker CSV.
# Replace the values below with the correct cell type for each cluster.
# Example:  "0": "Astrocyte",  "1": "Neuron",  "2": "Oligodendrocyte"
#
# Available labels (from config.CELL_TYPE_MARKERS):
#   Astrocyte, Neuron, Oligodendrocyte, OPC, Microglia, Endothelial, Pericyte
#   or "Unknown" for clusters you cannot confidently assign.

CLUSTER_TO_CELLTYPE = {
    # Cluster : Cell type
    # Fill these in after inspecting 02b_dotplot_markers.png
    # ↓ example entries (likely wrong for your data — replace them) ↓
}

# If the annotation dict is empty, assign all clusters as "Unassigned"
# so the rest of the pipeline can still run.
if not CLUSTER_TO_CELLTYPE:
    print(
        "\n  *** ACTION REQUIRED ***\n"
        "  Open 02b_dotplot_markers.png and 02_cluster_markers.csv,\n"
        "  then fill in the CLUSTER_TO_CELLTYPE dict in this script.\n"
        "  Re-run to produce the annotated UMAP and extract astrocytes.\n"
        "  Continuing with all clusters labelled 'Unassigned' for now.\n"
    )
    CLUSTER_TO_CELLTYPE = {c: "Unassigned"
                           for c in adata.obs["leiden"].cat.categories}

adata.obs["cell_type"] = (
    adata.obs["leiden"]
    .map(CLUSTER_TO_CELLTYPE)
    .fillna("Unknown")
    .astype("category")
)

# ── 9. UMAP coloured by cell type ────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color="cell_type", ax=ax3, show=False,
           title="Cell types (manual annotation)")
plt.tight_layout()
utils.save_fig(fig3, os.path.join(config.FIGURES_DIR, "02c_umap_celltypes.png"))

print("\nCell type distribution:")
print(adata.obs["cell_type"].value_counts().to_string())

# ── 10. Save annotated full dataset ──────────────────────────────
annotated_path = os.path.join(config.DATA_DIR, "02_annotated.h5ad")
adata.write_h5ad(annotated_path)
print(f"\n  Saved: {annotated_path}")

# ── 11. Extract astrocytes ────────────────────────────────────────
astro_mask = adata.obs["cell_type"] == "Astrocyte"
n_astro = astro_mask.sum()

if n_astro == 0:
    print(
        "\n  ERROR: No cells labelled 'Astrocyte' found.\n"
        "  Make sure at least one cluster is assigned 'Astrocyte'\n"
        "  in the CLUSTER_TO_CELLTYPE dict above."
    )
    sys.exit(1)

astro = adata[astro_mask].copy()
print(f"\n  Astrocytes extracted: {astro.n_obs:,} cells")

# Restore log-normalised values from the log_norm layer so step 3
# can re-run HVG selection on clean (unscaled) astrocyte data.
astro.X = astro.layers["log_norm"].copy()

astro_path = os.path.join(config.DATA_DIR, "02_astrocytes.h5ad")
astro.write_h5ad(astro_path)
print(f"  Saved: {astro_path}")

utils.print_section("Step 2 complete")
print("  Next step → run  03_astrocyte_subcluster.py")
