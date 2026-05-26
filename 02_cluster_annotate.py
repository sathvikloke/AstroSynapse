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

# ── 1b. Save any pre-existing cell_type column from the source data ──
# Allen Brain / CellXGene datasets often ship with their own annotations.
# We preserve them here so we can fall back to them if marker-gene matching fails.
orig_cell_type = None
if "cell_type" in adata.obs.columns:
    orig_cell_type = adata.obs["cell_type"].copy()
    print(f"\n  Pre-existing cell_type labels found ({orig_cell_type.nunique()} types):")
    print(orig_cell_type.value_counts().to_string())

# ── 2. Build neighbour graph ──────────────────────────────────────
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

# ── 5b. Plot UMAP coloured by original cell_type if available ────
if orig_cell_type is not None:
    fig_ct, ax_ct = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color="cell_type", ax=ax_ct, show=False,
               title="Cell types (dataset annotations)")
    plt.tight_layout()
    utils.save_fig(fig_ct, os.path.join(config.FIGURES_DIR, "02a2_umap_orig_celltypes.png"))

# ── 6. Marker gene dot plot ────────────────────────────────────────
markers = config.CELL_TYPE_MARKERS[config.SPECIES]
all_marker_genes = [g for genes in markers.values() for g in genes]
present_markers, _ = utils.find_genes(adata.var_names, all_marker_genes,
                                       label="cell-type markers")

if present_markers:
    markers_present = {}
    for cell_type, gene_list in markers.items():
        found, _ = utils.find_genes(adata.var_names, gene_list, verbose=False)
        if found:
            markers_present[cell_type] = found

    dp = sc.pl.dotplot(
        adata,
        var_names=markers_present,
        groupby="leiden",
        show=False,
        return_fig=True,
        title="Cell-type markers by Leiden cluster",
    )
    dotplot_path = os.path.join(config.FIGURES_DIR, "02b_dotplot_markers.png")
    dp_fig = dp.get_axes()["mainplot_ax"].get_figure()
    utils.save_fig(dp_fig, dotplot_path)
    plt.close("all")
    print(f"  Saved: {dotplot_path}")
    print("  Use this dot plot to fill in CLUSTER_TO_CELLTYPE below.")
else:
    print("  Warning: no marker genes found in dataset — skipping dot plot.")
    if orig_cell_type is not None:
        print("  Will use pre-existing dataset cell_type labels instead.")

# ── 7. Compute per-cluster marker genes ───────────────────────────
print("\nComputing per-cluster marker genes (Wilcoxon test) ...")
# FIX: must set use_raw=False when specifying a layer.
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                        n_genes=20, layer="log_norm", use_raw=False)

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

# ── 8. Cell-type annotation ───────────────────────────────────────
# Option A: If the dataset already has cell_type labels, use them directly.
# Option B: Fill in CLUSTER_TO_CELLTYPE manually after inspecting the dot plot.

CLUSTER_TO_CELLTYPE = {
    # Cluster : Cell type
    # Fill these in after inspecting 02b_dotplot_markers.png
    # ↓ example entries (likely wrong for your data — replace them) ↓
}

if not CLUSTER_TO_CELLTYPE:
    if orig_cell_type is not None:
        # Use the pre-existing labels that came with the dataset
        print(
            "\n  Using pre-existing cell_type annotations from the dataset.\n"
            "  (To override, fill in CLUSTER_TO_CELLTYPE in this script.)"
        )
        adata.obs["cell_type"] = orig_cell_type
    else:
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
else:
    adata.obs["cell_type"] = (
        adata.obs["leiden"]
        .map(CLUSTER_TO_CELLTYPE)
        .fillna("Unknown")
        .astype("category")
    )

# ── 9. UMAP coloured by cell type ────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color="cell_type", ax=ax3, show=False,
           title="Cell types (annotation)")
plt.tight_layout()
utils.save_fig(fig3, os.path.join(config.FIGURES_DIR, "02c_umap_celltypes.png"))

print("\nCell type distribution:")
print(adata.obs["cell_type"].value_counts().to_string())

# ── 10. Save annotated full dataset ──────────────────────────────
annotated_path = os.path.join(config.DATA_DIR, "02_annotated.h5ad")
adata.write_h5ad(annotated_path)
print(f"\n  Saved: {annotated_path}")

# ── 11. Extract astrocytes ────────────────────────────────────────
# Try multiple label variants since different datasets use different naming.
cell_type_labels = adata.obs["cell_type"].unique().tolist()
print(f"\n  Cell type labels in this dataset: {cell_type_labels}")

astro_labels = [l for l in cell_type_labels if "astro" in str(l).lower()]
if not astro_labels:
    # Fallback: try common abbreviations
    astro_labels = [l for l in cell_type_labels
                    if str(l).lower() in ("astro", "astrocyte", "astrocytes",
                                          "glia", "astro_pp", "astro_fb")]

if not astro_labels:
    print(
        "\n  ERROR: Could not find an astrocyte label in this dataset.\n"
        f"  Labels found: {cell_type_labels}\n"
        "  Either fill in CLUSTER_TO_CELLTYPE with 'Astrocyte' as a value,\n"
        "  or check what label the dataset uses for astrocytes above."
    )
    sys.exit(1)

print(f"  Astrocyte label(s) matched: {astro_labels}")
astro_mask = adata.obs["cell_type"].isin(astro_labels)
n_astro = astro_mask.sum()

astro = adata[astro_mask].copy()
print(f"\n  Astrocytes extracted: {astro.n_obs:,} cells")

astro.X = astro.layers["log_norm"].copy()

astro_path = os.path.join(config.DATA_DIR, "02_astrocytes.h5ad")
astro.write_h5ad(astro_path)
print(f"  Saved: {astro_path}")

utils.print_section("Step 2 complete")
print("  Next step → run  03_astrocyte_subcluster.py")
