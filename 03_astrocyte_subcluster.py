# 03_astrocyte_subcluster.py
# ─────────────────────────────────────────────────────────────────
# Step 3: Re-cluster astrocytes, annotate subtypes, produce UMAPs.
#
# Input  (ASTROCYTES_ONLY = True)  → data/01_preprocessed.h5ad
# Input  (ASTROCYTES_ONLY = False) → data/02_astrocytes.h5ad
# Output → data/03_astrocytes.h5ad
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

sc.settings.verbosity = 2
sc.settings.figdir = config.FIGURES_DIR
utils.ensure_dirs(config.DATA_DIR, config.FIGURES_DIR, config.RESULTS_DIR)

utils.print_section("STEP 3 — Astrocyte Sub-clustering")

# ── 1. Load the correct input file ───────────────────────────────
if config.ASTROCYTES_ONLY:
    in_path = os.path.join(config.DATA_DIR, "01_preprocessed.h5ad")
    print(f"\nASTROCYTES_ONLY = True — loading: {in_path}")
    adata = sc.read_h5ad(in_path)
    # adata.X currently holds scaled values from step 1.
    # Restore log-normalised values for a clean re-analysis.
    adata.X = adata.layers["log_norm"].copy()
else:
    in_path = os.path.join(config.DATA_DIR, "02_astrocytes.h5ad")
    print(f"\nASTROCYTES_ONLY = False — loading: {in_path}")
    adata = sc.read_h5ad(in_path)
    # Step 2 already restored log_norm to .X when saving 02_astrocytes.h5ad

utils.summarise_adata(adata, label="input")

# ── 2. Re-identify HVGs on astrocytes only ────────────────────────
# HVGs from the full brain may differ from those that are variable
# specifically within the astrocyte population.
print("\nIdentifying HVGs within astrocytes ...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=config.N_TOP_HVG,
    subset=False,
    flavor="seurat",
)
print(f"  HVGs: {adata.var['highly_variable'].sum():,}")

# Save the astrocyte log-norm layer again (needed by later scripts)
# adata.X is still log-norm at this point.
adata.layers["log_norm"] = adata.X.copy()

# ── 3. Scale and PCA (astrocyte-specific) ────────────────────────
# Scaling overwrites adata.X; log_norm layer is already saved above.
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=config.N_PCS_COMPUTE, use_highly_variable=True)
print(f"  PCA complete ({config.N_PCS_COMPUTE} components).")

# ── 4. Neighbour graph, UMAP, Leiden ─────────────────────────────
print(f"  Building neighbour graph (k={config.N_NEIGHBORS}, n_pcs={config.N_PCS_USE}) ...")
sc.pp.neighbors(adata, n_neighbors=config.N_NEIGHBORS, n_pcs=config.N_PCS_USE)

print("  Computing UMAP ...")
sc.tl.umap(adata)

print(f"  Leiden clustering (resolution = {config.LEIDEN_RES_ASTROCYTES}) ...")
sc.tl.leiden(adata, resolution=config.LEIDEN_RES_ASTROCYTES,
             key_added="astrocyte_cluster")
n_clusters = adata.obs["astrocyte_cluster"].nunique()
print(f"  Astrocyte sub-clusters found: {n_clusters}")

# ── 5. Plot raw UMAP (Figure 2) ───────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color="astrocyte_cluster", ax=ax, show=False,
           title="Astrocyte sub-clusters (Leiden)")
plt.tight_layout()
utils.save_fig(fig, os.path.join(config.FIGURES_DIR, "03a_umap_astrocyte_clusters.png"))
print("  → Figure 2 saved: 03a_umap_astrocyte_clusters.png")

# ── 6. Subtype marker dot plot ────────────────────────────────────
subtype_markers = config.ASTROCYTE_SUBTYPE_MARKERS[config.SPECIES]
markers_present = {}
for subtype, gene_list in subtype_markers.items():
    found, _ = utils.find_genes(adata.var_names, gene_list, verbose=False)
    if found:
        markers_present[subtype] = found

if markers_present:
    # BUG FIX: same dotplot layout issue as 02_cluster_annotate.py.
    # Use return_fig=True and save via the DotPlot object.
    dp = sc.pl.dotplot(
        adata,
        var_names=markers_present,
        groupby="astrocyte_cluster",
        show=False,
        return_fig=True,
        title="Astrocyte subtype markers",
    )
    dotplot_path = os.path.join(config.FIGURES_DIR, "03b_dotplot_subtype_markers.png")
    # Extract the backing figure through mainplot_ax for cross-version compatibility
    dp_fig = dp.get_axes()["mainplot_ax"].get_figure()
    utils.save_fig(dp_fig, dotplot_path)
    plt.close("all")
    print(f"  Saved: {dotplot_path}  — use to assign subtypes below.")

# ── 7. Per-cluster marker genes → CSV ────────────────────────────
print("\n  Computing per-cluster marker genes (Wilcoxon) ...")
# BUG FIX: adata.X is scaled at this point (after sc.pp.scale).
# layer="log_norm" ensures DE testing uses log-normalised expression values.
sc.tl.rank_genes_groups(adata, groupby="astrocyte_cluster",
                        method="wilcoxon", n_genes=20, layer="log_norm")

rows = []
for cluster in adata.obs["astrocyte_cluster"].cat.categories:
    genes  = adata.uns["rank_genes_groups"]["names"][cluster]
    scores = adata.uns["rank_genes_groups"]["scores"][cluster]
    pvals  = adata.uns["rank_genes_groups"]["pvals_adj"][cluster]
    for g, s, p in zip(genes, scores, pvals):
        rows.append({"cluster": cluster, "gene": g,
                     "score": float(s), "pval_adj": float(p)})

marker_df = pd.DataFrame(rows)
csv_path = os.path.join(config.RESULTS_DIR, "03_astrocyte_cluster_markers.csv")
marker_df.to_csv(csv_path, index=False)
print(f"  Marker genes saved: {csv_path}")

# ── 8. MANUAL SUBTYPE ANNOTATION — edit this dict ────────────────
# Inspect 03b_dotplot_subtype_markers.png and 03_astrocyte_cluster_markers.csv
# then assign each cluster to one of the subtype labels below.
#
# Available labels (from config.ASTROCYTE_SUBTYPE_MARKERS):
#   Grey_matter, White_matter, Bergmann_glia, Reactive, Perisynaptic
#   or "Unknown" if unsure.
#
# Example:
#   ASTRO_CLUSTER_TO_SUBTYPE = {"0": "Grey_matter", "1": "Perisynaptic", ...}

ASTRO_CLUSTER_TO_SUBTYPE = {
    # cluster : subtype
    # ↓ Fill in after inspecting the dot plot and marker CSV ↓
}

if not ASTRO_CLUSTER_TO_SUBTYPE:
    print(
        "\n  *** ACTION REQUIRED ***\n"
        "  Fill in ASTRO_CLUSTER_TO_SUBTYPE after inspecting:\n"
        "    03b_dotplot_subtype_markers.png\n"
        "    results/03_astrocyte_cluster_markers.csv\n"
        "  Continuing with numeric cluster IDs as placeholders.\n"
    )
    ASTRO_CLUSTER_TO_SUBTYPE = {
        c: f"Cluster_{c}"
        for c in adata.obs["astrocyte_cluster"].cat.categories
    }

adata.obs["astrocyte_subtype"] = (
    adata.obs["astrocyte_cluster"]
    .map(ASTRO_CLUSTER_TO_SUBTYPE)
    .fillna("Unknown")
    .astype("category")
)

# ── 9. Annotated UMAP (Figure 2 final) ───────────────────────────
fig3, ax3 = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color="astrocyte_subtype", ax=ax3, show=False,
           title="Astrocyte subtypes")
plt.tight_layout()
utils.save_fig(fig3, os.path.join(config.FIGURES_DIR, "03c_umap_astrocyte_subtypes.png"))
print("  → Annotated Figure 2 saved: 03c_umap_astrocyte_subtypes.png")

# ── 10. Optional: UMAP by brain region ───────────────────────────
# If the dataset contains a region column, plot UMAP coloured by region.
region_cols = [c for c in adata.obs.columns
               if any(k in c.lower() for k in ["region", "area", "tissue", "brain"])]
if region_cols:
    region_col = region_cols[0]
    print(f"\n  Region column detected: '{region_col}' — plotting regional UMAP ...")
    fig4, ax4 = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color=region_col, ax=ax4, show=False,
               title=f"Astrocytes coloured by {region_col}")
    plt.tight_layout()
    utils.save_fig(fig4, os.path.join(config.FIGURES_DIR, "03d_umap_astrocyte_region.png"))
else:
    print("\n  No brain region column found in metadata — skipping regional UMAP.")
    print("  (If your dataset has region metadata under a different column name,")
    print("   add it to the region_cols search list above.)")

# ── 11. Save ──────────────────────────────────────────────────────
out_path = os.path.join(config.DATA_DIR, "03_astrocytes.h5ad")
adata.write_h5ad(out_path)

utils.print_section("Step 3 complete")
print(f"  Output : {out_path}")
utils.summarise_adata(adata)
print("\n  Next step → run  04_neuroligin_expression.py")
