# 04_neuroligin_expression.py
# ─────────────────────────────────────────────────────────────────
# Step 4: Neuroligin isoform expression analysis across astrocyte
#         subtypes and brain regions.
#
# Produces Figures 3, 4, and 5 from the project plan:
#   Fig 3 — Feature plots (UMAP coloured by Nlgn1/2/3 expression)
#   Fig 4 — Violin plots (expression per astrocyte subtype)
#   Fig 5 — Heatmap (mean expression by subtype × region)
#
# Input  : data/03_astrocytes.h5ad
# Output : figures/04_*.png,  results/04_*.csv
# ─────────────────────────────────────────────────────────────────

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
os.chdir(_SCRIPT_DIR)
import config
import utils

sc.settings.verbosity = 1
sc.settings.figdir = config.FIGURES_DIR
utils.ensure_dirs(config.DATA_DIR, config.FIGURES_DIR, config.RESULTS_DIR)

utils.print_section("STEP 4 — Neuroligin Expression Analysis")

# ── 1. Load data ──────────────────────────────────────────────────
in_path = os.path.join(config.DATA_DIR, "03_astrocytes.h5ad")
print(f"\nLoading: {in_path}")
adata = sc.read_h5ad(in_path)
utils.summarise_adata(adata)

# ── 2. Locate neuroligin genes ────────────────────────────────────
target_genes = config.NLGN_GENES[config.SPECIES]
nlgn_found, nlgn_missing = utils.find_genes(
    adata.var_names, target_genes, label="Nlgn genes"
)

if not nlgn_found:
    print(
        "\n  ERROR: None of the neuroligin genes were found in this dataset.\n"
        f"  Looked for: {target_genes}\n"
        "  Check that SPECIES is set correctly in config.py, or that the\n"
        "  gene names match the naming convention used in your dataset."
    )
    sys.exit(1)

print(f"  Neuroligin genes found: {nlgn_found}")
if nlgn_missing:
    print(f"  Not found (dropout or absent): {nlgn_missing}")

# ── 3. Determine the groupby column ──────────────────────────────
# Use the annotated subtype column if available; fall back to raw cluster IDs.
if "astrocyte_subtype" in adata.obs.columns:
    groupby_col = "astrocyte_subtype"
elif "astrocyte_cluster" in adata.obs.columns:
    groupby_col = "astrocyte_cluster"
else:
    raise KeyError(
        "Neither 'astrocyte_subtype' nor 'astrocyte_cluster' found in adata.obs.\n"
        "Make sure 03_astrocyte_subcluster.py ran successfully."
    )
print(f"  Grouping by: '{groupby_col}'")

# ── 4. Restore log-normalised values for expression plots ─────────
# adata.X currently holds scaled values (from step 3 PCA prep).
# Feature plots and violin plots should use log-normalised, unscaled values
# so that expression levels are biologically interpretable.
if "log_norm" in adata.layers:
    adata.X = adata.layers["log_norm"].copy()
    print("  Restored log-normalised expression to adata.X for plotting.")
else:
    print("  Warning: 'log_norm' layer not found. Using adata.X as-is.")

# ── 5. FIGURE 3 — Feature plots ───────────────────────────────────
# UMAP coloured by expression level of each Nlgn isoform.
n_genes = len(nlgn_found)
fig3, axes = plt.subplots(1, n_genes, figsize=(5 * n_genes, 4.5))
if n_genes == 1:
    axes = [axes]   # ensure axes is always iterable

for ax, gene in zip(axes, nlgn_found):
    sc.pl.umap(
        adata,
        color=gene,
        ax=ax,
        show=False,
        title=gene,
        color_map="RdYlBu_r",
        vmin=0,
    )

fig3.suptitle("Figure 3 — Neuroligin isoform expression on astrocyte UMAP",
              fontsize=12, y=1.02)
plt.tight_layout()
fig3_path = os.path.join(config.FIGURES_DIR, "04a_feature_plots_nlgn.png")
utils.save_fig(fig3, fig3_path)
print(f"\n  → Figure 3 saved: {fig3_path}")

# ── 6. FIGURE 4 — Violin plots ────────────────────────────────────
# One panel per Nlgn isoform, showing expression across subtypes.
fig4, axes4 = plt.subplots(1, n_genes, figsize=(6 * n_genes, 5), sharey=False)
if n_genes == 1:
    axes4 = [axes4]

for ax, gene in zip(axes4, nlgn_found):
    sc.pl.violin(
        adata,
        keys=gene,
        groupby=groupby_col,
        ax=ax,
        show=False,
        rotation=45,
        stripplot=False,    # hide individual points for clarity on large datasets
        inner="box",        # show box-plot summary inside violins
    )
    ax.set_title(gene)
    ax.set_xlabel("")

fig4.suptitle("Figure 4 — Neuroligin expression by astrocyte subtype",
              fontsize=12, y=1.02)
plt.tight_layout()
fig4_path = os.path.join(config.FIGURES_DIR, "04b_violin_nlgn_by_subtype.png")
utils.save_fig(fig4, fig4_path)
print(f"  → Figure 4 saved: {fig4_path}")

# ── 7. FIGURE 5 — Heatmap of mean expression ─────────────────────
# Mean log-normalised expression of each isoform across subtypes.
# If region metadata is present, also produce a region × subtype heatmap.

# Compute mean expression per subtype
X_log = utils.get_layer(adata, "log_norm")
groups = adata.obs[groupby_col].values
unique_groups = sorted(adata.obs[groupby_col].unique())

# BUG FIX: precompute gene indices once outside the loop.
# Previously list(adata.var_names).index(g) was called inside the loop,
# rebuilding the full list on every iteration — O(n_groups × n_genes × n_nlgn).
gene_name_list = list(adata.var_names)
nlgn_indices   = [gene_name_list.index(g) for g in nlgn_found]

mean_expr = pd.DataFrame(
    index=unique_groups,
    columns=nlgn_found,
    dtype=float,
)
for grp in unique_groups:
    mask = groups == grp
    mean_expr.loc[grp] = X_log[mask][:, nlgn_indices].mean(axis=0)

# Save mean expression table
mean_csv = os.path.join(config.RESULTS_DIR, "04_mean_nlgn_expression.csv")
mean_expr.to_csv(mean_csv)
print(f"  Mean expression saved: {mean_csv}")

# Plot heatmap
fig5, ax5 = plt.subplots(figsize=(max(4, n_genes * 1.5),
                                   max(4, len(unique_groups) * 0.6)))
sns.heatmap(
    mean_expr.astype(float),
    ax=ax5,
    cmap="YlOrRd",
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "Mean log-normalised expression"},
)
ax5.set_title("Figure 5 — Mean Neuroligin expression per astrocyte subtype")
ax5.set_xlabel("Isoform")
ax5.set_ylabel("Astrocyte subtype")
plt.tight_layout()
fig5_path = os.path.join(config.FIGURES_DIR, "04c_heatmap_nlgn_by_subtype.png")
utils.save_fig(fig5, fig5_path)
print(f"  → Figure 5 saved: {fig5_path}")

# Optional: heatmap split by brain region
region_cols = [c for c in adata.obs.columns
               if any(k in c.lower() for k in ["region", "area", "tissue", "brain"])]
if region_cols:
    region_col = region_cols[0]
    print(f"\n  Region column found ('{region_col}') — generating region heatmap ...")
    regions = sorted(adata.obs[region_col].unique())

    # Build multi-index mean expression: rows = (subtype, region)
    rows_region = []
    # BUG FIX: nlgn_indices already computed above — reuse it here instead of
    # rebuilding list(adata.var_names) inside the double loop.
    for grp in unique_groups:
        for reg in regions:
            mask = (adata.obs[groupby_col] == grp) & (adata.obs[region_col] == reg)
            n_cells = mask.sum()
            if n_cells < 5:          # skip groups with < 5 cells (unreliable mean)
                continue
            vals = X_log[mask][:, nlgn_indices].mean(axis=0)
            row = {"subtype": grp, "region": reg, "n_cells": n_cells}
            for gene, val in zip(nlgn_found, vals):
                row[gene] = float(val)
            rows_region.append(row)

    region_df = pd.DataFrame(rows_region)
    region_csv = os.path.join(config.RESULTS_DIR, "04_mean_nlgn_by_region.csv")
    region_df.to_csv(region_csv, index=False)
    print(f"  Regional expression saved: {region_csv}")

    # Pivot for heatmap (use first Nlgn gene found for the region plot)
    if nlgn_found:
        pivot_gene = nlgn_found[-1]    # typically Nlgn3
        pivot = region_df.pivot(index="subtype", columns="region", values=pivot_gene)
        fig6, ax6 = plt.subplots(figsize=(max(5, len(regions) * 0.9),
                                           max(4, len(unique_groups) * 0.6)))
        sns.heatmap(pivot, ax=ax6, cmap="YlOrRd", annot=True, fmt=".2f",
                    linewidths=0.5,
                    cbar_kws={"label": "Mean log-norm expression"})
        ax6.set_title(f"{pivot_gene} mean expression — subtype × brain region")
        plt.tight_layout()
        fig6_path = os.path.join(config.FIGURES_DIR, "04d_heatmap_nlgn3_by_region.png")
        utils.save_fig(fig6, fig6_path)
        print(f"  Regional heatmap saved: {fig6_path}")
else:
    print("\n  No brain region column found — skipping regional heatmap.")

# ── 8. Differential expression of Nlgn genes across subtypes ──────
# Tests whether each Nlgn gene is significantly DE between any two subtypes.
print("\n  Running differential expression test for Nlgn genes ...")

de_rows = []
subtypes = adata.obs[groupby_col].cat.categories.tolist()

# Wilcoxon rank-sum: each subtype vs. all others (one-vs-rest)
X_log = utils.get_layer(adata, "log_norm")
gene_names = list(adata.var_names)

for gene in nlgn_found:
    gene_idx = gene_names.index(gene)
    expr = X_log[:, gene_idx]

    for subtype in subtypes:
        mask_group = (adata.obs[groupby_col] == subtype).values
        mask_rest  = ~mask_group

        group_expr = expr[mask_group]
        rest_expr  = expr[mask_rest]

        # Need at least 3 cells in each group for a meaningful test
        if len(group_expr) < 3 or len(rest_expr) < 3:
            continue

        stat, pval = mannwhitneyu(group_expr, rest_expr, alternative="two-sided")
        mean_group = float(group_expr.mean())
        mean_rest  = float(rest_expr.mean())
        log2fc     = float(np.log2((mean_group + 1e-9) / (mean_rest + 1e-9)))

        de_rows.append({
            "gene":       gene,
            "subtype":    subtype,
            "mean_group": round(mean_group, 4),
            "mean_rest":  round(mean_rest,  4),
            "log2fc":     round(log2fc, 4),
            "pval":       float(pval),
        })

if de_rows:
    de_df = pd.DataFrame(de_rows)
    # Multiple-testing correction (Benjamini-Hochberg) across all tests
    _, pvals_adj, _, _ = multipletests(de_df["pval"].values, method="fdr_bh")
    de_df["pval_adj"] = pvals_adj
    de_df = de_df.sort_values("pval_adj")

    de_csv = os.path.join(config.RESULTS_DIR, "04_nlgn_differential_expression.csv")
    de_df.to_csv(de_csv, index=False)
    print(f"  DE results saved: {de_csv}")

    sig = de_df[de_df["pval_adj"] < 0.05]
    print(f"  Significant results (FDR < 0.05): {len(sig)}")
    if not sig.empty:
        print(sig[["gene", "subtype", "log2fc", "pval_adj"]].to_string(index=False))
else:
    print("  Warning: no DE tests could be run (insufficient cells per group).")

utils.print_section("Step 4 complete")
print("  Figures: figures/04a_feature_plots_nlgn.png")
print("           figures/04b_violin_nlgn_by_subtype.png")
print("           figures/04c_heatmap_nlgn_by_subtype.png")
print("\n  Next step → run  05_coexpression.py")
