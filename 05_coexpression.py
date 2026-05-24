# 05_coexpression.py
# ─────────────────────────────────────────────────────────────────
# Step 5: Co-expression network analysis centred on Nlgn3.
#
# Approach:
#   1. Build a Pearson correlation matrix across the top HVGs in
#      astrocytes using log-normalised expression values.
#   2. Cluster genes hierarchically into co-expression modules.
#   3. Identify the module containing Nlgn3.
#   4. Run GO enrichment on that module (requires internet access).
#   5. Visualise the module as a network graph.
#
# Produces Figures 6, 7, and 8 from the project plan.
#
# Input  : data/03_astrocytes.h5ad
# Output : figures/05_*.png,  results/05_*.csv
# ─────────────────────────────────────────────────────────────────

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import networkx as nx
import scanpy as sc
# pearsonr removed — correlation is computed with np.corrcoef, not pearsonr
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from matplotlib.patches import Patch

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
os.chdir(_SCRIPT_DIR)
import config
import utils

utils.ensure_dirs(config.DATA_DIR, config.FIGURES_DIR, config.RESULTS_DIR)

utils.print_section("STEP 5 — Co-Expression Network Analysis")

# ── 1. Load data ──────────────────────────────────────────────────
sc.settings.verbosity = 1

in_path = os.path.join(config.DATA_DIR, "03_astrocytes.h5ad")
print(f"\nLoading: {in_path}")
adata = sc.read_h5ad(in_path)
utils.summarise_adata(adata)

# ── 2. Select genes for correlation matrix ────────────────────────
# Always include Nlgn genes regardless of expression level.
nlgn_genes = config.NLGN_GENES[config.SPECIES]
nlgn_found, _ = utils.find_genes(adata.var_names, nlgn_genes, label="Nlgn")

if not nlgn_found:
    print("ERROR: No neuroligin genes found. Cannot run co-expression analysis.")
    sys.exit(1)

# Filter genes to those expressed in >= MIN_EXPR_FRACTION of astrocytes.
# This removes very sparse genes whose correlation estimates are unreliable.
X_log = utils.get_layer(adata, "log_norm")          # shape: (cells, genes)
gene_names = np.array(list(adata.var_names))

frac_expressed = (X_log > 0).mean(axis=0)           # fraction of cells with > 0
expr_mask = frac_expressed >= config.MIN_EXPR_FRACTION

# Force-include Nlgn genes even if they fall below the threshold
for gene in nlgn_found:
    idx = np.where(gene_names == gene)[0]
    if len(idx):
        expr_mask[idx[0]] = True

expressed_genes = gene_names[expr_mask]
print(f"\n  Genes passing expression filter (>= {config.MIN_EXPR_FRACTION*100:.0f}% cells): "
      f"{expr_mask.sum():,}")

# Keep top N_TOP_COEXPR_GENES HVGs + Nlgn genes for the correlation matrix.
# Limiting gene count keeps the computation manageable on a laptop.
if "highly_variable" in adata.var.columns:
    hvg_mask  = adata.var["highly_variable"].values
    combined  = expr_mask & hvg_mask
    hvg_genes = gene_names[combined]

    # If fewer than 50 genes survive the HVG+expressed filter, fall back to
    # all expressed genes — but still cap at N_TOP_COEXPR_GENES so the
    # correlation matrix stays computationally tractable.
    if len(hvg_genes) < 50:
        hvg_genes = expressed_genes[:config.N_TOP_COEXPR_GENES]
        print(f"  Fewer than 50 HVG+expressed genes — using top "
              f"{len(hvg_genes):,} expressed genes.")
    elif len(hvg_genes) > config.N_TOP_COEXPR_GENES:
        # Rank by variance and take top N
        gene_var   = X_log[:, combined].var(axis=0)
        top_idx    = np.argsort(gene_var)[::-1][:config.N_TOP_COEXPR_GENES]
        hvg_genes  = hvg_genes[top_idx]
else:
    hvg_genes = expressed_genes[:config.N_TOP_COEXPR_GENES]

# Always include Nlgn genes
for gene in nlgn_found:
    if gene not in hvg_genes:
        hvg_genes = np.append(hvg_genes, gene)

print(f"  Genes in correlation matrix: {len(hvg_genes):,}")

# Extract expression matrix for selected genes
gene_idx_map = {g: i for i, g in enumerate(gene_names)}
sel_indices  = np.array([gene_idx_map[g] for g in hvg_genes])
X_sel        = X_log[:, sel_indices]          # shape: (cells, selected_genes)

# ── 3. Compute Pearson correlation matrix ─────────────────────────
print("\n  Computing Pearson correlation matrix ...")
# np.corrcoef takes rows as variables → transpose so genes are rows
corr_matrix = np.corrcoef(X_sel.T)            # shape: (genes, genes)

# Replace any NaN values (can arise from constant-expression genes)
corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
# Set diagonal to 1.0 explicitly
np.fill_diagonal(corr_matrix, 1.0)

# Save correlation matrix
corr_df = pd.DataFrame(corr_matrix, index=hvg_genes, columns=hvg_genes)
corr_csv = os.path.join(config.RESULTS_DIR, "05_correlation_matrix.csv")
corr_df.to_csv(corr_csv)
print(f"  Correlation matrix saved: {corr_csv}")

# ── 4. Hierarchical clustering of genes (Figure 6) ────────────────
print("\n  Running hierarchical clustering on genes ...")
# Distance = 1 - |r|  (so highly correlated genes, positive or negative, cluster together)
dist_matrix = 1.0 - np.abs(corr_matrix)
np.fill_diagonal(dist_matrix, 0.0)             # genes perfectly correlated with themselves
dist_condensed = squareform(dist_matrix, checks=False)

linkage_matrix = linkage(dist_condensed, method="average")

# Cut the dendrogram at distance = 0.7  (i.e., |r| > 0.3 within clusters)
CUT_DISTANCE = 1.0 - config.CORR_THRESHOLD
gene_cluster_labels = fcluster(linkage_matrix, t=CUT_DISTANCE, criterion="distance")
n_modules = gene_cluster_labels.max()
print(f"  Gene modules identified: {n_modules}  (cut at distance {CUT_DISTANCE:.2f})")

# Plot dendrogram (Figure 6)
fig6, ax6 = plt.subplots(figsize=(14, 5))
dendrogram(
    linkage_matrix,
    labels=hvg_genes,
    ax=ax6,
    leaf_font_size=4,
    color_threshold=CUT_DISTANCE,
    above_threshold_color="grey",
)
ax6.axhline(y=CUT_DISTANCE, color="red", linestyle="--",
            label=f"Cut at distance {CUT_DISTANCE:.2f}  (|r| > {config.CORR_THRESHOLD})")
ax6.set_title("Figure 6 — Gene co-expression dendrogram (astrocytes)")
ax6.set_xlabel("Gene")
ax6.set_ylabel("Distance  (1 − |r|)")
ax6.legend(fontsize=9)
plt.tight_layout()
fig6_path = os.path.join(config.FIGURES_DIR, "05a_coexpr_dendrogram.png")
utils.save_fig(fig6, fig6_path)
print(f"  → Figure 6 saved: {fig6_path}")

# ── 5. Identify Nlgn3's co-expression module ──────────────────────
# Determine which module contains each Nlgn gene
gene_to_module = {gene: int(lbl)
                  for gene, lbl in zip(hvg_genes, gene_cluster_labels)}

nlgn3_candidates = [g for g in nlgn_found if "3" in g or "nlgn3" in g.lower()
                    or "NLGN3" in g]
focus_gene = nlgn3_candidates[0] if nlgn3_candidates else nlgn_found[0]
focus_module_id = gene_to_module[focus_gene]
module_genes = [g for g, m in gene_to_module.items() if m == focus_module_id]

print(f"\n  Focus gene  : {focus_gene}")
print(f"  Module ID   : {focus_module_id}")
print(f"  Module size : {len(module_genes)} genes")

# Save module membership
module_df = pd.DataFrame({
    "gene":   list(gene_to_module.keys()),
    "module": list(gene_to_module.values()),
})
module_csv = os.path.join(config.RESULTS_DIR, "05_gene_modules.csv")
module_df.to_csv(module_csv, index=False)
print(f"  Module assignments saved: {module_csv}")

focus_module_df = module_df[module_df["module"] == focus_module_id]
focus_csv = os.path.join(config.RESULTS_DIR, f"05_{focus_gene}_module_genes.csv")
focus_module_df.to_csv(focus_csv, index=False)
print(f"  {focus_gene} module genes saved: {focus_csv}")

# ── 6. FIGURE 7 — Co-expression network graph ─────────────────────
# Build a networkx graph for the focus module.
# Edges are drawn only where |r| > CORR_THRESHOLD to keep the graph readable.

# Guard: a module with 0 or 1 gene cannot form a meaningful network.
if len(module_genes) < 2:
    print(f"\n  Warning: {focus_gene} module contains only {len(module_genes)} gene(s) "
          f"— skipping network graph.\n"
          "  Try lowering CORR_THRESHOLD in config.py to merge more genes into modules.")
else:
    print(f"\n  Building network graph for {focus_gene} module ...")

    # Index positions within corr_matrix for module genes
    mod_gene_list = sorted(module_genes, key=lambda g: g != focus_gene)   # focus first
    mod_indices   = np.array([np.where(hvg_genes == g)[0][0] for g in mod_gene_list])
    mod_corr      = corr_matrix[np.ix_(mod_indices, mod_indices)]

    G = nx.Graph()
    G.add_nodes_from(mod_gene_list)

    for i in range(len(mod_gene_list)):
        for j in range(i + 1, len(mod_gene_list)):
            r = mod_corr[i, j]
            if abs(r) >= config.CORR_THRESHOLD:
                G.add_edge(mod_gene_list[i], mod_gene_list[j], weight=abs(r))

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Spring layout — focus_gene placed at centre
    pos = nx.spring_layout(G, seed=42, weight="weight")

    # Node colours: focus_gene = red, Nlgn family = orange, others = steelblue
    node_colors = []
    for node in G.nodes():
        if node == focus_gene:
            node_colors.append("#e74c3c")
        elif node in nlgn_found:
            node_colors.append("#e67e22")
        else:
            node_colors.append("#2980b9")

    node_sizes = [800 if n == focus_gene else
                  400 if n in nlgn_found else 150
                  for n in G.nodes()]

    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    # edge_alphas removed — was computed but never passed to any draw call

    fig7, ax7 = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_nodes(G, pos, ax=ax7, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9)
    nx.draw_networkx_edges(G, pos, ax=ax7, alpha=0.3,
                           width=[w * 2 for w in edge_weights],
                           edge_color="grey")

    # Label only high-connectivity nodes and Nlgn genes to avoid clutter
    degree_dict = dict(G.degree())
    top_degree  = sorted(degree_dict, key=degree_dict.get, reverse=True)[:15]
    label_nodes = set(top_degree) | set(nlgn_found)
    labels = {n: n for n in G.nodes() if n in label_nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax7, font_size=7)

    ax7.set_title(f"Figure 7 — Co-expression network: {focus_gene} module\n"
                  f"Edges: |Pearson r| ≥ {config.CORR_THRESHOLD}  "
                  f"({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)",
                  fontsize=11)
    ax7.axis("off")

    # Legend (Patch imported at top of file)
    legend_elements = [
        Patch(facecolor="#e74c3c", label=f"{focus_gene} (focus)"),
        Patch(facecolor="#e67e22", label="Other Nlgn isoforms"),
        Patch(facecolor="#2980b9", label="Co-expressed genes"),
    ]
    ax7.legend(handles=legend_elements, loc="lower left", fontsize=9)

    plt.tight_layout()
    fig7_path = os.path.join(config.FIGURES_DIR, "05b_coexpr_network.png")
    utils.save_fig(fig7, fig7_path)
    print(f"  → Figure 7 saved: {fig7_path}")

# ── 7. FIGURE 8 — GO Enrichment ───────────────────────────────────
print(f"\n  Running GO enrichment for {focus_gene} module ({len(module_genes)} genes) ...")

if len(module_genes) < 5:
    print(f"  Warning: module has only {len(module_genes)} genes — skipping GO enrichment.")
    print("  Try lowering CORR_THRESHOLD in config.py to get larger modules.")
else:
    try:
        from gprofiler import GProfiler
        organism = config.GPROFILER_ORGANISM[config.SPECIES]
        gp = GProfiler(return_dataframe=True)
        go_results = gp.profile(
            organism=organism,
            query=module_genes,
            sources=["GO:BP", "GO:MF", "KEGG"],
            significance_threshold_method="fdr",
            user_threshold=0.05,
            no_evidences=True,
        )

        if go_results.empty:
            print("  No significant GO terms found (FDR < 0.05).")
            print("  The module may be too small or the genes too diverse.")
        else:
            go_results = go_results.sort_values("p_value")
            go_csv = os.path.join(config.RESULTS_DIR,
                                  f"05_{focus_gene}_GO_enrichment.csv")
            go_results.to_csv(go_csv, index=False)
            print(f"  GO results saved: {go_csv}")
            print(f"  Significant terms: {len(go_results)}")

            # Plot top 20 terms
            top_go = go_results.head(20).copy()
            top_go["-log10(FDR)"] = -np.log10(top_go["p_value"].clip(1e-300))
            top_go = top_go.sort_values("-log10(FDR)")

            # Colour bars by source (GO:BP, GO:MF, KEGG)
            source_palette = {"GO:BP": "#3498db", "GO:MF": "#2ecc71", "KEGG": "#e74c3c"}
            bar_colors = [source_palette.get(s, "#95a5a6")
                          for s in top_go["source"]]

            fig8, ax8 = plt.subplots(figsize=(10, max(6, len(top_go) * 0.4)))
            # bars= assignment removed — return value was stored but never used
            ax8.barh(top_go["name"], top_go["-log10(FDR)"],
                     color=bar_colors, edgecolor="white", linewidth=0.5)
            ax8.axvline(x=-np.log10(0.05), color="red", linestyle="--")
            ax8.set_xlabel("-log₁₀(FDR-adjusted p-value)")
            ax8.set_title(f"Figure 8 — GO enrichment: {focus_gene} co-expression module\n"
                          f"({len(module_genes)} genes)", fontsize=11)

            # Single combined legend: source colours + FDR threshold line.
            # (Patch imported at top of file — no duplicate import needed here.)
            source_legend = [Patch(facecolor=c, label=s)
                             for s, c in source_palette.items()
                             if s in top_go["source"].values]
            ax8.legend(handles=source_legend + [
                mlines.Line2D([], [], color="red", linestyle="--", label="FDR = 0.05")
            ], fontsize=8, loc="lower right")

            plt.tight_layout()
            fig8_path = os.path.join(config.FIGURES_DIR, "05c_GO_enrichment.png")
            utils.save_fig(fig8, fig8_path)
            print(f"  → Figure 8 saved: {fig8_path}")

    except ImportError:
        print("  gprofiler-official not installed.")
        print("  Run:  pip install gprofiler-official  then re-run this script.")
    except Exception as e:
        print(f"  GO enrichment failed: {e}")
        print("  Check your internet connection and try again.")

# ── 8. Summary correlation table for Nlgn genes ──────────────────
# For each Nlgn gene, show the top 20 most correlated genes.
print("\n  Computing top co-expressed genes for each Nlgn isoform ...")
for gene in nlgn_found:
    if gene not in gene_to_module:
        continue
    gene_pos = np.where(hvg_genes == gene)[0][0]
    corr_with_gene = corr_matrix[gene_pos, :]

    top_idx = np.argsort(np.abs(corr_with_gene))[::-1][1:21]   # skip self (rank 0)
    top_genes  = hvg_genes[top_idx]
    top_corrs  = corr_with_gene[top_idx]

    top_df = pd.DataFrame({"gene": top_genes, "pearson_r": top_corrs.round(4)})
    top_csv = os.path.join(config.RESULTS_DIR, f"05_top_coexpr_{gene}.csv")
    top_df.to_csv(top_csv, index=False)
    print(f"  Top co-expressed genes for {gene} saved: {top_csv}")

utils.print_section("Step 5 complete — all figures generated")
print("  figures/05a_coexpr_dendrogram.png   (Figure 6)")
print("  figures/05b_coexpr_network.png      (Figure 7)")
print("  figures/05c_GO_enrichment.png       (Figure 8)")
print()
print("  Pipeline complete. Review all figures and CSV files in:")
print("  results/  and  figures/")
