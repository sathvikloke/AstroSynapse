# config.py
# ══════════════════════════════════════════════════════════════════
#  Neuroligin Isoform Expression Project — Central Configuration
#  Edit the USER SETTINGS section before running any script.
# ══════════════════════════════════════════════════════════════════

# ── USER SETTINGS ─────────────────────────────────────────────────

# Full path to your downloaded .h5ad file.
# See 00_download_instructions.txt for exact download steps.
DATA_PATH = "data/astrocytes.h5ad"

# "mouse" or "human"
SPECIES = "mouse"

# True  → your file already contains only astrocytes (e.g. Bayraktar 2020).
#          Script 02_cluster_annotate.py is NOT needed; start at 03.
# False → your file contains all brain cell types (e.g. Allen Brain Cell Atlas).
#          Run 02_cluster_annotate.py to identify and extract astrocytes first.
ASTROCYTES_ONLY = True

# ── DIRECTORIES (created automatically) ───────────────────────────
DATA_DIR    = "data"
FIGURES_DIR = "figures"
RESULTS_DIR = "results"

# ── QC THRESHOLDS ─────────────────────────────────────────────────
MIN_GENES_PER_CELL = 200    # cells expressing < 200 genes → likely empty droplets
MAX_GENES_PER_CELL = 6000   # cells expressing > 6000 genes → likely doublets
MIN_CELLS_PER_GENE = 3      # genes detected in < 3 cells → removed
MAX_MITO_PERCENT   = 10.0   # cells with > 10 % mitochondrial reads → likely damaged

# ── PREPROCESSING ─────────────────────────────────────────────────
N_TOP_HVG     = 2000   # highly variable genes used for PCA
N_PCS_COMPUTE = 50     # total PCs computed
N_PCS_USE     = 30     # PCs used to build the neighbour graph
N_NEIGHBORS   = 15     # k for the k-nearest-neighbour graph

# ── CLUSTERING ────────────────────────────────────────────────────
LEIDEN_RES_FULL       = 0.5   # resolution for whole-brain clustering
LEIDEN_RES_ASTROCYTES = 0.4   # resolution for astrocyte sub-clustering
                               # Tip: lower value → fewer, broader clusters

# ── GENES OF INTEREST ─────────────────────────────────────────────
NLGN_GENES = {
    "mouse": ["Nlgn1", "Nlgn2", "Nlgn3"],
    "human": ["NLGN1", "NLGN2", "NLGN3"],
}

NRXN_GENES = {
    "mouse": ["Nrxn1", "Nrxn2", "Nrxn3"],
    "human": ["NRXN1", "NRXN2", "NRXN3"],
}

# ── CELL TYPE MARKER GENES ────────────────────────────────────────
CELL_TYPE_MARKERS = {
    "mouse": {
        "Astrocyte":       ["Gfap", "Aqp4", "S100b", "Aldh1l1", "Slc1a2", "Slc1a3"],
        "Neuron":          ["Snap25", "Syt1", "Rbfox3", "Map2"],
        "Oligodendrocyte": ["Mbp", "Mog", "Plp1"],
        "OPC":             ["Pdgfra", "Cspg4", "Vcan"],
        "Microglia":       ["Aif1", "Cx3cr1", "Tmem119", "P2ry12"],
        "Endothelial":     ["Cldn5", "Pecam1", "Ly6c1"],
        "Pericyte":        ["Pdgfrb", "Acta2", "Rgs5"],
    },
    "human": {
        "Astrocyte":       ["GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "SLC1A3"],
        "Neuron":          ["SNAP25", "SYT1", "RBFOX3", "MAP2"],
        "Oligodendrocyte": ["MBP", "MOG", "PLP1"],
        "OPC":             ["PDGFRA", "CSPG4", "VCAN"],
        "Microglia":       ["AIF1", "CX3CR1", "TMEM119", "P2RY12"],
        "Endothelial":     ["CLDN5", "PECAM1"],
        "Pericyte":        ["PDGFRB", "ACTA2", "RGS5"],
    },
}

# Astrocyte subtype markers (Bayraktar et al. 2020, Zeisel et al. 2018)
ASTROCYTE_SUBTYPE_MARKERS = {
    "mouse": {
        "Grey_matter":      ["Aldoc", "Slc1a2", "Gja1", "Kcnj10"],
        "White_matter":     ["Gfap", "Id3", "Clu", "Gjb6"],
        "Bergmann_glia":    ["Hopx", "Ptprz1", "Ttyh1", "Fabp7"],
        "Reactive":         ["Lcn2", "Steap4", "C3", "Serpina3n"],
        "Perisynaptic":     ["Ezr", "Slc7a10", "Mfge8", "Apoe"],
    },
    "human": {
        "Grey_matter":      ["ALDOC", "SLC1A2", "GJA1", "KCNJ10"],
        "White_matter":     ["GFAP", "ID3", "CLU", "GJB6"],
        "Bergmann_glia":    ["HOPX", "PTPRZ1", "TTYH1", "FABP7"],
        "Reactive":         ["LCN2", "STEAP4", "C3", "SERPINA3"],
        "Perisynaptic":     ["EZR", "SLC7A10", "MFGE8", "APOE"],
    },
}

# ── CO-EXPRESSION ─────────────────────────────────────────────────
N_TOP_COEXPR_GENES = 500    # top HVGs used for the correlation matrix
MIN_EXPR_FRACTION  = 0.10   # gene must be expressed in >= 10 % of astrocytes
CORR_THRESHOLD     = 0.30   # |Pearson r| cutoff for drawing network edges

# ── G:PROFILER ORGANISM CODES ─────────────────────────────────────
GPROFILER_ORGANISM = {
    "mouse": "mmusculus",
    "human": "hsapiens",
}
