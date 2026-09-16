"""
Analyse SQANTI3 junctions.txt files across 4 platforms (3 replicates each).

For each platform, junctions detected in ALL 3 replicates are kept (matched
on chrom + genomic_start_coord + genomic_end_coord + strand) and counted once.

Two grouped barplots are produced side by side:
  1. % canonical vs non-canonical junctions
  2. % known vs novel junctions (junction_category column)
"""

__version__ = "1.0.0"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ─────────────────────────────────────────────
# USER CONFIGURATION — edit these values
# ─────────────────────────────────────────────

# For each platform, list the full paths to the 3 replicate junctions files
# Replace with your own paths (relative or absolute).
PLATFORMS = {
    "10X-3PRIME": [
        "data/sqanti3/10X-3PRIME_REP1_isoforms_junctions.txt",
        "data/sqanti3/10X-3PRIME_REP2_isoforms_junctions.txt",
        "data/sqanti3/10X-3PRIME_REP3_isoforms_junctions.txt",
    ],
    "10X-5PRIME": [
        "data/sqanti3/10X_5PRIME_REP1_isoforms_junctions.txt",
        "data/sqanti3/10X_5PRIME_REP2_isoforms_junctions.txt",
        "data/sqanti3/10X_5PRIME_REP3_isoforms_junctions.txt",
    ],
    "ARGENTAG": [
        "data/sqanti3/ARGENTAG_REP1_isoforms_junctions.txt",
        "data/sqanti3/ARGENTAG_REP2_isoforms_junctions.txt",
        "data/sqanti3/ARGENTAG_REP3_isoforms_junctions.txt",
    ],
    "PARSE": [
        "data/sqanti3/PARSE_REP1_isoforms_junctions.txt",
        "data/sqanti3/PARSE_REP2_isoforms_junctions.txt",
        "data/sqanti3/PARSE_REP3_isoforms_junctions.txt",
    ],
}

PLATFORM_COLORS = {
    "10X-3PRIME":"#658b66",
    "10X-5PRIME": "#8aa783",
    "ARGENTAG": "#2680a6",
    "PARSE": "#875670",
}

# Output figure path
OUTPUT_FIG = "junction_barplots.pdf"   # also accepts .png, .svg

# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────
 
def load_junctions_per_replicate(filepath: str) -> pd.DataFrame:
    """
    Read one SQANTI3 junctions.txt file and return a deduplicated DataFrame
    indexed by a composite junction_id (chrom|strand|start|end).
    Keeps one row per unique junction (same junction can appear in multiple isoforms).
    """
    df = pd.read_csv(filepath, sep="\t", low_memory=False)
    df.columns = df.columns.str.strip()
 
    for col in ["canonical", "junction_category"]:
        df[col] = df[col].str.strip().str.lower()
 
    df["junction_id"] = (
        df["chrom"].astype(str) + "|" +
        df["strand"].astype(str) + "|" +
        df["genomic_start_coord"].astype(str) + "|" +
        df["genomic_end_coord"].astype(str)
    )
 
    df = df.drop_duplicates(subset="junction_id")
    return df[["junction_id", "canonical", "junction_category"]].set_index("junction_id")
 
 
def load_junction_data(filepaths: list) -> pd.DataFrame:
    """
    Load junctions from all replicates, intersect on junction_id,
    and return a single deduplicated DataFrame counted once per junction.
    canonical and junction_category are deterministic for a given splice site,
    so values are taken from the first replicate.
    """
    rep_dfs = []
    for f in filepaths:
        rep_df = load_junctions_per_replicate(f)
        print(f"  {f}: {len(rep_df):,} unique junctions")
        rep_dfs.append(rep_df)
 
    common_ids = set(rep_dfs[0].index)
    for rep_df in rep_dfs[1:]:
        common_ids &= set(rep_df.index)
    print(f"  --> {len(common_ids):,} junctions detected in all {len(filepaths)} replicates\n")
 
    consensus = rep_dfs[0].loc[list(common_ids)].reset_index(drop=True)
    return consensus
 
 
platform_data = {}
for name, paths in PLATFORMS.items():
    print(f"Loading {name} ...")
    platform_data[name] = load_junction_data(paths)
 
# ─────────────────────────────────────────────
# COMPUTE PERCENTAGES
# ─────────────────────────────────────────────
 
canonical_pct      = {}  # % canonical / non-canonical out of all junctions
canon_knownovel    = {}  # % known / novel within canonical junctions
noncanon_knownovel = {}  # % known / novel within non-canonical junctions
 
for name, df in platform_data.items():
    n = len(df)
 
    # canonical / non-canonical
    cc = df["canonical"].value_counts()
    canonical_pct[name] = {
        "canonical":     cc.get("canonical",     0) / n * 100,
        "non_canonical": cc.get("non_canonical", 0) / n * 100,
    }
 
    # known / novel within canonical
    df_canon  = df[df["canonical"] == "canonical"]
    n_canon   = len(df_canon)
    if n_canon > 0:
        c = df_canon["junction_category"].value_counts()
        canon_knownovel[name] = {
            "known": c.get("known", 0) / n_canon * 100,
            "novel": c.get("novel", 0) / n_canon * 100,
        }
    else:
        canon_knownovel[name] = {"known": 0.0, "novel": 0.0}
 
    # known / novel within non-canonical
    df_nc   = df[df["canonical"] == "non_canonical"]
    n_nc    = len(df_nc)
    if n_nc > 0:
        nc = df_nc["junction_category"].value_counts()
        noncanon_knownovel[name] = {
            "known": nc.get("known", 0) / n_nc * 100,
            "novel": nc.get("novel", 0) / n_nc * 100,
        }
    else:
        noncanon_knownovel[name] = {"known": 0.0, "novel": 0.0}
 
platforms = list(PLATFORMS.keys())
colors    = [PLATFORM_COLORS[p] for p in platforms]
x         = np.arange(len(platforms))
bar_width = 0.35
 
# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────
 
def label_bars(ax, bars):
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.5,
                f"{h:.1f}%",
                ha="center", va="bottom", fontsize=8,
            )
 
def annotate_n(ax, n_dict):
    for i, p in enumerate(platforms):
        ax.text(x[i], -8, f"n={n_dict[p]:,}", ha="center", va="top",
                fontsize=8, color="gray", transform=ax.get_xaxis_transform())
 
def style_ax(ax, title):
    ax.set_xticks(x)
    ax.set_xticklabels(platforms, fontsize=10)
    ax.set_ylabel("% of junctions", fontsize=12)
    ax.set_ylim(0, 110)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_title(title, fontsize=12, pad=10)
    ax.legend(fontsize=9, framealpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
 
# ─────────────────────────────────────────────
# PLOT — three panels
# ─────────────────────────────────────────────
 
fig, axes = plt.subplots(1, 3, figsize=(20, 5))
 
# ── Panel 1: canonical vs non-canonical ──────
ax = axes[0]
bars1 = ax.bar(x - bar_width / 2,
               [canonical_pct[p]["canonical"]     for p in platforms],
               bar_width, color=colors, alpha=0.90, label="Canonical", edgecolor="white")
bars2 = ax.bar(x + bar_width / 2,
               [canonical_pct[p]["non_canonical"] for p in platforms],
               bar_width, color=colors, alpha=0.45, label="Non-canonical",
               edgecolor="white", hatch="///")
label_bars(ax, bars1)
label_bars(ax, bars2)
style_ax(ax, "Canonical vs Non-canonical junctions\n(detected in all replicates)")
annotate_n(ax, {p: len(platform_data[p]) for p in platforms})
 
# ── Panel 2: known/novel within canonical ────
ax = axes[1]
# n shown below = number of canonical junctions per platform
n_canon_dict = {
    p: int(round(len(platform_data[p]) * canonical_pct[p]["canonical"] / 100))
    for p in platforms
}
bars3 = ax.bar(x - bar_width / 2,
               [canon_knownovel[p]["known"] for p in platforms],
               bar_width, color=colors, alpha=0.90, label="Known", edgecolor="white")
bars4 = ax.bar(x + bar_width / 2,
               [canon_knownovel[p]["novel"] for p in platforms],
               bar_width, color=colors, alpha=0.45, label="Novel",
               edgecolor="white", hatch="///")
label_bars(ax, bars3)
label_bars(ax, bars4)
style_ax(ax, "Known vs Novel — canonical junctions\n(detected in all replicates)")
annotate_n(ax, n_canon_dict)
 
# ── Panel 3: known/novel within non-canonical
ax = axes[2]
n_noncanon_dict = {
    p: int(round(len(platform_data[p]) * canonical_pct[p]["non_canonical"] / 100))
    for p in platforms
}
bars5 = ax.bar(x - bar_width / 2,
               [noncanon_knownovel[p]["known"] for p in platforms],
               bar_width, color=colors, alpha=0.90, label="Known", edgecolor="white")
bars6 = ax.bar(x + bar_width / 2,
               [noncanon_knownovel[p]["novel"] for p in platforms],
               bar_width, color=colors, alpha=0.45, label="Novel",
               edgecolor="white", hatch="///")
label_bars(ax, bars5)
label_bars(ax, bars6)
style_ax(ax, "Known vs Novel — non-canonical junctions\n(detected in all replicates)")
annotate_n(ax, n_noncanon_dict)
 
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300)
print(f"Figure saved -> {OUTPUT_FIG}")
plt.show()