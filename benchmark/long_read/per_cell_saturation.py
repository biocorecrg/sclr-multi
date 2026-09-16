#!/usr/bin/env python3
"""
per_cell_saturation.py — sequencing saturation from count matrices, per cell

Uses the same YAML config as sequencing_saturation.py.  Each category must
additionally define a 'matrix' field pointing to its count matrix.

Algorithm per cell:
  1. Expand UMI counts into a pseudo-read pool:
       feature j with count c[j]  →  c[j] copies of index j in the pool
  2. Shuffle the pool once (random read order) and subsample at increasing depths.
  3. Count unique features seen at each depth.
  4. Saturation = 1 - unique_features / UMIs_sampled.

Category curves are the mean ± 1 SD across N randomly selected cells.
"""

__version__ = "1.0.0"
All categories are subsampled to the same UMI depth so the curves are
directly comparable regardless of the original sequencing depth.

Note: this is a rarefaction-based saturation estimate computed from the count
matrix (deduplicated UMIs).  It answers "how many new features does each
additional UMI reveal?"  It differs from PCR-duplicate-based saturation
(1 - unique_reads / total_reads), which requires the raw BAM.

Supported matrix formats:
  *.h5ad          AnnData (scanpy)
  *.h5            10x Genomics HDF5 (scanpy)
  *.csv / *.tsv   delimited text; first row = header, first column = index;
                  orientation auto-detected: if n_rows > n_cols the matrix
                  is assumed to be features × cells and transposed
  directory       10x MTX  (matrix.mtx[.gz] + barcodes.tsv[.gz]
                  + features/genes.tsv[.gz])

Example config addition (add to existing sequencing_saturation.py config):

  categories:
    10X-3PRIME:
      color: "#658b66"
      matrix: /path/to/10x3prime_matrix.h5ad
      replicates: [...]   # kept for sequencing_saturation.py; ignored here
    PARSE:
      color: "#875670"
      matrix: /path/to/parse_matrix.h5ad
      replicates: [...]

  # New keys for this script:
  n_cells:           100    # cells to draw per category (same for all)
  min_umis_per_cell: 200    # exclude cells below this UMI threshold
  subsample_to:      null   # null = auto (min total UMIs across selected cells)
  output:            saturation_cells.png
"""

import argparse
import colorsys
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import scipy.sparse as sp
import yaml


# ---------------------------------------------------------------------------
# Matrix I/O
# ---------------------------------------------------------------------------

def load_matrix(path: str):
    """
    Load a count matrix and return a (cells × features) CSR matrix.

    Keeps data sparse throughout to minimise RAM; per-cell dense extraction
    happens later in get_cell_counts().
    """
    p = Path(path)
    ext = p.suffix.lower()

    if ext == ".h5ad":
        import scanpy as sc
        adata = sc.read_h5ad(path)
        X = adata.X
    elif ext in (".h5", ".hdf5"):
        import scanpy as sc
        adata = sc.read_10x_h5(path)
        X = adata.X
    elif ext in (".csv", ".tsv"):
        import pandas as pd
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(path, sep=sep, index_col=0)
        # Heuristic: more rows than columns usually means features × cells
        if df.shape[0] > df.shape[1]:
            df = df.T
        X = df.values
    elif p.is_dir():
        import scanpy as sc
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        X = adata.X
    else:
        raise ValueError(
            f"Unsupported format: {path!r}. "
            "Expected .h5ad, .h5, .csv, .tsv, or a 10x MTX directory."
        )

    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    return X.astype(np.float32)


def cell_totals(X) -> np.ndarray:
    """Total UMIs per cell (row sums)."""
    return np.asarray(X.sum(axis=1)).ravel()


def get_cell_counts(X, i: int) -> np.ndarray:
    """Extract cell i as a flat int32 dense array."""
    row = X[i]
    if sp.issparse(row):
        row = row.toarray()
    return np.asarray(row, dtype=np.int32).ravel()


# ---------------------------------------------------------------------------
# Per-cell saturation curve
# ---------------------------------------------------------------------------

def per_cell_curve(count_vector: np.ndarray, n_points: int,
                   subsample_max: int, seed: int):
    """
    Rarefaction-based saturation curve for a single cell.

    Parameters
    ----------
    count_vector  : 1D int32 array of UMI counts, one entry per feature
    subsample_max : UMI depth to subsample to (capped at cell's total)
    seed          : per-cell RNG seed for reproducibility

    Returns
    -------
    (depths, unique_counts, saturation)  or  None if the cell is empty.

    Uses the same O(n log n) first-occurrence + searchsorted algorithm as
    sequencing_saturation.py: shuffle once, find where each unique feature
    first appears, then answer all depth queries with binary search.
    """
    nz_mask = count_vector > 0
    nonzero_counts = count_vector[nz_mask]
    total = int(nonzero_counts.sum())
    if total == 0:
        return None

    # Expand: feature index j repeated count[j] times
    nz_idx = np.where(nz_mask)[0].astype(np.int32)
    reads = np.repeat(nz_idx, nonzero_counts)

    cap = min(subsample_max, total)
    rng = np.random.default_rng(seed)
    rng.shuffle(reads)
    reads = reads[:cap]

    # For each unique feature, the index in `reads` where it first appears
    _, first_occ = np.unique(reads, return_index=True)
    first_occ_sorted = np.sort(first_occ)

    depths = np.unique(
        np.round(np.linspace(0, cap, n_points + 1)).astype(int)
    )
    depths = depths[depths > 0]

    # searchsorted(..., side='left') counts first occurrences strictly < d,
    # which equals the number of unique features in reads[:d]
    unique_counts = np.searchsorted(first_occ_sorted, depths, side="left")
    saturation = 1.0 - float(unique_counts[-1]) / float(depths[-1])

    return depths.astype(float), unique_counts.astype(float), saturation


# ---------------------------------------------------------------------------
# Colour helper  (shared with sequencing_saturation.py)
# ---------------------------------------------------------------------------

def _lighten_color(color, amount=0.5):
    """Lighten a matplotlib color by blending with white."""
    try:
        c = mc.cnames[color]
    except (KeyError, TypeError):
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_saturation(results: dict, output: str, n_cells: int):
    """Median saturation curves per category (no variability bands)."""
    fig, ax = plt.subplots(figsize=(9, 6))

    for cat, d in results.items():
        xs = d["depths"] / 1e3          # thousands of UMIs per cell
        sat_pct = d["saturation"] * 100
        label = f"{cat}  (saturation: {sat_pct:.1f}%)"

        ax.plot(xs, d["median_uniq"], color=d["color"], linewidth=2.5,
                marker="o", markersize=4, label=label, zorder=3)

    ax.set_xlabel("UMIs sampled per cell (thousands)", fontsize=13)
    ax.set_ylabel("Unique features detected per cell", fontsize=12)
    ax.set_title(
        f"Sequencing Saturation — Per Cell  (n = {n_cells} cells / category)",
        fontsize=15, fontweight="bold"
    )
    ax.legend(frameon=True, fontsize=11, loc="lower right")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

    plt.tight_layout()
    plt.savefig(output, dpi=180, bbox_inches="tight")
    print(f"\nPlot saved → {output}")
    plt.close(fig)


def plot_saturation_iqr_subplots(results: dict, output: str, n_cells: int):
    """
    One subplot per category showing median + IQR band (25th–75th percentile),
    arranged in a 2-column grid. All subplots share axis limits for comparison.
    """
    cats = list(results.keys())
    n_cat = len(cats)
    n_cols = 2
    n_rows = (n_cat + 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(7 * n_cols, 5 * n_rows),
        sharex=True, sharey=True,
    )
    axes_flat = np.array(axes).ravel()

    for ax, cat in zip(axes_flat, cats):
        d = results[cat]
        xs = d["depths"] / 1e3
        sat_pct = d["saturation"] * 100

        ax.fill_between(
            xs, d["q25_uniq"], d["q75_uniq"],
            color=d["color"], alpha=0.20, zorder=2, label="IQR (25–75%)"
        )
        ax.plot(xs, d["median_uniq"], color=d["color"], linewidth=2.5,
                marker="o", markersize=4, zorder=3,
                label=f"median  (saturation: {sat_pct:.1f}%)")

        ax.set_title(cat, fontsize=13, fontweight="bold", color=d["color"])
        ax.legend(frameon=True, fontsize=9, loc="lower right")
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

    for ax in axes_flat[n_cat:]:
        ax.set_visible(False)

    fig.supxlabel("UMIs sampled per cell (thousands)", fontsize=13, y=0.02)
    fig.supylabel("Unique features detected per cell", fontsize=12, x=0.02)
    fig.suptitle(
        f"Sequencing Saturation — Per Category  (n = {n_cells} cells / category)",
        fontsize=15, fontweight="bold", y=1.01
    )

    plt.tight_layout()
    plt.savefig(output, dpi=180, bbox_inches="tight")
    print(f"IQR subplots saved → {output}")
    plt.close(fig)


def plot_per_cell(results: dict, output: str):
    """
    Individual per-cell curves (light, thin, dashed) overlaid with the
    category mean (bold, solid).
    """
    fig, ax = plt.subplots(figsize=(11, 7))

    for cat, d in results.items():
        light = _lighten_color(d["color"], amount=0.55)
        xs = d["depths"] / 1e3

        for cell_uniq in d["cell_curves"]:
            ax.plot(xs, cell_uniq, color=light, linewidth=0.8,
                    linestyle="--", alpha=0.4, zorder=2)

        sat_pct = d["saturation"] * 100
        label = f"{cat} [median]  (saturation: {sat_pct:.1f}%)"
        ax.plot(xs, d["median_uniq"], color=d["color"], linewidth=2.8,
                marker="o", markersize=5, label=label, zorder=3)

    ax.set_xlabel("UMIs sampled per cell (thousands)", fontsize=13)
    ax.set_ylabel("Unique features detected per cell", fontsize=12)
    ax.set_title(
        "Sequencing Saturation — Individual Cells",
        fontsize=15, fontweight="bold"
    )
    n_cat = len(results)
    ax.legend(frameon=True, fontsize=10, loc="lower right",
              ncol=1 + (n_cat > 4))
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

    plt.tight_layout()
    plt.savefig(output, dpi=180, bbox_inches="tight")
    print(f"Per-cell plot saved → {output}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# TSV report
# ---------------------------------------------------------------------------

def write_table(results: dict, output: str):
    """Write median + IQR saturation table."""
    with open(output, "w") as fh:
        fh.write(
            "# Per-cell sequencing saturation (rarefaction on count matrix)\n"
            "# saturation = 1 - median_unique_features / umis_sampled\n"
            "#\n"
        )
        fh.write(
            "category\tn_cells\tumis_sampled\t"
            "median_unique_features\tq25_unique_features\tq75_unique_features\t"
            "duplication_rate\tsaturation_pct\n"
        )
        for cat, d in results.items():
            n = len(d["cell_curves"])
            for depth, med, q25, q75 in zip(
                d["depths"], d["median_uniq"], d["q25_uniq"], d["q75_uniq"]
            ):
                dup = 1.0 - med / depth if depth > 0 else 0.0
                fh.write(
                    f"{cat}\t{n}\t{int(depth)}\t"
                    f"{med:.2f}\t{q25:.2f}\t{q75:.2f}\t"
                    f"{dup:.6f}\t{dup * 100:.2f}\n"
                )
    print(f"Table saved → {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Per-cell sequencing saturation from count matrices. "
            "Uses the same YAML config as sequencing_saturation.py; "
            "each category must have a 'matrix' field."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True,
                        help="Path to YAML config file")
    parser.add_argument("--n-cells", type=int, default=None,
                        help="Cells to draw per category (overrides config n_cells)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed (overrides config seed)")
    args = parser.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    n_cells      = args.n_cells if args.n_cells is not None else cfg.get("n_cells", 100)
    n_points     = cfg.get("n_points", 20)
    seed         = args.seed   if args.seed   is not None else cfg.get("seed", 42)
    min_umis     = cfg.get("min_umis_per_cell", 200)
    subsample_to = cfg.get("subsample_to", None)
    output       = cfg.get("output", "saturation_cells.png")
    output_table = cfg.get(
        "output_table", os.path.splitext(output)[0] + "_table.tsv"
    )
    categories   = cfg["categories"]

    stem, ext = os.path.splitext(output)
    output_per_cell   = stem + "_per_cell" + ext
    output_sd_subplots = stem + "_sd_subplots" + ext

    print(f"\n=== Per-cell sequencing saturation (matrix-based) ===")
    print(f"    n_cells={n_cells}  n_points={n_points}  seed={seed}")
    print(f"    min_umis_per_cell={min_umis}\n")

    # ------------------------------------------------------------------
    # Load matrices and subsample cells
    # ------------------------------------------------------------------
    rng = np.random.default_rng(seed)
    category_cells = {}

    for cat, info in categories.items():
        if "matrix" not in info:
            sys.exit(
                f"ERROR: category '{cat}' has no 'matrix' field in config.\n"
                "Add a 'matrix:' entry pointing to its count matrix."
            )

        t0 = time.time()
        print(f"  Loading {cat}: {info['matrix']}")
        X = load_matrix(info["matrix"])
        totals = cell_totals(X)

        valid = np.where(totals >= min_umis)[0]
        print(f"    {X.shape[0]:,} total cells  |  "
              f"{len(valid):,} with >= {min_umis} UMIs")

        if len(valid) < n_cells:
            print(f"    WARNING: only {len(valid)} valid cells — "
                  f"using all (requested {n_cells})")
            chosen = valid
        else:
            chosen = rng.choice(valid, size=n_cells, replace=False)

        chosen_totals = totals[chosen]
        elapsed = time.time() - t0
        print(f"    Selected {len(chosen)} cells  |  "
              f"median {int(np.median(chosen_totals)):,} UMIs/cell  |  "
              f"min {int(chosen_totals.min()):,}  |  "
              f"max {int(chosen_totals.max()):,}  |  "
              f"{elapsed:.1f}s")

        category_cells[cat] = {
            "X":      X[chosen],
            "totals": chosen_totals,
            "color":  info["color"],
        }

    # ------------------------------------------------------------------
    # Determine common subsampling depth
    # ------------------------------------------------------------------
    if subsample_to is None:
        subsample_to = int(
            min(d["totals"].min() for d in category_cells.values())
        )
        print(f"\nAuto subsample depth: {subsample_to:,} UMIs/cell "
              f"(= min total UMIs across all selected cells)")
    else:
        print(f"\nSubsample depth (from config): {subsample_to:,} UMIs/cell")
        for cat, d in category_cells.items():
            under = int((d["totals"] < subsample_to).sum())
            if under:
                print(f"    WARNING: {under} cells in '{cat}' have fewer "
                      f"than {subsample_to:,} UMIs — capped at their total.")

    # ------------------------------------------------------------------
    # Compute per-cell saturation curves
    # ------------------------------------------------------------------
    results = {}

    for cat, d in category_cells.items():
        X_sub = d["X"]
        n     = X_sub.shape[0]
        print(f"\nComputing saturation curves: {cat}  ({n} cells) ...")
        t0 = time.time()

        cell_curves  = []
        common_depths = None

        for i in range(n):
            counts = get_cell_counts(X_sub, i)
            curve  = per_cell_curve(counts, n_points, subsample_to, seed + i)
            if curve is None:
                continue
            depths, unique_counts, _ = curve
            cell_curves.append(unique_counts)
            if common_depths is None:
                common_depths = depths

        if not cell_curves:
            print(f"  WARNING: no valid cells for '{cat}', skipping.")
            continue

        cell_matrix  = np.stack(cell_curves)     # (n_valid, n_points)
        median_uniq  = np.median(cell_matrix, axis=0)
        q25_uniq     = np.percentile(cell_matrix, 25, axis=0)
        q75_uniq     = np.percentile(cell_matrix, 75, axis=0)
        saturation   = 1.0 - float(median_uniq[-1]) / float(common_depths[-1])

        elapsed = time.time() - t0
        iqr_at_max = q75_uniq[-1] - q25_uniq[-1]
        print(f"  -> Saturation @ {int(common_depths[-1]):,} UMIs/cell: "
              f"{saturation * 100:.1f}%  |  "
              f"median unique features: {median_uniq[-1]:.0f}  "
              f"IQR [{q25_uniq[-1]:.0f}, {q75_uniq[-1]:.0f}]  |  "
              f"{elapsed:.1f}s")

        results[cat] = {
            "color":        d["color"],
            "depths":       common_depths,
            "median_uniq":  median_uniq,
            "q25_uniq":     q25_uniq,
            "q75_uniq":     q75_uniq,
            "saturation":   saturation,
            "cell_curves":  cell_curves,
        }

    if not results:
        sys.exit("ERROR: no valid results to plot.")

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------
    plot_saturation(results, output, n_cells)
    plot_saturation_iqr_subplots(results, output_sd_subplots, n_cells)
    plot_per_cell(results, output_per_cell)
    write_table(results, output_table)


if __name__ == "__main__":
    main()
