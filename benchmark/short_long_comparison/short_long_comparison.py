"""
Compare Short-read and Long-read Single-cell RNA-seq Datasets
==============================================================

This script performs comprehensive comparison of same-sample sequencing
using different modalities, including:
- Barcode filtering FIRST (extracting barcode from obs_names with format: barcode-additional_data)
- Pseudobulk aggregation using decoupler (on filtered data only)
- CPM normalization
- Correlation analysis (shared genes and all genes)
- Venn diagram visualization
"""

import argparse
import re
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from matplotlib_venn import venn2
import warnings
warnings.filterwarnings('ignore')


def load_and_prepare_data(short_read_path, long_read_path):
    """
    Load count matrices from both modalities.
    
    Parameters:
    -----------
    short_read_path : str
        Path to short-read count matrix (h5ad, h5, csv, etc.)
    long_read_path : str
        Path to long-read count matrix (h5ad, h5, csv, etc.)
        
    Returns:
    --------
    adata_short : AnnData object
    adata_long : AnnData object
    """
    print("=" * 80)
    print("LOADING DATA")
    print("=" * 80)
    
    def load_file(filepath):
        """Helper function to load various file formats."""
        if filepath.endswith('.h5ad'):
            return sc.read_h5ad(filepath)
        elif filepath.endswith('.h5'):
            return sc.read_10x_h5(filepath)
        elif filepath.endswith('.csv'):
            df = pd.read_csv(filepath, index_col=0)
            return sc.AnnData(X=df.T.values, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        else:
            raise ValueError(f"Unsupported format: {filepath}. Use .h5ad, .h5, or .csv")
    
    # Load short-read data
    adata_short = load_file(short_read_path)

    # Load long-read data
    adata_long = load_file(long_read_path)
    
    # Calculate per-cell count statistics
    short_cell_counts = np.asarray(adata_short.X.sum(axis=1)).flatten()
    long_cell_counts = np.asarray(adata_long.X.sum(axis=1)).flatten()

    print(f"\nShort-read data loaded:")
    print(f"  Shape: {adata_short.n_obs} cells × {adata_short.n_vars} genes")
    print(f"  Barcodes: {adata_short.n_obs}")
    print(f"  Sample obs_names: {list(adata_short.obs_names[:3])}")
    print(f"  Median counts per cell: {np.median(short_cell_counts):.0f}")
    print(f"  Mean counts per cell: {short_cell_counts.mean():.0f}")

    print(f"\nLong-read data loaded:")
    print(f"  Shape: {adata_long.n_obs} cells × {adata_long.n_vars} genes")
    print(f"  Barcodes: {adata_long.n_obs}")
    print(f"  Sample obs_names: {list(adata_long.obs_names[:3])}")
    print(f"  Median counts per cell: {np.median(long_cell_counts):.0f}")
    print(f"  Mean counts per cell: {long_cell_counts.mean():.0f}")
    
    return adata_short, adata_long


def extract_barcodes(obs_names, barcode_sep='-'):
    """
    Extract barcodes from obs_names.

    Handles different barcode formats:
    - Numeric: XXXX-XXXX-XXXX-suffix (extracts first 3 numeric parts)
    - Nucleotide: ACGTACGT...-suffix (extracts all before the suffix)

    Parameters:
    -----------
    obs_names : pd.Index or list
        Observation names (format: barcode-additional_data)
    barcode_sep : str
        Separator to split barcode from additional data (default: '-')

    Returns:
    --------
    barcodes : dict
        Dictionary mapping barcode -> list of original indices with this barcode
    """
    # Barcode patterns anchored at the start of obs_names:
    #   10X      (nucleotide): AAACCAAAGGGCTTAC-1-10X-3PRIME_REP1 / AAACCCGCAAACCGCA-10X-3PRIME_REP1
    #   Argentag (numeric '-'): 0140-0117-0307-1 / 0076-0000-0481-ARGENTAG_REP1
    #   Parse    (numeric '_'): 05_01_91__s1 / 01_02_13-PARSE_REP1
    patterns = {
        '10X': re.compile(r'^([ACGTN]{6,})(?![A-Za-z])'),
        'numeric': re.compile(r'^(\d+[-_]\d+[-_]\d+)(?!\d)'),
    }

    barcodes = {}
    n_matched = {key: 0 for key in patterns}
    n_fallback = 0
    for idx, name in enumerate(obs_names):
        name = str(name)
        for key, pattern in patterns.items():
            match = pattern.match(name)
            if match:
                barcode = match.group(1)
                if key == 'numeric':
                    # Canonicalise separators so '-' and '_' variants still match
                    barcode = re.sub(r'[-_]', '_', barcode)
                n_matched[key] += 1
                break
        else:
            # Unknown format: fall back to the first field before the separator
            barcode = name.split(barcode_sep)[0]
            n_fallback += 1

        barcodes.setdefault(barcode, []).append(idx)

    print(f"  Barcode formats detected: 10X={n_matched['10X']}, "
          f"numeric (Argentag/Parse)={n_matched['numeric']}, fallback={n_fallback}")
    n_dup = sum(1 for idx_list in barcodes.values() if len(idx_list) > 1)
    if n_dup:
        print(f"  WARNING: {n_dup} barcodes map to more than one obs_name "
              f"(e.g. same barcode in different sublibraries/replicates)")

    return barcodes


def filter_by_shared_barcodes(adata_short, adata_long, barcode_sep_short='-', barcode_sep_long='-'):
    """
    FIRST STEP: Filter both datasets to keep only shared barcodes.

    Parameters:
    -----------
    adata_short : AnnData
    adata_long : AnnData
    barcode_sep_short : str
        Separator used in short-read obs_names (default: '-')
    barcode_sep_long : str
        Separator used in long-read obs_names (default: '-')
        
    Returns:
    --------
    adata_short_filtered : AnnData
    adata_long_filtered : AnnData
    report : dict with filtering statistics
    """
    print("\n" + "=" * 80)
    print("STEP 1: BARCODE FILTERING (IDENTIFY SHARED BARCODES)")
    print("=" * 80)
    
    print(f"  Short-read barcode separator: {barcode_sep_short!r}")
    print(f"  Long-read  barcode separator: {barcode_sep_long!r}")

    # Extract barcodes from both datasets
    barcodes_short = extract_barcodes(adata_short.obs_names, barcode_sep_short)
    barcodes_long = extract_barcodes(adata_long.obs_names, barcode_sep_long)
    
    short_barcode_set = set(barcodes_short.keys())
    long_barcode_set = set(barcodes_long.keys())
    
    # Find common barcodes
    common_barcodes = sorted(short_barcode_set & long_barcode_set)
    
    print(f"\nBarcode statistics:")
    print(f"  Short-read unique barcodes: {len(short_barcode_set)}")
    print(f"  Long-read unique barcodes: {len(long_barcode_set)}")
    print(f"  Common barcodes: {len(common_barcodes)}")
    print(f"  Overlap (% of short-read): {100 * len(common_barcodes) / len(short_barcode_set):.1f}%")
    print(f"  Overlap (% of long-read): {100 * len(common_barcodes) / len(long_barcode_set):.1f}%")
    
    # Filter to keep only common barcodes
    short_indices = []
    for bc in common_barcodes:
        short_indices.extend(barcodes_short[bc])
    
    long_indices = []
    for bc in common_barcodes:
        long_indices.extend(barcodes_long[bc])
    
    # Create filtered datasets
    adata_short_filtered = adata_short[short_indices].copy()
    adata_long_filtered = adata_long[long_indices].copy()

    # Calculate per-cell counts for filtered data
    short_filtered_counts = np.asarray(adata_short_filtered.X.sum(axis=1)).flatten()
    long_filtered_counts = np.asarray(adata_long_filtered.X.sum(axis=1)).flatten()

    print(f"\nAfter filtering to common barcodes:")
    print(f"  Short-read: {len(barcodes_short)} cells → {len(adata_short_filtered)} cells")
    print(f"    Median counts per cell: {np.median(short_filtered_counts):.0f}")
    print(f"  Long-read: {len(barcodes_long)} cells → {len(adata_long_filtered)} cells")
    print(f"    Median counts per cell: {np.median(long_filtered_counts):.0f}")
    
    report = {
        'short_read_original': len(barcodes_short),
        'long_read_original': len(barcodes_long),
        'common_barcodes': len(common_barcodes),
        'overlap_pct_short': 100 * len(common_barcodes) / len(short_barcode_set),
        'overlap_pct_long': 100 * len(common_barcodes) / len(long_barcode_set),
    }
    
    return adata_short_filtered, adata_long_filtered, report


def filter_zero_count_genes(adata_short, adata_long):
    """
    Report gene counts per modality and filter out genes with no counts.

    Parameters:
    -----------
    adata_short : AnnData
    adata_long : AnnData

    Returns:
    --------
    adata_short_filtered : AnnData
    adata_long_filtered : AnnData
    report : dict with gene filtering statistics
    """
    print("\n" + "=" * 80)
    print("STEP 2: GENE FILTERING (REMOVE ZERO-COUNT GENES)")
    print("=" * 80)

    print(f"\nGenes per modality (before filtering):")
    print(f"  Short-read: {adata_short.n_vars} genes")
    print(f"  Long-read:  {adata_long.n_vars} genes")

    adata_short_filtered = adata_short.copy()
    adata_long_filtered  = adata_long.copy()

    sc.pp.filter_genes(adata_short_filtered, min_counts=1)
    sc.pp.filter_genes(adata_long_filtered,  min_counts=1)

    print(f"\nGenes per modality (after removing zero-count genes):")
    print(f"  Short-read: {adata_short.n_vars} → {adata_short_filtered.n_vars} genes "
          f"({adata_short.n_vars - adata_short_filtered.n_vars} removed)")
    print(f"  Long-read:  {adata_long.n_vars} → {adata_long_filtered.n_vars} genes "
          f"({adata_long.n_vars - adata_long_filtered.n_vars} removed)")

    report = {
        'short_genes_before':  adata_short.n_vars,
        'short_genes_after':   adata_short_filtered.n_vars,
        'short_genes_removed': adata_short.n_vars - adata_short_filtered.n_vars,
        'long_genes_before':   adata_long.n_vars,
        'long_genes_after':    adata_long_filtered.n_vars,
        'long_genes_removed':  adata_long.n_vars - adata_long_filtered.n_vars,
    }

    return adata_short_filtered, adata_long_filtered, report


def preprocess_gene_names(var_names, dataset_name=""):
    """
    Preprocess gene names for standardization.
    Currently converts underscores to hyphens.
    
    Parameters:
    -----------
    var_names : pd.Index or list
        Gene names
    dataset_name : str
        Name of dataset for logging
        
    Returns:
    --------
    processed_names : pd.Index
        Processed gene names
    """
    processed = [str(name).replace('_', '-') for name in var_names]
    
    # Check if any genes were modified
    n_modified = sum(1 for orig, proc in zip(var_names, processed) if orig != proc)
    if n_modified > 0:
        print(f"\n  Gene name preprocessing ({dataset_name}):")
        print(f"    Modified {n_modified}/{len(var_names)} gene names (replaced '_' with '-')")
        print(f"    Examples:")
        modified_examples = [(orig, proc) for orig, proc in zip(var_names, processed) if orig != proc][:5]
        for orig, proc in modified_examples:
            print(f"      {orig} → {proc}")
    
    return pd.Index(processed)


def compute_pseudobulk(adata, preprocess=False, dataset_name=""):
    """
    Compute pseudobulk by summing raw counts across all cells per gene.
    Produces one aggregated profile for the entire dataset (one vector per matrix).

    Parameters:
    -----------
    adata : AnnData
        Single-cell data (ALREADY FILTERED TO SHARED BARCODES)
    preprocess : bool
        If True, preprocess gene names (convert _ to -)
    dataset_name : str
        Name of dataset for logging

    Returns:
    --------
    pseudobulk_counts : pd.Series
        Gene counts (summed across all cells)
    """
    print(f"  Summing counts across all cells...")
    counts = np.asarray(adata.X.sum(axis=0)).flatten()
    pseudobulk_series = pd.Series(counts, index=adata.var_names)

    if preprocess:
        gene_names = preprocess_gene_names(pseudobulk_series.index, dataset_name)
        pseudobulk_series.index = gene_names

    return pseudobulk_series


def normalize_cpm(counts_series):
    """
    Normalize counts to CPM (Counts Per Million).
    
    Parameters:
    -----------
    counts_series : pd.Series
        Gene counts
        
    Returns:
    --------
    cpm_series : pd.Series
        Log2-transformed CPM values
    """
    total = counts_series.sum()
    cpm = (counts_series / total) * 1e6
    cpm_log = np.log2(cpm + 1)  # Add pseudocount to avoid log(0)
    
    return cpm_log


def diagnose_gene_mismatch(short_norm, long_norm):
    """
    Diagnose why genes might not be matching between datasets.
    
    Parameters:
    -----------
    short_norm : pd.Series
    long_norm : pd.Series
    """
    print("\n" + "=" * 80)
    print("GENE NAME DIAGNOSTIC")
    print("=" * 80)
    
    short_genes = set(short_norm.index)
    long_genes = set(long_norm.index)
    
    print(f"\nShort-read genes: {len(short_genes)}")
    print(f"  Sample names: {list(short_genes)[:5]}")
    print(f"  First gene type: {type(list(short_genes)[0])}")
    
    print(f"\nLong-read genes: {len(long_genes)}")
    print(f"  Sample names: {list(long_genes)[:5]}")
    print(f"  First gene type: {type(list(long_genes)[0])}")
    
    shared = short_genes & long_genes
    print(f"\nShared genes: {len(shared)}")
    if len(shared) > 0:
        print(f"  Sample shared genes: {list(shared)[:5]}")
    
    # Check for potential issues
    print(f"\n--- Potential Issues ---")
    
    # Check if gene names are formatted differently
    short_sample = list(short_genes)[:3]
    long_sample = list(long_genes)[:3]
    
    print(f"\nChecking for common patterns:")
    print(f"  Short-read sample: {short_sample}")
    print(f"  Long-read sample: {long_sample}")
    
    # Try case-insensitive matching
    short_genes_lower = {str(g).lower() for g in short_genes}
    long_genes_lower = {str(g).lower() for g in long_genes}
    shared_lower = short_genes_lower & long_genes_lower
    print(f"\n  Case-insensitive match: {len(shared_lower)} genes")
    
    # Try without version numbers (e.g., ENS123.1 -> ENS123)
    import re
    short_genes_noversion = {re.sub(r'\.\d+$', '', str(g)) for g in short_genes}
    long_genes_noversion = {re.sub(r'\.\d+$', '', str(g)) for g in long_genes}
    shared_noversion = short_genes_noversion & long_genes_noversion
    print(f"  Match (ignoring version numbers): {len(shared_noversion)} genes")


def handle_and_visualize_duplicates(short_pb, long_pb, output_dir="./results"):
    """
    Identify, report, and handle duplicate gene names in raw pseudobulk counts.
    Must be called BEFORE normalization so duplicates are summed on raw counts.
    Aggregates duplicates using SUM (appropriate for pseudobulk counts).
    
    Parameters:
    -----------
    short_pb : pd.Series
        Raw pseudobulk counts for short-read data
    long_pb : pd.Series
        Raw pseudobulk counts for long-read data
    output_dir : str
        Directory where output files will be saved
        
    Returns:
    --------
    short_pb_agg : pd.Series
        Aggregated short-read raw counts (duplicates summed)
    long_pb_agg : pd.Series
        Aggregated long-read raw counts (duplicates summed)
    duplicates_report : dict
    """
    print("\n" + "=" * 80)
    print("STEP 3: DUPLICATE GENE HANDLING (on raw pseudobulk counts, before normalization)")
    print("=" * 80)

    # Use local aliases to avoid renaming throughout the function body
    short_norm = short_pb
    long_norm = long_pb
    
    # Identify duplicates
    short_duplicates_mask = short_norm.index.duplicated(keep=False)
    long_duplicates_mask = long_norm.index.duplicated(keep=False)
    
    short_duplicates = short_norm[short_duplicates_mask].sort_index()
    long_duplicates = long_norm[long_duplicates_mask].sort_index()
    
    short_dup_genes = short_norm.index[short_norm.index.duplicated()].unique()
    long_dup_genes = long_norm.index[long_norm.index.duplicated()].unique()
    
    print(f"\nShort-read duplicates:")
    print(f"  Unique duplicate genes: {len(short_dup_genes)}")
    if len(short_dup_genes) > 0:
        print(f"  Total duplicate entries: {short_duplicates_mask.sum()}")
        print(f"  Examples: {list(short_dup_genes[:5])}")
        print(f"\n  Duplicate entries:")
        print(short_duplicates.to_string())
    else:
        print(f"  None found")
    
    print(f"\nLong-read duplicates:")
    print(f"  Unique duplicate genes: {len(long_dup_genes)}")
    if len(long_dup_genes) > 0:
        print(f"  Total duplicate entries: {long_duplicates_mask.sum()}")
        print(f"  Examples: {list(long_dup_genes[:5])}")
        print(f"\n  Duplicate entries:")
        print(long_duplicates.to_string())
    else:
        print(f"  None found")
    
    # Aggregate duplicates by SUMMING (appropriate for counts/pseudobulk)
    def aggregate_by_sum(series):
        """Aggregate duplicate indices by taking the sum."""
        if series.index.duplicated().any():
            print(f"\n  Aggregating duplicates by SUM...")
            return series.groupby(level=0).sum()
        return series
    
    print(f"\nAggregating duplicates...")
    short_pb_agg = aggregate_by_sum(short_norm)
    long_pb_agg = aggregate_by_sum(long_norm)
    
    print(f"\nAfter aggregation:")
    print(f"  Short-read genes: {len(short_norm)} → {len(short_pb_agg)}")
    print(f"  Long-read genes: {len(long_norm)} → {len(long_pb_agg)}")
    
    duplicates_report = {
        'short_read_duplicate_genes': len(short_dup_genes),
        'short_read_duplicate_entries': short_duplicates_mask.sum(),
        'long_read_duplicate_genes': len(long_dup_genes),
        'long_read_duplicate_entries': long_duplicates_mask.sum(),
    }
    
    return short_pb_agg, long_pb_agg, duplicates_report


def _density_scatter(ax, x, y, bins=200, cmap='viridis', s=5, **kwargs):
    """Scatter plot coloured by local point density using a 2D histogram."""
    x = np.asarray(x)
    y = np.asarray(y)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    x_idx = np.clip(np.searchsorted(x_edges, x) - 1, 0, counts.shape[0] - 1)
    y_idx = np.clip(np.searchsorted(y_edges, y) - 1, 0, counts.shape[1] - 1)
    density = counts[x_idx, y_idx]
    order = density.argsort()
    sc = ax.scatter(x[order], y[order], c=density[order], s=s, cmap=cmap, **kwargs)
    sc.set_rasterized(True)
    return sc


def analysis_correlations_combined(short_norm, long_norm, output_dir="./results", min_log2cpm=0):
    """
    Combined Analysis 1 & 2: Correlation plots side-by-side.
    - Shared genes: genes expressed in both modalities
    - All expressed genes: union of genes expressed in at least one modality
    - Excludes genes with zero expression in both modalities
    
    Parameters:
    -----------
    short_norm : pd.Series
        Normalized short-read pseudobulk
    long_norm : pd.Series
        Normalized long-read pseudobulk
    output_dir : str
        Directory where output files will be saved
        
    Returns:
    --------
    results : dict
    """
    print("\n" + "=" * 80)
    print("STEP 4: CORRELATION ANALYSIS")
    print("=" * 80)
    
    # ===== SHARED GENES ANALYSIS =====
    print("\n--- Analysis 1: Shared Genes ---")
    
    # Find shared genes
    shared_genes = list(set(short_norm.index) & set(long_norm.index))
    shared_genes_sorted = sorted(shared_genes)
    
    # Handle case with no/very few shared genes
    if len(shared_genes) < 2:
        print(f"\n⚠️  WARNING: Only {len(shared_genes)} shared genes found!")
        print("\nThis likely means gene names are formatted differently.")
        print("Running diagnostic to identify the issue...\n")
        diagnose_gene_mismatch(short_norm, long_norm)
        
        # Return placeholder results
        return {
            'shared_genes': len(shared_genes),
            'shared_genes_expressed': 0,
            'pearson_r_shared': np.nan,
            'pearson_p_shared': np.nan,
            'spearman_r_shared': np.nan,
            'spearman_p_shared': np.nan,
            'all_genes_expressed': 0,
            'pearson_r_all': np.nan,
            'pearson_p_all': np.nan,
            'spearman_r_all': np.nan,
            'spearman_p_all': np.nan,
            'data': pd.DataFrame(),
            'skipped': True
        }
    
    short_shared = short_norm[shared_genes_sorted]
    long_shared = long_norm[shared_genes_sorted]
    
    # Calculate correlations for shared genes
    pearson_r_shared, pearson_p_shared = pearsonr(short_shared, long_shared)
    spearman_r_shared, spearman_p_shared = spearmanr(short_shared, long_shared)
    
    print(f"Shared genes (in both datasets): {len(shared_genes)}")
    r2_shared = pearson_r_shared ** 2
    print(f"Correlation metrics:")
    print(f"  Pearson r: {pearson_r_shared:.4f} (p-value: {pearson_p_shared:.2e})")
    print(f"  R²: {r2_shared:.4f}")
    print(f"  Spearman ρ: {spearman_r_shared:.4f} (p-value: {spearman_p_shared:.2e})")
    
    # ===== ALL EXPRESSED GENES ANALYSIS =====
    print("\n--- Analysis 2: All Expressed Genes ---")
    
    # Union of all genes (expressed in at least one modality)
    all_genes_raw = sorted(set(short_norm.index) | set(long_norm.index))
    
    # Filter to exclude genes with zero expression in both modalities
    all_genes_expressed = []
    for gene in all_genes_raw:
        short_val = short_norm.get(gene, 0)
        long_val = long_norm.get(gene, 0)
        # Keep gene if expressed in at least one modality (not both zero)
        if short_val > 0 or long_val > 0:
            all_genes_expressed.append(gene)
    
    all_genes_expressed = sorted(all_genes_expressed)
    
    # Reindex with 0 for missing genes (in log space: log2(1) = 0 means not expressed)
    short_all = short_norm.reindex(all_genes_expressed, fill_value=0)
    long_all = long_norm.reindex(all_genes_expressed, fill_value=0)
    
    # Calculate correlations for all expressed genes
    pearson_r_all, pearson_p_all = pearsonr(short_all, long_all)
    spearman_r_all, spearman_p_all = spearmanr(short_all, long_all)
    
    n_genes_both_zero = len(all_genes_raw) - len(all_genes_expressed)
    print(f"Total genes (union): {len(all_genes_raw)}")
    print(f"Genes expressed in at least one modality: {len(all_genes_expressed)}")
    print(f"Genes with zero expression in both: {n_genes_both_zero}")
    r2_all = pearson_r_all ** 2
    print(f"Correlation metrics:")
    print(f"  Pearson r: {pearson_r_all:.4f} (p-value: {pearson_p_all:.2e})")
    print(f"  R²: {r2_all:.4f}")
    print(f"  Spearman ρ: {spearman_r_all:.4f} (p-value: {spearman_p_all:.2e})")
    
    # ===== COMBINED PLOT =====
    print("\nGenerating combined correlation plots...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Left plot: Shared genes
    sc0 = _density_scatter(axes[0], short_shared, long_shared, alpha=0.8)
    plt.colorbar(sc0, ax=axes[0], label='Density')
    axes[0].set_xlabel('Short-read (log2 CPM)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Long-read (log2 CPM)', fontsize=12, fontweight='bold')
    axes[0].set_title(f'Shared Genes Correlation\n(n={len(shared_genes)})', 
                     fontsize=14, fontweight='bold')
    
    # Add correlation text for shared genes
    textstr_shared = f'Pearson r = {pearson_r_shared:.4f}\nR² = {r2_shared:.4f}\nSpearman ρ = {spearman_r_shared:.4f}'
    axes[0].text(0.05, 0.95, textstr_shared, transform=axes[0].transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add diagonal line for shared genes
    lims_shared = [
        np.min([axes[0].get_xlim(), axes[0].get_ylim()]),
        np.max([axes[0].get_xlim(), axes[0].get_ylim()]),
    ]
    axes[0].plot(lims_shared, lims_shared, 'k--', alpha=0.3, zorder=0)
    axes[0].set_xlim(lims_shared)
    axes[0].set_ylim(lims_shared)
    axes[0].set_aspect('equal')
    axes[0].grid(alpha=0.3)
    
    # Right plot: All expressed genes
    sc1 = _density_scatter(axes[1], short_all, long_all, alpha=0.8)
    plt.colorbar(sc1, ax=axes[1], label='Density')
    axes[1].set_xlabel('Short-read (log2 CPM)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Long-read (log2 CPM)', fontsize=12, fontweight='bold')
    axes[1].set_title(f'All Expressed Genes Correlation\n(n={len(all_genes_expressed)})', 
                     fontsize=14, fontweight='bold')
    
    # Add correlation text for all expressed genes
    textstr_all = f'Pearson r = {pearson_r_all:.4f}\nR² = {r2_all:.4f}\nSpearman ρ = {spearman_r_all:.4f}'
    axes[1].text(0.05, 0.95, textstr_all, transform=axes[1].transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add diagonal line for all expressed genes
    lims_all = [
        np.min([axes[1].get_xlim(), axes[1].get_ylim()]),
        np.max([axes[1].get_xlim(), axes[1].get_ylim()]),
    ]
    axes[1].plot(lims_all, lims_all, 'k--', alpha=0.3, zorder=0)
    axes[1].set_xlim(lims_all)
    axes[1].set_ylim(lims_all)
    axes[1].set_aspect('equal')
    axes[1].grid(alpha=0.3)
    
    plt.suptitle('Correlation Analysis: Short-read vs Long-read', 
                fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    import os
    output_path = os.path.join(output_dir, 'correlation_plots_combined.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nCombined plot saved: {output_path}")
    plt.close()
    
    return {
        'shared_genes': len(shared_genes),
        'shared_genes_expressed': len(shared_genes),
        'pearson_r_shared': pearson_r_shared,
        'pearson_p_shared': pearson_p_shared,
        'r2_shared': r2_shared,
        'spearman_r_shared': spearman_r_shared,
        'spearman_p_shared': spearman_p_shared,
        'all_genes_total': len(all_genes_raw),
        'all_genes_expressed': len(all_genes_expressed),
        'all_genes_both_zero': n_genes_both_zero,
        'pearson_r_all': pearson_r_all,
        'pearson_p_all': pearson_p_all,
        'r2_all': r2_all,
        'spearman_r_all': spearman_r_all,
        'spearman_p_all': spearman_p_all,
        'data_shared': pd.DataFrame({
            'short_read': short_shared,
            'long_read': long_shared
        }),
        'data_all': pd.DataFrame({
            'short_read': short_all,
            'long_read': long_all
        }),
        'skipped': False
    }


def analysis_venn_diagram(short_norm, long_norm, output_dir="./results"):
    """
    Analysis 3: Venn diagram of gene overlap between modalities.
    Excludes genes with zero expression in both modalities.

    Parameters:
    -----------
    short_norm : pd.Series
        Normalized short-read pseudobulk
    long_norm : pd.Series
        Normalized long-read pseudobulk
    output_dir : str
        Directory where output files will be saved

    Returns:
    --------
    results : dict
    """
    print("\n" + "=" * 80)
    print("STEP 5: GENE OVERLAP VENN DIAGRAM")
    print("=" * 80)

    # Handle duplicates - get unique gene names
    short_genes_all = set(short_norm.index.unique())
    long_genes_all = set(long_norm.index.unique())

    # Filter to keep only genes expressed in at least one modality
    short_genes = set()
    for gene in short_genes_all:
        if short_norm[gene] > 0:
            short_genes.add(gene)

    long_genes = set()
    for gene in long_genes_all:
        if long_norm[gene] > 0:
            long_genes.add(gene)

    shared = short_genes & long_genes
    only_short = short_genes - long_genes
    only_long = long_genes - short_genes

    # Count genes with zero expression in both
    all_genes_total = short_genes_all | long_genes_all
    genes_both_zero = all_genes_total - short_genes - long_genes - shared

    print(f"\nGene overlap statistics (expressed genes only):")
    print(f"  Short-read only: {len(only_short)}")
    print(f"  Long-read only: {len(only_long)}")
    print(f"  Shared (in both): {len(shared)}")
    print(f"  Total expressed: {len(short_genes | long_genes)}")
    print(f"\nGenes with zero expression in both modalities: {len(genes_both_zero)}")

    # Plot Venn diagram
    all_expressed = short_genes | long_genes
    overlap_pct = 100 * len(shared) / len(all_expressed) if all_expressed else 0
    total_expressed = len(all_expressed)

    fig, ax = plt.subplots(figsize=(8, 8))
    venn = venn2([short_genes, long_genes],
                  set_labels=('Short-read', 'Long-read'),
                  ax=ax)

    # Modify text labels to show percentages
    # subset_labels contains [only_short, only_long, shared]
    counts = [len(only_short), len(only_long), len(shared)]
    for i, text in enumerate(venn.subset_labels):
        if text is not None:
            count = counts[i]
            percentage = 100 * count / total_expressed if total_expressed > 0 else 0
            text.set_text(f'{count}\n({percentage:.1f}%)')

    plt.title(
        f'Gene Detection Overlap Between Modalities\n(Expressed genes only)\n'
        f'Overlap: {overlap_pct:.1f}% of total expressed genes',
        fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    import os
    output_path = os.path.join(output_dir, 'analysis3_venn_diagram.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved: {output_path}")
    plt.close()

    # Save gene lists to txt files (one gene per line, sorted)
    gene_lists = {
        'genes_only_short_read.txt': sorted(only_short),
        'genes_only_long_read.txt':  sorted(only_long),
        'genes_shared.txt':          sorted(shared),
    }
    for filename, gene_list in gene_lists.items():
        txt_path = os.path.join(output_dir, filename)
        with open(txt_path, 'w') as f:
            f.write('\n'.join(gene_list))
        print(f"  Gene list saved: {txt_path} ({len(gene_list)} genes)")

    return {
        'short_read_only': len(only_short),
        'long_read_only': len(only_long),
        'shared': len(shared),
        'total_expressed': len(short_genes | long_genes),
        'total_unique': len(all_genes_total),
        'genes_both_zero': len(genes_both_zero),
        'short_genes': only_short,
        'long_genes': only_long,
        'shared_genes': shared
    }


def generate_summary_report(barcode_report, combined_analysis, analysis3, output_dir="./results"):
    """
    Generate a comprehensive summary report.
    
    Parameters:
    -----------
    barcode_report : dict
    combined_analysis : dict
    analysis3 : dict
    output_dir : str
        Directory where output files will be saved
        
    Returns:
    --------
    report_df : pd.DataFrame
    """
    print("\n" + "=" * 80)
    print("STEP 6: SUMMARY REPORT")
    print("=" * 80)
    
    # Format values with handling for skipped analyses
    def format_value(val):
        if isinstance(val, float) and np.isnan(val):
            return "N/A (insufficient data)"
        elif isinstance(val, float):
            return f"{val:.4f}"
        return str(val)
    
    not_skipped = not combined_analysis.get('skipped', False)
    shared_pearson  = format_value(combined_analysis['pearson_r_shared'])  if not_skipped else "N/A (insufficient data)"
    shared_r2       = format_value(combined_analysis['r2_shared'])          if not_skipped else "N/A (insufficient data)"
    shared_spearman = format_value(combined_analysis['spearman_r_shared']) if not_skipped else "N/A (insufficient data)"
    all_pearson     = format_value(combined_analysis['pearson_r_all'])      if not_skipped else "N/A (insufficient data)"
    all_r2          = format_value(combined_analysis['r2_all'])              if not_skipped else "N/A (insufficient data)"
    all_spearman    = format_value(combined_analysis['spearman_r_all'])     if not_skipped else "N/A (insufficient data)"
    
    report_dict = {
        'Metric': [
            'Short-read cells (original)',
            'Long-read cells (original)',
            'Common barcodes',
            'Overlap (% short-read)',
            'Overlap (% long-read)',
            '',
            'SHARED GENES ANALYSIS',
            'Shared genes (detected in both)',
            'Pearson r (shared)',
            'R² (shared)',
            'Spearman ρ (shared)',
            '',
            'ALL EXPRESSED GENES ANALYSIS',
            'Total genes (union)',
            'Genes with zero expression in both',
            'Genes expressed in at least one',
            'Pearson r (all expressed)',
            'R² (all expressed)',
            'Spearman ρ (all expressed)',
            '',
            'GENE OVERLAP (VENN DIAGRAM)',
            'Short-read only genes',
            'Long-read only genes',
            'Shared genes (expressed)',
            'Total unique genes'
        ],
        'Value': [
            f"{barcode_report['short_read_original']}",
            f"{barcode_report['long_read_original']}",
            f"{barcode_report['common_barcodes']}",
            f"{barcode_report['overlap_pct_short']:.1f}%",
            f"{barcode_report['overlap_pct_long']:.1f}%",
            '',
            '',
            f"{combined_analysis['shared_genes']}",
            shared_pearson,
            shared_r2,
            shared_spearman,
            '',
            '',
            f"{combined_analysis['all_genes_total']}",
            f"{combined_analysis['all_genes_both_zero']}",
            f"{combined_analysis['all_genes_expressed']}",
            all_pearson,
            all_r2,
            all_spearman,
            '',
            '',
            f"{analysis3['short_read_only']}",
            f"{analysis3['long_read_only']}",
            f"{analysis3['shared']}",
            f"{analysis3['total_unique']}"
        ]
    }
    
    report_df = pd.DataFrame(report_dict)
    
    print("\n")
    print(report_df.to_string(index=False))
    
    # Save report
    import os
    output_path = os.path.join(output_dir, 'summary_report.csv')
    report_df.to_csv(output_path, index=False)
    print(f"\nReport saved: {output_path}")
    
    return report_df


def main(short_read_path, long_read_path, output_dir="./results",
         barcode_sep_short='-', barcode_sep_long=None):
    """
    Main analysis pipeline - BARCODE-FIRST APPROACH.

    Workflow:
    1. Load data
    2. Filter to shared barcodes FIRST
    3. Compute pseudobulk on filtered data
    4. Handle duplicates (sum duplicate gene names on raw counts)
    5. Normalize to CPM (on clean, deduplicated counts)
    6. Perform correlation analyses
    7. Generate report

    Parameters:
    -----------
    short_read_path : str
        Path to short-read matrix
    long_read_path : str
        Path to long-read matrix
    output_dir : str
        Directory where output files will be saved (default: "./results")
    barcode_sep_short : str
        Separator used in short-read obs_names (default: '-')
    barcode_sep_long : str or None
        Separator used in long-read obs_names. If None, uses barcode_sep_short.
    """
    if barcode_sep_long is None:
        barcode_sep_long = barcode_sep_short
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 78 + "║")
    print("║" + "SINGLE-CELL RNA-SEQ MODALITY COMPARISON (BARCODE-FIRST APPROACH)".center(78) + "║")
    print("║" + "Short-read vs. Long-read Analysis".center(78) + "║")
    print("║" + " " * 78 + "║")
    print("╚" + "=" * 78 + "╝")
    
    # Step 1: Load data
    adata_short, adata_long = load_and_prepare_data(short_read_path, long_read_path)
    
    # Step 2: Filter to shared barcodes FIRST
    adata_short_filtered, adata_long_filtered, barcode_report = filter_by_shared_barcodes(
        adata_short, adata_long,
        barcode_sep_short=barcode_sep_short,
        barcode_sep_long=barcode_sep_long,
    )

    # Step 2b: Filter out zero-count genes per modality
    adata_short_filtered, adata_long_filtered, gene_filter_report = filter_zero_count_genes(
        adata_short_filtered, adata_long_filtered
    )

    # Step 3: Pseudobulk on FILTERED data
    print("\n" + "=" * 80)
    print("STEP 3: PSEUDOBULK AGGREGATION (ON FILTERED DATA)")
    print("=" * 80)
    short_pb = compute_pseudobulk(adata_short_filtered, preprocess=True, dataset_name="Short-read")
    long_pb = compute_pseudobulk(adata_long_filtered, preprocess=False, dataset_name="Long-read")
    print(f"\nShort-read pseudobulk genes: {len(short_pb)}")
    print(f"Long-read pseudobulk genes: {len(long_pb)}")
    
    # Step 3: Handle duplicates on raw pseudobulk counts (BEFORE normalization)
    short_pb_agg, long_pb_agg, duplicates_report = handle_and_visualize_duplicates(
        short_pb, long_pb, output_dir
    )

    # Step 4: Normalize to CPM (on deduplicated counts)
    print("\n" + "=" * 80)
    print("STEP 4: CPM NORMALIZATION (on deduplicated counts)")
    print("=" * 80)
    short_norm = normalize_cpm(short_pb_agg)
    long_norm = normalize_cpm(long_pb_agg)
    print(f"\nNormalization complete")
    print(f"Short-read total counts (before normalization): {short_pb_agg.sum():.0f}")
    print(f"Long-read total counts (before normalization): {long_pb_agg.sum():.0f}")

    # Step 5: Correlation analyses
    combined_analysis = analysis_correlations_combined(short_norm, long_norm, output_dir)
    analysis3 = analysis_venn_diagram(short_norm, long_norm, output_dir)
    
    # Step 6: Summary report
    report = generate_summary_report(barcode_report, combined_analysis, analysis3, output_dir)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nOutput files saved to: {output_dir}")
    print("  - summary_report.csv")
    print("  - correlation_plots_combined.png")
    print("  - analysis3_venn_diagram.png")
    print("  - genes_only_short_read.txt")
    print("  - genes_only_long_read.txt")
    print("  - genes_shared.txt")
    
    return {
        'adata_short_full': adata_short,
        'adata_long_full': adata_long,
        'adata_short_filtered': adata_short_filtered,
        'adata_long_filtered': adata_long_filtered,
        'short_pb': short_pb_agg,
        'long_pb': long_pb_agg,
        'short_norm': short_norm,
        'long_norm': long_norm,
        'barcode_report': barcode_report,
        'duplicates_report': duplicates_report,
        'combined_analysis': combined_analysis,
        'analysis3': analysis3,
        'summary_report': report
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare Short-read and Long-read scRNA-seq datasets. "
                    "Supported formats: .h5ad, .h5, or .csv"
    )
    parser.add_argument("--short", required=True, metavar="PATH",
                        help="Path to short-read count matrix")
    parser.add_argument("--long", required=True, metavar="PATH",
                        help="Path to long-read count matrix")
    parser.add_argument("--output-dir", default="./short_long_output", metavar="DIR",
                        help="Output directory (default: ./short_long_output)")
    parser.add_argument("--barcode-sep-short", default="-", metavar="SEP",
                        help="Separator for short-read obs_names (default: -)")
    parser.add_argument("--barcode-sep-long", default=None, metavar="SEP",
                        help="Separator for long-read obs_names (default: same as --barcode-sep-short)")
    args = parser.parse_args()

    results = main(args.short, args.long, output_dir=args.output_dir,
                   barcode_sep_short=args.barcode_sep_short,
                   barcode_sep_long=args.barcode_sep_long)
