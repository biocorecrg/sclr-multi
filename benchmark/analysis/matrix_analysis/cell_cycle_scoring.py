#!/usr/bin/env python3
"""
Cell Cycle Scoring and UMAP Visualization
===========================================

Scores cells based on cell cycle genes and visualizes the scores on UMAP plots.
Uses scanpy's built-in cell cycle gene sets (S phase and G2/M phase genes).

Usage
-----
python cell_cycle_scoring.py \
    --input matrix1.h5ad matrix2.h5ad \
    --output_dir ./cell_cycle_results

python cell_cycle_scoring.py \
    --input short_read.h5ad long_read.h5ad \
    --output_dir ./results
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("cell_cycle_scoring")


def load_file(filepath):
    """Load various file formats into AnnData object."""
    filepath = Path(filepath)

    if filepath.suffix == ".h5ad":
        return sc.read_h5ad(filepath)
    elif filepath.suffix == ".h5":
        return sc.read_10x_h5(filepath)
    elif filepath.suffix == ".csv":
        df = pd.read_csv(filepath, index_col=0)
        return sc.AnnData(
            X=df.T.values, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index)
        )
    else:
        raise ValueError(f"Unsupported format: {filepath.suffix}. Use .h5ad, .h5, or .csv")


def score_cell_cycle(adata, sample_name=None):
    """
    Score cells for cell cycle phase using mouse cell cycle genes.
    Mouse genes are used as a universal mammalian set across all organisms.

    Parameters
    ----------
    adata : AnnData
        Input count matrix. Must have gene names in .var_names.
    sample_name : str, optional
        Name of the sample for logging.

    Returns
    -------
    adata : AnnData
        Updated matrix with S_score, G2M_score, and phase columns in .obs.
    gene_info : dict
        Information about genes used for scoring.
    """
    if sample_name is None:
        sample_name = "Unknown"

    logger.info("Scoring cell cycle for %s ...", sample_name)

    # Log original shape
    logger.info("  Shape before processing: %d cells × %d genes", adata.n_obs, adata.n_vars)

    # Check for raw count layers
    if adata.layers:
        logger.info("  Available layers: %s", list(adata.layers.keys()))
        raw_layer = None
        for layer_name in ['raw_counts', 'counts', 'X_raw', 'raw']:
            if layer_name in adata.layers:
                raw_layer = layer_name
                break
        if raw_layer:
            logger.info("  Using raw counts from layer: %s", raw_layer)
            adata.X = adata.layers[raw_layer]
        elif adata.X.max() < 100:
            logger.warning("  Warning: Expression values seem normalized. Best results with raw counts.")
    elif adata.X.max() < 100:
        logger.warning("  Warning: Expression values seem normalized. Best results with raw counts.")

    # Mouse cell cycle genes (Ensembl IDs)
    # S phase and G2/M phase genes based on Tirosh et al. (2016)
    # "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq" (Science)
    # Modified: 2017-09-13
    s_genes = [
        "ENSMUSG00000000028", "ENSMUSG00000001228", "ENSMUSG00000002870", "ENSMUSG00000004642", "ENSMUSG00000005410",
        "ENSMUSG00000006678", "ENSMUSG00000006715", "ENSMUSG00000017499", "ENSMUSG00000020649", "ENSMUSG00000022360",
        "ENSMUSG00000022422", "ENSMUSG00000022673", "ENSMUSG00000022945", "ENSMUSG00000023104", "ENSMUSG00000024151",
        "ENSMUSG00000024742", "ENSMUSG00000025001", "ENSMUSG00000025395", "ENSMUSG00000025747", "ENSMUSG00000026355",
        "ENSMUSG00000027242", "ENSMUSG00000027323", "ENSMUSG00000027342", "ENSMUSG00000028212", "ENSMUSG00000028282",
        "ENSMUSG00000028560", "ENSMUSG00000028693", "ENSMUSG00000028884", "ENSMUSG00000029591", "ENSMUSG00000030346",
        "ENSMUSG00000030528", "ENSMUSG00000030726", "ENSMUSG00000030978", "ENSMUSG00000031629", "ENSMUSG00000031821",
        "ENSMUSG00000032397", "ENSMUSG00000034329", "ENSMUSG00000037474", "ENSMUSG00000039748", "ENSMUSG00000041712",
        "ENSMUSG00000042489", "ENSMUSG00000046179", "ENSMUSG00000055612",
        ]

    g2m_genes = [
        "ENSMUSG00000001403", "ENSMUSG00000004880", "ENSMUSG00000005698", "ENSMUSG00000006398", "ENSMUSG00000009575",
        "ENSMUSG00000012443", "ENSMUSG00000015749", "ENSMUSG00000017716", "ENSMUSG00000019942", "ENSMUSG00000019961",
        "ENSMUSG00000020330", "ENSMUSG00000020737", "ENSMUSG00000020808", "ENSMUSG00000020897", "ENSMUSG00000020914",
        "ENSMUSG00000022385", "ENSMUSG00000022391", "ENSMUSG00000023505", "ENSMUSG00000024056", "ENSMUSG00000024795",
        "ENSMUSG00000026605", "ENSMUSG00000026622", "ENSMUSG00000026683", "ENSMUSG00000027306", "ENSMUSG00000027379",
        "ENSMUSG00000027469", "ENSMUSG00000027496", "ENSMUSG00000027699", "ENSMUSG00000028044", "ENSMUSG00000028678",
        "ENSMUSG00000028873", "ENSMUSG00000029177", "ENSMUSG00000031004", "ENSMUSG00000032218", "ENSMUSG00000032254",
        "ENSMUSG00000034349", "ENSMUSG00000035293", "ENSMUSG00000036752", "ENSMUSG00000036777", "ENSMUSG00000037313",
        "ENSMUSG00000037544", "ENSMUSG00000037725", "ENSMUSG00000038252", "ENSMUSG00000038379", "ENSMUSG00000040549",
        "ENSMUSG00000044201", "ENSMUSG00000044783", "ENSMUSG00000045328", "ENSMUSG00000048327", "ENSMUSG00000048922",
        "ENSMUSG00000054717", "ENSMUSG00000062248", "ENSMUSG00000068744", "ENSMUSG00000074802",
        ]

    logger.info("  Cell cycle gene sets (mouse + human gene names for cross-organism compatibility):")
    logger.info("    S phase genes: %d total", len(set(s_genes)))
    logger.info("    G2/M phase genes: %d total", len(set(g2m_genes)))

    # Use gene_name column if it exists, otherwise fall back to var_names
    if "gene_id" in adata.var.columns:
        logger.info("  Using gene_id column for cell cycle gene matching")
        gene_names_set = set(adata.var["gene_id"].dropna())
    else:
        logger.info("  Using var_names for cell cycle gene matching")
        gene_names_set = set(adata.var_names)

    # Find which genes are in the dataset
    s_genes_found = [g for g in s_genes if g in gene_names_set]
    g2m_genes_found = [g for g in g2m_genes if g in gene_names_set]

    # Remove duplicates while preserving order
    s_genes_found = list(dict.fromkeys(s_genes_found))
    g2m_genes_found = list(dict.fromkeys(g2m_genes_found))

    logger.info("  Genes found in matrix:")
    logger.info("    S phase: %d genes found", len(s_genes_found))
    logger.info("    G2/M phase: %d genes found", len(g2m_genes_found))

    if len(s_genes_found) == 0 or len(g2m_genes_found) == 0:
        logger.error("  ERROR: Not enough cell cycle genes found in matrix!")
        logger.error("  Please check that gene names match standard symbols (e.g., 'MKI67' or 'Mki67')")
        raise ValueError("Insufficient cell cycle genes found for scoring")

    # Create fresh AnnData with gene_id as var_names for scoring
    if "gene_id" in adata.var.columns:
        adata_for_scoring = sc.AnnData(
            X=adata.X,
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=adata.var["gene_id"].astype(str).values)
        )
        adata_for_scoring.var_names_make_unique()

        sc.pp.normalize_total(adata_for_scoring, target_sum=1e4)
        sc.pp.log1p(adata_for_scoring)

        logger.info("  Created temporary AnnData with gene_id as var_names")
        logger.info("  Using %d S-phase genes and %d G2M genes",
                    len(s_genes_found), len(g2m_genes_found))

        # Score using scanpy's cell cycle function on the fresh object
        sc.tl.score_genes_cell_cycle(adata_for_scoring, s_genes=s_genes_found, g2m_genes=g2m_genes_found)

        # Copy scores back to original adata
        adata.obs["S_score"] = adata_for_scoring.obs["S_score"].values
        adata.obs["G2M_score"] = adata_for_scoring.obs["G2M_score"].values

        # Apply absolute threshold: cells only labeled S or G2M if score > median + 2.5*SD
        s_median = adata.obs["S_score"].median()
        s_std = adata.obs["S_score"].std()
        g2m_median = adata.obs["G2M_score"].median()
        g2m_std = adata.obs["G2M_score"].std()

        # Establish threshold, if median lower than 0, threshold is 0.
        if s_median < 0:
            s_threshold = 0 + s_std
        else:
            s_threshold = s_median + s_std
        
        if g2m_median < 0:
            g2m_threshold = 0 + g2m_std
        else:
            g2m_threshold = g2m_median + g2m_std

        logger.info("  Applying absolute thresholds:")
        logger.info("    S_score threshold:   %.3f", s_threshold)
        logger.info("    G2M_score threshold: %.3f", g2m_threshold)

        # Assign phase based on absolute thresholds
        phase = []
        for idx in adata.obs.index:
            s_score = adata.obs.loc[idx, "S_score"]
            g2m_score = adata.obs.loc[idx, "G2M_score"]

            if s_score > s_threshold and s_score > g2m_score:
                phase.append("S")
            elif g2m_score > g2m_threshold and g2m_score > s_score:
                phase.append("G2M")
            else:
                phase.append("G1")

        adata.obs["phase"] = phase

        logger.info("  Cell cycle phase assignment complete (with absolute thresholds)")
    else:
        raise ValueError("gene_id column not found in .var")

    logger.info("  Cell cycle scoring complete")
    logger.info("    Mean S score: %.3f ± %.3f", adata.obs["S_score"].mean(),
                adata.obs["S_score"].std())
    logger.info("    Mean G2M score: %.3f ± %.3f", adata.obs["G2M_score"].mean(),
                adata.obs["G2M_score"].std())
    logger.info("    Phases: %s", dict(adata.obs["phase"].value_counts()))

    gene_info = {
        "s_genes_found": s_genes_found,
        "g2m_genes_found": g2m_genes_found,
        "s_genes_total": len(set(s_genes)),
        "g2m_genes_total": len(set(g2m_genes)),
    }

    return adata, gene_info


def compute_umap_if_needed(adata, sample_name=None):
    """
    Compute PCA and UMAP if not already present.

    Parameters
    ----------
    adata : AnnData
        Input matrix.
    sample_name : str, optional
        Name of the sample for logging.

    Returns
    -------
    adata : AnnData
        Updated matrix with PCA and UMAP computed.
    """
    if sample_name is None:
        sample_name = "Unknown"

    if "X_pca" not in adata.obsm:
        logger.info("  Computing PCA for %s ...", sample_name)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

    if "X_umap" not in adata.obsm:
        logger.info("  Computing UMAP for %s ...", sample_name)
        sc.tl.umap(adata)

    return adata


def plot_cell_cycle_umap(adata, sample_name, output_path):
    """
    Create UMAP plots colored by cell cycle scores and assigned phase.

    Parameters
    ----------
    adata : AnnData
        Scored matrix with UMAP coordinates.
    sample_name : str
        Name of the sample for plot title.
    output_path : Path
        Directory to save plots.
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Get matching color scale for S and G2M scores
    score_min = min(adata.obs["S_score"].min(), adata.obs["G2M_score"].min())
    score_max = max(adata.obs["S_score"].max(), adata.obs["G2M_score"].max())

    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Cell Cycle Scoring: {sample_name}", fontsize=14, fontweight="bold")

    # Plot 1: S phase score
    sc.pl.umap(adata, color="S_score", ax=axes[0], show=False, cmap="Reds",
               vmin=score_min, vmax=score_max)
    axes[0].set_title("S Phase Score")

    # Plot 2: G2M phase score
    sc.pl.umap(adata, color="G2M_score", ax=axes[1], show=False, cmap="Reds",
               vmin=score_min, vmax=score_max)
    axes[1].set_title("G2/M Phase Score")

    # Plot 3: Assigned phase
    phase_colors = {"G1": "#dadada", "S": "#5a53ce", "G2M": "#7fb228"}
    sc.pl.umap(adata, color="phase", ax=axes[2], show=False, palette=phase_colors)
    axes[2].set_title("Assigned Cell Cycle Phase")

    plt.tight_layout()
    output_file = output_path / f"{sample_name}_cell_cycle_umap.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    logger.info("  Saved plot: %s", output_file)
    plt.close()


def save_gene_info(gene_info_dict, sample_names, output_path):
    """
    Save detailed gene information to file.

    Parameters
    ----------
    gene_info_dict : dict
        Dictionary mapping sample names to gene info.
    sample_names : list
        List of sample names.
    output_path : Path
        Directory to save files.
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save summary as text
    summary_file = output_path / "cell_cycle_genes_summary.txt"
    with open(summary_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("CELL CYCLE GENES USED FOR SCORING\n")
        f.write("=" * 80 + "\n\n")

        for sample_name in sample_names:
            if sample_name not in gene_info_dict:
                continue

            info = gene_info_dict[sample_name]
            f.write(f"\nSample: {sample_name}\n")
            f.write("-" * 80 + "\n")
            f.write(f"\nS PHASE GENES ({len(info['s_genes_found'])} / {info['s_genes_total']}):\n")
            f.write(", ".join(sorted(info['s_genes_found'])) + "\n")
            f.write(f"\nG2/M PHASE GENES ({len(info['g2m_genes_found'])} / {info['g2m_genes_total']}):\n")
            f.write(", ".join(sorted(info['g2m_genes_found'])) + "\n")
            f.write("\n")

    logger.info("Saved gene information: %s", summary_file)

    # Save as CSV for easy access
    csv_file = output_path / "cell_cycle_genes.csv"
    rows = []
    for sample_name in sample_names:
        if sample_name not in gene_info_dict:
            continue

        info = gene_info_dict[sample_name]
        max_len = max(len(info['s_genes_found']), len(info['g2m_genes_found']))

        for i in range(max_len):
            s_gene = info['s_genes_found'][i] if i < len(info['s_genes_found']) else ""
            g2m_gene = info['g2m_genes_found'][i] if i < len(info['g2m_genes_found']) else ""
            rows.append({"Sample": sample_name, "S_Phase_Gene": s_gene, "G2M_Phase_Gene": g2m_gene})

    pd.DataFrame(rows).to_csv(csv_file, index=False)
    logger.info("Saved gene list: %s", csv_file)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="cell_cycle_scoring",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Path(s) to input count matrices (.h5ad, .h5, .csv)",
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Directory to save output plots and gene lists",
    )

    args = parser.parse_args(argv)

    logger.info("=" * 80)
    logger.info("CELL CYCLE SCORING AND VISUALIZATION")
    logger.info("=" * 80)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_info_dict = {}
    sample_names = []

    # Process each input file
    for input_file in args.input:
        input_path = Path(input_file)
        sample_name = input_path.stem  # Filename without extension

        logger.info("\n[%s] Loading data...", sample_name)
        try:
            adata = load_file(input_path)
            logger.info("  Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)
        except Exception as e:
            logger.error("  ERROR loading file: %s", e)
            continue

        try:
            # Score cell cycle
            adata, gene_info = score_cell_cycle(adata, sample_name)
            gene_info_dict[sample_name] = gene_info
            sample_names.append(sample_name)

            # Compute dimensionality reduction
            logger.info("[%s] Computing UMAP...", sample_name)
            adata = compute_umap_if_needed(adata, sample_name)

            # Plot results
            logger.info("[%s] Creating plots...", sample_name)
            plot_cell_cycle_umap(adata, sample_name, output_dir)

            # Save annotated matrix
            logger.info("[%s] Saving annotated matrix...", sample_name)
            output_h5ad = output_dir / f"{sample_name}_cell_cycle_scored.h5ad"
            adata.write(output_h5ad)
            logger.info("  Saved: %s", output_h5ad)

        except Exception as e:
            logger.error("  ERROR processing %s: %s", sample_name, e)
            continue

    # Save gene information
    if gene_info_dict:
        logger.info("\n" + "=" * 80)
        logger.info("SAVING GENE INFORMATION")
        logger.info("=" * 80)
        save_gene_info(gene_info_dict, sample_names, output_dir)

        logger.info("\n" + "=" * 80)
        logger.info("COMPLETE!")
        logger.info("=" * 80)
        logger.info("Output directory: %s", output_dir)
        return 0
    else:
        logger.error("No samples processed successfully")
        return 1


if __name__ == "__main__":
    sys.exit(main())
