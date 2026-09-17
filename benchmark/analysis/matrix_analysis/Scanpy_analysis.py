#!/usr/bin/env python3
"""
Single-cell RNA-seq analysis pipeline.

Subcommands (mutually exclusive):
  annotate_features      Annotate features from input matrices
  filter                 Filter cells and features
  cluster                Cluster cells
  diffexp                Differential expression analysis
  merge                  Merge N input h5ad matrixes

If no subcommand is provided, all steps run in order:
  annotate → filter → cluster → diffexp
"""

__version__ = "1.0.0"

import argparse
import logging
import sys
from pathlib import Path
import os
from scipy.stats import median_abs_deviation
from scipy import sparse
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import anndata as ad
import pybiomart as bm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
from sklearn.metrics import silhouette_score
from itertools import combinations
from pybiomart import Server
import math 
import decoupler as dc
import re
from dataclasses import dataclass
import pertpy as pt

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("Scanpy-based single cell analysis")


# ---------------------------------------------------------------------------
# Additional functions
# ---------------------------------------------------------------------------

def features_preprocessing(raw_matrix, genomes, feature_type, condition, replicate):

    #Extract feature information:
    features = raw_matrix.var.copy()

    # Check if features already have the required columns (from cell ranger or similar)
    has_genome_col = 'genome' in features.columns
    has_gene_ids_col = 'gene_ids' in features.columns

    if has_genome_col and has_gene_ids_col:
        # Features already have proper structure; rename columns for consistency
        features['organism'] = features['genome']
        features['features_id'] = features['gene_ids']
        logger.info("Features already have genome and gene_ids columns; using them directly.")
    else:
        # Parse from first column (expected format: organism-feature_id or organism_feature_id)
        split_result = features.iloc[:, 0].str.split(r'[-_]+', n=1, expand=True)

        if split_result.shape[1] == 2:
            features[['organism', 'features_id']] = split_result
        elif split_result.shape[1] == 1:
            features['organism'] = features.iloc[:, 0].apply(lambda x: genomes[0] if genomes else 'unknown')
            features['features_id'] = split_result[0]
            logger.warning("Could not parse organism from feature names; using first genome as default.")
        else:
            raise ValueError(f"Unexpected split result with {split_result.shape[1]} columns")

    _BIOMART_HOSTS = [
        "oct2024.archive.ensembl.org",
        "www.ensembl.org",
        "apr2024.archive.ensembl.org",
    ]

    def _connect_biomart():
        for host in _BIOMART_HOSTS:
            try:
                s = Server(host=host)
                _ = s.marts['ENSEMBL_MART_ENSEMBL']
                logger.info("BioMart: connected to %s", host)
                return s
            except Exception as e:
                logger.warning("BioMart: %s failed (%s), trying next host…", host, e)
        raise RuntimeError("All BioMart hosts failed. Check network or Ensembl status.")

    def _fetch_organism_features(organism, feature_type):
        """Fetch features for an organism, retrying on transient errors."""
        for host in _BIOMART_HOSTS:
            try:
                server = Server(host=host)
                mart = server.marts['ENSEMBL_MART_ENSEMBL']
                organism_db = mart.datasets[organism]

                if feature_type == 'transcript':
                    organism_features = organism_db.query(attributes=["ensembl_transcript_id", "external_transcript_name", "ensembl_gene_id", "external_gene_name", "chromosome_name"]).rename(
                        columns={"Transcript stable ID": "transcript_id", "Transcript name": "transcript_name", "Gene stable ID": "gene_id", "Gene name": "gene_name", "Chromosome/scaffold name": "chrom"})
                else:
                    organism_features = organism_db.query(attributes=["ensembl_gene_id", "external_gene_name", "chromosome_name"]).rename(
                        columns={"Gene stable ID": "gene_id", "Gene name": "gene_name", "Chromosome/scaffold name": "chrom"}
                    )
                logger.info("Successfully fetched features for %s from %s", organism, host)
                return organism_features
            except Exception as e:
                logger.warning("Failed to fetch %s from %s (%s), trying next host…", organism, host, e)
        raise RuntimeError(f"Could not fetch features for {organism} from any BioMart host.")

    initial_genome = True
    for organism in genomes:
        organism_features = _fetch_organism_features(organism, feature_type)

        #Concatenate results if needed:
        if initial_genome:
            all_features = organism_features
            initial_genome = False
        else:
            all_features = pd.concat([all_features, organism_features])
    #Add information to the AnnData's features(var):
    if feature_type == 'transcript':
        raw_matrix.var = features.merge(all_features, how="left", right_on=("transcript_id"), left_on=('features_id'))
    else:
        raw_matrix.var = features.merge(all_features, how="left", right_on=("gene_id"), left_on=('features_id'))
    
    # For long read datasets, will be features, for cell ranger h5, it will be gene_ids:
    index_col = 'features' if 'features' in raw_matrix.var.columns else 'gene_ids'
    raw_matrix.var.set_index(index_col, inplace=True)
    raw_matrix.var.index.name = None
    raw_matrix.var_names_make_unique()
    raw_matrix.var["chrom"] = raw_matrix.var["chrom"].astype(str)

    #Add column to determine if the feature is mitochondrial or not for QC purposes:
    raw_matrix.var["mt"] = raw_matrix.var["chrom"] == "MT"

    #Add column to determine if the feature is ribosomal or not for QC purposes:
    raw_matrix.var["ribosomal"] = raw_matrix.var["gene_name"].str.match(r"^(RP[SL]\d|MRP[SL]\d)", case=False, na=False)

    #Make CB unique:
    raw_matrix.obs_names_make_unique()

    #Include sample name and condition:
    raw_matrix.obs['sample'] = '{}_{}'.format(condition, replicate)
    raw_matrix.obs['condition'] = condition
    raw_matrix.obs['replicate'] = replicate

    return raw_matrix

def cell_specie_assignment(raw_matrix, specie_threshold, sample_id, output, feature_type):
    
    #Extract species to assign:
    species_to_assign = raw_matrix.var['organism'].values
    X = raw_matrix.X

    if sp.issparse(X):
        X = X.tocsr()

    cell_labels = []
    df_proportions = []

    dominance_threshold = float(specie_threshold)
    
    for i in range(X.shape[0]):
        if sp.issparse(X):
            gene_indices = X[i].indices
            gene_counts = X[i].data
        else:
            gene_indices = np.where(X[i] > 0)[0]
            gene_counts = X[i, gene_indices]
            
        if len(gene_indices) == 0:
            cell_labels.append('Empty')
            continue
        
        # Map gene indices to genomes and sum counts per genome
        genome_counts = pd.Series(gene_counts, index=species_to_assign[gene_indices]).groupby(level=0, observed=True).sum()
        
        # Normalize to get proportions
        total_counts = genome_counts.sum()
        proportions = genome_counts / total_counts

        # Find the dominant genome
        dominant_genome = proportions.idxmax()
        dominance = proportions.max()

        # Save proportion data to generate violin plots:
        df_proportions.append({"species": dominant_genome, "prop": float(dominance)})
        
        # Assign cell to specie if dominant proportion is above the threshold:
        if dominance >= dominance_threshold:
            cell_labels.append(dominant_genome)
        else:
            cell_labels.append('mixed')
        
    #Attach labels to the AnnData object
    raw_matrix.obs['specie_assignment'] = cell_labels

    #Plot distribution of the proportion:
    df_proportions = pd.DataFrame(df_proportions)
    fig_specie, ax = plt.subplots(figsize=(10, 6))

    # Alphabetical order of categories
    order = sorted(df_proportions.iloc[:,0].unique())

    sns.violinplot(
        data = df_proportions,
        x = df_proportions.iloc[:,0],      
        y = df_proportions.iloc[:,1],      
        cut = 0, 
        order = order, 
        ax = ax
    )
    
    ax.set_ylim(0,1)
    ax.set_title('Highest specie proportion per cell - {}'.format(sample_id))
    ax.set_ylabel("Value")
    ax.set_xlabel("")
    
    fig_specie.savefig('{}/{}/unfiltered/{}_specie_proportion_threshold-{}.pdf'.format(output, feature_type, sample_id, specie_threshold), format='pdf', bbox_inches='tight')
    plt.close(fig_specie)

    return raw_matrix

def feature_counts_stats(ann_data):

    #Declare variables to store results:
    number_cells = list()
    number_features = list()
    median_counts_per_cell = list()
    median_features_per_cell = list()

    for individual_sample in ann_data.obs['sample'].unique():
        #Subset based on sample:
        subset_ann_data = ann_data[ann_data.obs['sample'] == individual_sample]

        #Get the number of cells and features:
        number_cells.append(len(subset_ann_data.obs))
        number_features.append(len(subset_ann_data.var))

        #If X is sparse, use efficient operations
        if sp.issparse(subset_ann_data.X):
            # Sum of counts per cell (row-wise sum)
            counts_per_cell = np.array(subset_ann_data.X.sum(axis=1)).flatten()
        
            # Number of expressed features per cell (non-zeros per row)
            features_per_cell = np.array((subset_ann_data.X > 0).sum(axis=1)).flatten()
        else:
            # Dense matrix handling
            counts_per_cell = subset_ann_data.X.sum(axis=1)
            features_per_cell = (subset_ann_data.X > 0).sum(axis=1)
        
        #Compute medians
        median_counts_per_cell.append(np.median(counts_per_cell))
        median_features_per_cell.append(np.median(features_per_cell))

    return number_cells, number_features, median_counts_per_cell, median_features_per_cell

def calculate_species_stats(ann_data):
    counts = ann_data.obs["specie_assignment"].value_counts(dropna=False)
    percents = counts / counts.sum() * 100

    assignment_percentages = pd.DataFrame()

    for k in counts.index:
        assignment_percentages[f"{k}_counts"] = [counts.get(k, 0)]
        assignment_percentages[f"{k}_percentage"] = [round(percents.get(k, 0), 2)]
    
    return assignment_percentages

def specie_assignment_stats(ann_data):
    
    #Generate stats for individual samples:
    assignment_stats = pd.DataFrame()

    for individual_sample in ann_data.obs['sample'].unique():
        subset_ann_data = ann_data[ann_data.obs['sample'] == individual_sample]

        individual_stats = calculate_species_stats(subset_ann_data)
        individual_stats.insert(0, 'sample', individual_sample)
        assignment_stats = pd.concat([assignment_stats, individual_stats])
    
    #Generate stats for the whole condition:
    condition_stats = calculate_species_stats(ann_data)
    condition_stats.insert(0, 'sample', 'NA')
    assignment_stats = pd.concat([assignment_stats, condition_stats])
    
    return assignment_stats

def detect_knee(totals_sorted: np.ndarray) -> int:
    """
    Detect the knee point in a ranked count distribution using the maximum
    curvature (second derivative) method on log-log scale.

    Parameters
    ----------
    totals_sorted : np.ndarray
        Counts per barcode, sorted in descending order.

    Returns
    -------
    int
        1-based rank of the detected knee point.
    """
    # Work in log-log space, as the knee is much clearer there
    x = np.log10(np.arange(1, len(totals_sorted) + 1))
    y = np.log10(totals_sorted + 1)  # +1 to avoid log(0)

    # Smooth y with a rolling window to reduce noise
    window = max(5, len(y) // 100)
    y_smooth = np.convolve(y, np.ones(window) / window, mode="same")

    # Second derivative: high curvature = knee
    dy2 = np.diff(y_smooth, n=2)

    # The knee is where the second derivative is most negative (sharpest drop)
    knee_idx = int(np.argmin(dy2)) + 1  # +1 to correct for diff offset
    knee_rank = knee_idx + 1            # convert to 1-based rank

    return knee_rank

def blaze_cell_calling(input_matrix, sample_name, output, feature_type, quantile=0.95, fraction=0.05):
    
    #This function applies blaze-like cell calling for datasets whose tagged reads are not filtered by the demultiplexing software.

    #Calculate the total number of counts per barcode:
    if sp.issparse(input_matrix.X):
        # Sum of counts per cell (row-wise sum)
        counts_per_cell = np.array(input_matrix.X.sum(axis=1)).flatten()
    
    else:
        # Dense matrix handling
        counts_per_cell = input_matrix.X.sum(axis=1)
    
    input_matrix.obs["total_counts"] = counts_per_cell
    
    #Rank cells (barcodes) by counts with descending order:
    totals_sorted = np.sort(counts_per_cell)[::-1]

    # Auto-detect expected_cells from knee if not provided
    expected_cells = detect_knee(totals_sorted)

    #Get the counts of cell barcode in rank expected_cells * quantile:
    target_rank_1based = int(expected_cells * quantile)
    if target_rank_1based < 1:
        target_rank_1based = 1

    ##Now, transform the rank (1-based) into 0-based:
    target_rank = target_rank_1based -1
    target_rank = min(target_rank, len(totals_sorted) - 1)  # clamp if expected_cells > n_obs

    c = float(totals_sorted[target_rank])
    
    #Calculate the threshold:
    threshold = c * fraction
    
    #Filter:
    filtered_matrix = input_matrix[input_matrix.obs["total_counts"] >= threshold].copy()

    #Save the 'true' cell barcodes into a txt file:
    with open("{}/{}/unfiltered/{}_accepted_cell_barcodes.txt".format(output, feature_type, sample_name), "w") as f:
        for bc in filtered_matrix.obs_names:
            f.write(bc + "\n")

    return filtered_matrix

def annotate_merge_matrix(input_matrixes, feature_type, genomes, specie_threshold, conditions, replicates, platform_type, output):

    #Initialize variables to store stats:
    annotated_stats = []

    #Initialize dicts for grouping based on condition:
    grouped = defaultdict(dict)

    #Process individually:
    for idx, single_matrix in enumerate(input_matrixes):
        if single_matrix.endswith('.h5ad'):
            raw_matrix = sc.read_h5ad(single_matrix)
        elif single_matrix.endswith('.h5'):
            raw_matrix = sc.read_10x_h5(single_matrix)

        # Get condition and replicate for this matrix
        current_condition = conditions[idx]
        current_replicate = replicates[idx]

        sample_id = '{}_{}'.format(current_condition, current_replicate)

        #Preprocess features and obs:
        preprocessed_matrix = features_preprocessing(raw_matrix, genomes, feature_type, current_condition, current_replicate)

        #If sample is Argentag/Parse, apply blaze cell calling to obtain raw matrix:
        if platform_type != '10X' and feature_type == 'gene':
            preprocessed_matrix = blaze_cell_calling(preprocessed_matrix, sample_id, output, feature_type)

        elif platform_type != '10X' and feature_type == 'transcript':
            #Read accepted cells file, generated at gene level:
            with open('{}/gene/unfiltered/{}_accepted_cell_barcodes.txt'.format(output, sample_id)) as f:
                cell_ids = [line.strip() for line in f]
            
            #Filter the transcript matrix:
            preprocessed_matrix = preprocessed_matrix[preprocessed_matrix.obs_names.isin(cell_ids)].copy()

        #If needed, perform cell-specie assignment:
        if len(genomes)>=2:
            annotated_matrix = cell_specie_assignment(preprocessed_matrix, specie_threshold, sample_id, output, feature_type)
        else:
            annotated_matrix = preprocessed_matrix

        #Get statistics for the annotated matrix:
        number_cells, number_features, median_counts_per_cell, median_features_per_cell = feature_counts_stats(annotated_matrix)

        # Convert numpy/list values to scalars for clean CSV output
        def to_scalar(val, dtype=float):
            """Convert numpy array, list, or scalar to a Python scalar."""
            if isinstance(val, (list, np.ndarray)):
                val = val[0] if len(val) > 0 else 0
            return dtype(val) if not isinstance(val, (int, float)) else val

        annotated_stats.append({
            "condition": current_condition,
            "replicate": current_replicate,
            "number_cells": int(to_scalar(number_cells, int)),
            "number_features": int(to_scalar(number_features, int)),
            "median_counts_per_cell": float(to_scalar(median_counts_per_cell, float)),
            "median_feature_per_cell": float(to_scalar(median_features_per_cell, float)),
        })

        #Add matrix to the dictionary for grouping later on:
        grouped[sample_id] = annotated_matrix
    
    #Generate the stats dataframe:
    condition_name_for_output = '-'.join(dict.fromkeys(conditions))  # join all unique conditions with hyphen
    pd.DataFrame(annotated_stats).to_csv("{}/{}/unfiltered/{}_{}_unfiltered_stats.csv".format(output, feature_type, condition_name_for_output, feature_type), index=False)

    #Merge matrixes:
    #Extracting metadata:
    all_var = [x.var for x in grouped.values()]

    #Concatenate them
    all_var = pd.concat(all_var, join="outer")

    #Remove duplicates
    all_var = all_var[~all_var.index.duplicated(keep='first')]

    merged_data = sc.concat(grouped, join="outer", label="sample_id", index_unique="-")
    merged_data.var = all_var.loc[merged_data.var_names]

    return merged_data

def generate_QC_plots(AnnData_object, axs, idx_plot1, cluster_by, scatter):
    sc.pl.violin(AnnData_object, ['n_genes_by_counts'], groupby=cluster_by, ax=axs[idx_plot1], show=False, stripplot=False)
    sc.pl.violin(AnnData_object, ['total_counts'], groupby=cluster_by, ax=axs[idx_plot1+1], show=False, stripplot=False)
    sc.pl.violin(AnnData_object, ['pct_counts_mt'], groupby=cluster_by, ax=axs[idx_plot1+2], show=False, stripplot=False)
    sc.pl.violin(AnnData_object, ['pct_counts_ribosomal'], groupby=cluster_by, ax=axs[idx_plot1+3], show=False, stripplot=False)

    if scatter:
        sc.pl.scatter(AnnData_object, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', ax=axs[idx_plot1+4], show=False)

    # Format x-axis labels for 'sample' grouping (condition_replicate → Condition \n Replicate)
    if cluster_by == 'sample':
        for ax in axs[:4]:  # Skip scatter plot
            labels = [label.get_text() for label in ax.get_xticklabels()]
            new_labels = []
            for label in labels:
                if '_' in label:
                    parts = label.split('_')
                    new_label = f"{parts[0]}\n{parts[1]}"
                    new_labels.append(new_label)
                else:
                    new_labels.append(label)
            ax.set_xticklabels(new_labels)

def edit_axis(axs, idx_plot1, num_obs):
    
    # Metrics corresponding to each subplot
    metric_labels = [
        "Genes per cell",
        "Total counts per cell",
        "Percentage MT per cell",
        "Percentage ribosomal per cell",
        None  # scatter plot (counts vs genes)
    ]

    #Format axis labels:
    for i, ax in enumerate(axs[idx_plot1:idx_plot1+5]):
        ax.set_xlabel('')
        ax.set_ylabel(metric_labels[i])

        if i == 4:
            ax.set_xlabel('Total counts per cell')
            ax.set_ylabel('Genes per cell')

    # Add number of observations to the first violin plot
    axs[idx_plot1].text(
        0.01, 0.98,  # x, y (in axis fraction coords)
        f'n={num_obs}', 
        transform=axs[idx_plot1].transAxes,
        verticalalignment='top',
        horizontalalignment='left',
        fontsize=12
    )

def edit_colorbar_placement(subfig_axs, axs):
    
    # Get the figure
    subfig = subfig_axs.get_figure()

    #Remove scatter title:
    subfig_axs.set_title("")

    #Find colorbar axis:
    colorbar_ax = None

    # Find colorbar axes: those not in the main plot axes list
    candidate_colorbars = [ax for ax in subfig.axes if ax not in axs]

    if not candidate_colorbars:
        return    

    #Set the coordinates:
    coordinates = [0.902, 0.11, 0.01, 0.78]
    colorbar_ax = candidate_colorbars[0]

    #Add label MT for colorbar:
    subfig_axs.text(1.01, 1.08, "MT(%)", transform=axs[4].transAxes, rotation=45, va='center', ha='left', fontsize=12)
        
    #Resize/reposition the colorbar axis
    if colorbar_ax is not None:
        colorbar_ax.set_position(coordinates)

def plotting_QC(matrix, sample_name, feature_type, filtering_status, output, cluster_by='replicate', scatter=True):

    #Set up figure for QC plot:
    if scatter:
        num_subplots = 5
        figwidth = 35
    else:
        num_subplots = 4
        figwidth = 28

    fig_qc, axs_qc = plt.subplots(1, num_subplots, figsize=(figwidth, 6))
    axs_qc = axs_qc.flatten()
    
    #Generate plots:
    generate_QC_plots(matrix, axs_qc, 0, cluster_by, scatter)

    #Add sample name as a title:
    fig_qc.suptitle('{} - {} - {}'.format(sample_name, feature_type, filtering_status), fontsize=18, y=0.925)

    #Modify axis of the plots:
    edit_axis(axs_qc, 0, matrix.shape[0])
    
    #Modify colorbar:
    edit_colorbar_placement(axs_qc[3], axs_qc)

    #Save QC figure:
    if cluster_by != 'replicate':
        fig_qc.savefig('{}/{}/{}/{}_{}_{}_{}_QC.pdf'.format(output, feature_type, filtering_status, sample_name, feature_type, filtering_status, cluster_by), dpi='figure', format='pdf', bbox_inches='tight')
    else:
        fig_qc.savefig('{}/{}/{}/{}_{}_{}_QC.pdf'.format(output, feature_type, filtering_status, sample_name, feature_type, filtering_status), dpi='figure', format='pdf', bbox_inches='tight')

    plt.close(fig_qc)

def is_outlier(adata, metric: str, nmads: int):
  M = adata.obs[metric]
  
  outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
      np.median(M) + nmads * median_abs_deviation(M) < M
  )
  return outlier

def scrublet_per_sample(adata, sample_key: str = "sample", n_top_genes: int = 3000, hvg_flavor: str = "seurat_v3", expected_doublet_rate: float = 0.06, sim_doublet_ratio: float = 2.0, n_prin_comps: int = 30, random_state: int = 0, score_key: str = "doublet_score", pred_key: str = "predicted_doublet"):

    #Create variables to store scrubblet results:
    scores = np.full(adata.n_obs, np.nan, dtype=float)
    preds = np.zeros(adata.n_obs, dtype=bool)

    #Run scrubblet sample-wise:
    for sample, idx in adata.obs.groupby(sample_key, observed=True).indices.items():
        
        print(f"[sc.pp.scrublet] sample={sample!r} n_cells={len(idx)}")

        ad = adata[idx].copy()

        # Based scrubblet results on HVGs per sample:
        tmp = ad.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(
            tmp,
            flavor=hvg_flavor,
            n_top_genes=n_top_genes,
            inplace=True,
        )
        hvg_mask = tmp.var["highly_variable"].values
        if hvg_mask.sum() == 0:
            raise RuntimeError(
                f"No HVGs selected for sample {sample!r}. "
                "Try increasing n_top_genes."
            )
        ad = ad[:, hvg_mask].copy()

        # Filter cells with 0 counts (ie: all counts were in non-HVG):
        #rs = np.asarray(ad.X.sum(axis=1)).ravel() if sp.issparse(ad.X) else ad.X.sum(axis=1)
        #ad = ad[rs > 0].copy()

        # Determine a valid number of PCs for this sample
        n_cells = ad.n_obs
        n_genes = ad.n_vars
        max_pcs = min(n_cells - 1, n_genes - 1)

        n_pcs_use = min(n_prin_comps, max_pcs)

        # Run scrubblet:
        sc.pp.scrublet(
            ad,
            expected_doublet_rate=expected_doublet_rate,
            sim_doublet_ratio=sim_doublet_ratio,
            n_prin_comps=n_pcs_use,
            random_state=random_state,
        )

        #Store results:
        rep_scores = ad.obs["doublet_score"].to_numpy()
        rep_preds = ad.obs["predicted_doublet"].to_numpy(dtype=bool)

        scores[idx] = rep_scores
        preds[idx] = rep_preds

    # Write back to the full object (you can rename if you want)
    adata.obs[score_key] = scores
    adata.obs[pred_key] = preds

    summary = (
        adata.obs
        .groupby("sample")["predicted_doublet"]
        .agg(
            n_doublets="sum",
            n_cells="count"
        )
    )

    summary["pct_doublets"] = 100 * summary["n_doublets"] / summary["n_cells"]
    print(summary)


    return adata

def QC_filtering(AnnData_object, feature_type, output, dynamic = False, nmad = 5, minimal_genes = 1000, genes_by_counts_upper = 5000, total_counts = 20000, pct_mito = 15, min_cells_gene = 15):

    #Extract sample name:
    sample_name = AnnData_object.obs['condition'].unique()[0]
    
    ## Filtering based on QC metrics:
    #Filter genes:
    sc.pp.filter_genes(AnnData_object, min_cells = min_cells_gene)

    # Hardcoded filters:
    if not dynamic:

        #Cell filtering based on QC metrics:
        sc.pp.filter_cells(AnnData_object, min_genes = minimal_genes)
        print(minimal_genes)
        print(pct_mito)
        print(genes_by_counts_upper)
        print(total_counts)
        # One mask, vectorized
        mask = (
            (AnnData_object.obs["n_genes_by_counts"] >= minimal_genes) &
            (AnnData_object.obs["pct_counts_mt"] <= pct_mito) &
            (AnnData_object.obs["n_genes_by_counts"] <= genes_by_counts_upper) &
            (AnnData_object.obs["total_counts"] <= total_counts)
            )

        # Single subset + single copy
        filtered_AnnData = AnnData_object[mask].copy()

    #Dynamic filtering thresholds, based on the number of Median Absolute Deviations (MADs) away from the median:
    else:
        AnnData_object.obs["counts_outlier"] = is_outlier(AnnData_object, "total_counts", nmad)
        AnnData_object.obs["genes_outlier"] = is_outlier(AnnData_object, "n_genes_by_counts", nmad)
        AnnData_object.obs["mito_outlier"] = is_outlier(AnnData_object, "pct_counts_mt", nmad) | (AnnData_object.obs["pct_counts_mt"] > pct_mito)

        removed_mito = sum(AnnData_object.obs["mito_outlier"] == True)

        #Filtering:
        AnnData_object.obs["outlier"] = (
            AnnData_object.obs["counts_outlier"]  
            | AnnData_object.obs["genes_outlier"] 
            | AnnData_object.obs["mito_outlier"]
        )
        # Remove outliers
        filtered_AnnData = AnnData_object[~AnnData_object.obs["outlier"]].copy()

    #Remove doublets with scrubblet:
    if feature_type == 'gene':
        filtered_AnnData = scrublet_per_sample(filtered_AnnData)

    #Set up figure for QC plot for filtered data:
    plotting_QC(filtered_AnnData, sample_name, feature_type, 'filtered', output)
    plotting_QC(filtered_AnnData, sample_name, feature_type, 'filtered', output, 'specie_assignment', False)

    return(filtered_AnnData)

def PCA_clustering(data, data_type, study, normalization_method, feature_type, output, by_condition, number_pcs = 50, number_neighbours = 15, plotting = True):

    #Data normalization:
    if normalization_method == 'log-normalized' and data_type == 'uncorrected':

        # Create layer for raw counts and another for normalized counts:
        data.layers["raw_counts"] = data.X.copy()

        # Log-normalization:
        sc.pp.normalize_total(data, target_sum=10000)
        sc.pp.log1p(data)
    
    else:
        next
 
    #Determine highly variable genes (excluding MT):
    data_without_mt = data[:, ~data.var["mt"]].copy()
    sc.pp.highly_variable_genes(data_without_mt, n_top_genes=int(0.1 * data_without_mt.n_vars), flavor='seurat', batch_key='sample')

    #Transfer HVG mask back to the full object:
    data.var["highly_variable"] = False
    data.var.loc[data_without_mt.var_names, "highly_variable"] = data_without_mt.var["highly_variable"].values

    #Run PCA:
    sc.tl.pca(data, mask_var="highly_variable", n_comps=100)

    if plotting:
        sc.pl.pca_variance_ratio(data, n_pcs=number_pcs, log=False, show=False) 
        fig = plt.gcf()                                 
        output_path = '{}/{}/clustering/{}_{}_{}_VarianceExplained.pdf'.format(output, feature_type, study, data_type, feature_type)
        fig.savefig(output_path, bbox_inches='tight')
        plt.show()                                        
        plt.close(fig)

    #Calculate the number of highly variable genes and extract the top10:
    hvgs = data.var.index[data.var.highly_variable]
    print(f"Number of HVGs: {len(hvgs)}")
    print(f"First 10 HVGs: {hvgs[:10].tolist()}")
    
    #Calculate the KNN graph on a lower-dimensional gene expression representation
    if data_type == 'corrected':
        print('Running analysis for corrected data')
        sc.pp.neighbors(data, use_rep='X_pca_harmony', key_added = 'harmony', n_pcs=number_pcs, n_neighbors=number_neighbours, random_state=0)
        sc.tl.umap(data, neighbors_key = 'harmony', random_state=0)
        
    else:
        sc.pp.neighbors(data, n_pcs=number_pcs, n_neighbors=number_neighbours, random_state=0)
        sc.tl.umap(data, random_state=0)
    
    # Store embedding under a data_type-specific key
    umap_key = f'X_umap_{data_type}'
    data.obsm[umap_key] = data.obsm['X_umap'].copy()
    
    #Establish the different resolutions for leiden:
    resolutions = [0.05, 0.1, 0.3, 0.5, 1.0, 1.5]
    key_list = list()
    
    #Perform Leiden clustering for each resolution and store the results with unique keys
    for res in resolutions:
        key = '{}_leiden_res{}'.format(data_type,res)
        key_list.append(key)

        if data_type == 'corrected':
            sc.tl.leiden(data, resolution=res, key_added=key, neighbors_key='harmony', flavor='leidenalg')
        else:
            sc.tl.leiden(data, resolution=res, key_added=key, flavor='leidenalg')
    
    #Plotting the result:
    if plotting:
        sc.pl.umap(data, color=key_list, wspace=0.2, frameon=False, ncols=len(resolutions), show=False) 
        fig = plt.gcf()                                 
        output_path = '{}/{}/clustering/umap_{}_{}_{}_Leiden.pdf'.format(output, feature_type, study, data_type, feature_type)
        fig.savefig(output_path, bbox_inches='tight')
        plt.show()                                        
        plt.close(fig)
    
    #To determine which resolution is better, we use silhouette score:
    
    #Loop through each Leiden resolution that we have computed and calculate the silhouette score:
    silhouette_scores = list()
    for res in key_list:
        n_clusters = data.obs[res].nunique()
        if n_clusters < 2:
            print(f"Silhouette Score for resolution {res}: skipped (only {n_clusters} cluster found)")
            silhouette_scores.append(-1)
            continue
        score = silhouette_score(data.obsm[umap_key], data.obs[res])
        silhouette_scores.append(score)
        print(f"Silhouette Score for resolution {res}: {score:.4f}")
    
    #Select the resolution with the highest silhouette score:
    selected_resolution = key_list[silhouette_scores.index(max(silhouette_scores))]
    
    #Plot the results:
    if plotting:

        if by_condition:
            sc.pl.umap(data, color=['sample', 'condition', 'specie_assignment', selected_resolution], wspace=0.2, frameon=False, ncols=4, show=False) 
        else:
            sc.pl.umap(data, color=['sample', 'specie_assignment', selected_resolution], wspace=0.2, frameon=False, ncols=3, show=False)   
        fig = plt.gcf()                                 
        output_path = '{}/{}/clustering/umap_{}_{}_{}_Leiden_BySample.pdf'.format(output, feature_type, study, data_type, feature_type)
        fig.savefig(output_path, bbox_inches='tight')
        plt.show()                                        
        plt.close(fig)
        
    return data

def harmonpy_correction(data, by_condition):
    
    #Run harmony via scanpy external API - harmony will overwrite the X_umap slot.
    if by_condition:
        sc.external.pp.harmony_integrate(data, key='condition', basis = 'X_pca')
    else:
        sc.external.pp.harmony_integrate(data, key='sample', basis = 'X_pca')

    return data

def _rank_genes_to_df(data: sc.AnnData, key: str, n_genes: int, cluster_a, cluster_b) -> pd.DataFrame:
    """
    Extract rank_genes_groups results into a tidy DataFrame.
 
    Parameters
    ----------
    data    : AnnData
    key     : str   — key in data.uns containing rank_genes_groups results
    n_genes : int   — number of top genes to extract per group
 
    Returns
    -------
    pd.DataFrame with columns:
        cluster, rank, gene, score, log2fc, pval, pval_adj, pts
    """
    result = data.uns[key]
    groups = result["names"].dtype.names
 
    # pts can be a dict-of-arrays (older scanpy) or a DataFrame (newer scanpy)
    pts_data = result.get("pts", None)
    pts_rest_data = result.get("pts_rest", None)

    rows = []
    for group in groups:
        n = min(n_genes, len(result["names"][group]))
        for i in range(n):
            gene = result["names"][group][i]
 
            # Resolve pts (expression in group A) safely regardless of storage format
            pts_val = None
            if pts_data is not None:
                if isinstance(pts_data, pd.DataFrame):
                    pts_val = pts_data.loc[gene, group] if (
                        gene in pts_data.index and group in pts_data.columns
                    ) else None
                elif isinstance(pts_data, dict) and group in pts_data:
                    arr = pts_data[group]
                    pts_val = arr[i] if i < len(arr) else None
            
            # Resolve pts_rest (expression in group B - reference) safely regardless of storage format
            pts_rest_val = None
            if pts_rest_data is not None:
                if isinstance(pts_rest_data, pd.DataFrame):
                    pts_rest_val = pts_rest_data.loc[gene, group] if (
                        gene in pts_rest_data.index and group in pts_rest_data.columns
                    ) else None
                elif isinstance(pts_rest_data, dict) and group in pts_rest_data:
                    arr = pts_rest_data[group]
                    pts_rest_val = arr[i] if i < len(arr) else None
 
            rows.append({
                "cluster":  group,
                "rank":     i + 1,
                "gene":     gene,
                "score":    result["scores"][group][i],
                "log2fc":   result["logfoldchanges"][group][i],
                "pval":     result["pvals"][group][i],
                "pval_adj": result["pvals_adj"][group][i],
                f"pts_{cluster_a}":      pts_val,
                f"pts_{cluster_b}": pts_rest_val,
            })
 
    return pd.DataFrame(rows)

def Global_DE_PerCluster(
    data,
    condition,
    feature_type,
    selected_resolution,
    output,
    number_top_DE_genes,
):
    """
    One-vs-rest Wilcoxon marker detection for every cluster in one platform.
    Plots a dotplot with var_names (display_name) on the x-axis.
    """

    out_dir = Path(output) / feature_type / "marker_genes"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Generate display name ────────────────────────────────────────────────
    if feature_type == 'transcript':
        has_name = data.var['transcript_name'].notna()
        display_name = (
            data.var['organism'].astype(str) + '---' +
            data.var['transcript_id'].astype(str) + ' (' +
            data.var['transcript_name'].astype(str) + ')'
        ).where(
            has_name,
            data.var['organism'].astype(str) + '---' + data.var['features_id'].astype(str)
        )
    else:
        has_name = data.var['gene_name'].notna()
        display_name = (
            data.var['organism'].astype(str) + '---' +
            data.var['features_id'].astype(str) + ' (' +
            data.var['gene_name'].astype(str) + ')'
        ).where(
            has_name,
            data.var['organism'].astype(str) + '---' + data.var['features_id'].astype(str)
        )

    data.var['display_name'] = display_name.values
    data.var.index = display_name.values
    data.raw = data

    # ── DE test ──────────────────────────────────────────────────────────────
    sc.tl.rank_genes_groups(
        data,
        groupby   = selected_resolution,
        method    = "wilcoxon",
        key_added = "rank_genes_groups_wilcoxon",
        pts       = True,
    )

    # ── Extract results into a DataFrame using var_names directly ────────────
    result  = data.uns["rank_genes_groups_wilcoxon"]
    groups  = result["names"].dtype.names
    pts_data = result.get("pts", None)
    
    rows = []
    for group in groups:

        n = min(number_top_DE_genes, len(result["names"][group]))
        for i in range(n):
            raw_id = result["names"][group][i]   # whatever scanpy stored

            # ── Resolve to var_name ──────────────────────────────────────────
            # Try 1: already a valid var_name
            if raw_id in data.var_names:
                gene_name = raw_id

            # Try 2: it's a features_id
            elif 'features_id' in data.var.columns and raw_id in data.var['features_id'].astype(str).values:
                mask      = data.var['features_id'].astype(str) == raw_id
                gene_name = data.var_names[mask][0]

            # Try 3: it's a positional integer
            else:
                try:
                    gene_name = data.var_names[int(raw_id)]
                except (ValueError, IndexError):
                    logger.warning(f"      Could not resolve '{raw_id}' — keeping as-is")
                    gene_name = raw_id

            # ── pts ──────────────────────────────────────────────────────────
            pts_val = None
            if pts_data is not None:
                if isinstance(pts_data, pd.DataFrame):
                    pts_val = pts_data.loc[gene_name, group] if (
                        gene_name in pts_data.index and group in pts_data.columns
                    ) else None
                elif isinstance(pts_data, dict) and group in pts_data:
                    arr     = pts_data[group]
                    pts_val = arr[i] if i < len(arr) else None

            rows.append({
                "cluster":  group,
                "rank":     i + 1,
                "feature":     gene_name,
                "score":    result["scores"][group][i],
                "log2fc":   result["logfoldchanges"][group][i],
                "pval":     result["pvals"][group][i],
                "pval_adj": result["pvals_adj"][group][i],
                "pts":      pts_val,
            })

    marker_df = pd.DataFrame(rows)
    
    # ── Save table ───────────────────────────────────────────────────────────
    table_path = out_dir / f"{condition}_{feature_type}_{selected_resolution}_TopDE_table.tsv"
    marker_df.to_csv(table_path, sep="\t", index=False)
    logger.info(f"    Saved table  → {table_path}")

    # ── Dotplot ──────────────────────────────────────────────────────────────
    n_clusters = len(marker_df.iloc[0].unique())
    fig_width  = max(10, number_top_DE_genes * n_clusters * 0.45)
    fig_height = max(4,  n_clusters * 0.65)
    
    dp = sc.pl.rank_genes_groups_dotplot(
        data,
        key            = 'rank_genes_groups_wilcoxon',
        groupby        = selected_resolution,
        standard_scale = "var",
        figsize        = (fig_width, fig_height),
        show = False
    )

    fig      = plt.gcf()
    plot_path = out_dir / f"{condition}_{feature_type}_{selected_resolution}_TopDE.pdf"
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"    Saved dotplot → {plot_path}")

    return data

def decoupler_pseudobulk(data, grouping_variable, sample_column, min_cells):
    
    # Ensure that the grouping variable is filled:
    data.obs[grouping_variable] = data.obs[grouping_variable].cat.add_categories(['Unassigned'])
    data.obs[grouping_variable] = data.obs[grouping_variable].fillna('Unassigned')

    # Run pseudobulk:
    pdata = dc.get_pseudobulk(
        data,
        sample_col=sample_column,
        groups_col=grouping_variable,
        mode="sum",
        min_cells=min_cells,
        skip_checks=True
    )

    return pdata

def filter_by_expression(adata, group):
    """
    Filter lowly expressed genes using edgeR filterByExpr function.

    :adata: anndata object
    :group: column name in the anndata obs layer containing sample grouping
    
    :return: anndata object with lowly expressed genes removed
    """
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    import rpy2.robjects.numpy2ri
    
    rpy2.robjects.numpy2ri.activate()
    edger = importr("edgeR")

    # Ensure adata.X is a NumPy array
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
    
    keep = edger.filterByExpr(X.T, ro.FactorVector(adata.obs[group]))

    return adata[:, list(keep)]

def pseudo_bulk_dge_engine(
    data,
    label,                  # column name in the anndata obs layer containing cell grouping 
    sample_col="sample",    # column with replicate information
    condition_col=None,     # e.g. "Condition"
    baseline=None,          # reference level, e.g. "IgG" or "ControlCluster"
    comparisons=None,       # explicit list of (A, B) tuples, e.g. [("T_cells", "B_cells")]
    pairwise=False,         # auto-generate all pairwise comparisons
):
    """
    Production-grade pseudobulk DE engine using edgeR via pertpy.

    Comparison resolution (in order of priority):
        1. `comparisons` — explicit list of (group_to_compare, baseline) tuples
        2. `baseline`    — auto-generate (g, baseline) for all other groups
        3. `pairwise=True` — auto-generate all pairwise combinations
    """

    @dataclass
    class dgeResults:
        edger_object: pt.tools.EdgeR
        de_result: pd.core.frame.DataFrame

        def plot_volcano(self, **kwargs):
            self.edger_object.plot_volcano(self.de_result, **kwargs)
 
    # ------------------------------------------------------------------
    # RESOLVE COMPARISONS (explicit > baseline > pairwise)
    # ------------------------------------------------------------------
    if comparisons is not None:
        # Validate that all specified groups exist in the data
        known_groups = data.obs[condition_col].dropna().unique().tolist()
        for a, b in comparisons:
            for g in (a, b):
                if g not in known_groups:
                    raise ValueError(
                        f"Group '{g}' in comparisons not found in "
                        f"adata.obs['{condition_col}']. Available: {known_groups}"
                    )
        logger.info(f"Using {len(comparisons)} explicitly provided comparisons")
    else:
        groups = data.obs[condition_col].dropna().unique().tolist()
 
        if pairwise:
            comparisons = list(combinations(groups, 2))
        else:
            if baseline not in groups:
                raise ValueError(
                    f"baseline='{baseline}' not found in adata.obs['{condition_col}']. "
                    f"Available: {groups}"
                )
            comparisons = [(g, baseline) for g in groups if g != baseline]
 
    logger.info(f"Running {len(comparisons)} comparisons: {comparisons}")

    # Go through each comparison:
    results = {}

    for comparison in comparisons:
        try:
            design = '~{} + {}'.format(sample_col, condition_col)
            
            # Determine grouping labels
            if label is not None:
                grouping_labels = data.obs[label].unique()
            else:
                grouping_labels = ['']
            
            # Process each group
            for l in grouping_labels:
                try:
                    # If cell grouping:
                    if label is not None:
                        subset = data[data.obs[label] == l]
                    else:
                        subset = data.copy()
                    
                    # Filter based on expression:
                    subset = filter_by_expression(subset, condition_col)
                    
                    # Fit model:
                    edgr = pt.tl.EdgeR(subset, design=design)
                    edgr.fit(robust=True)
                    result = edgr.test_contrasts(edgr.contrast(column=condition_col, baseline=comparison[1], group_to_compare=comparison[0]))
                    result.index = result["variable"]
                    
                    # Store results with composite key (comparison + label)
                    key = f"{comparison[0]}_vs_{comparison[1]}_{l}"
                    results[key] = dgeResults(edgr, result)
                    
                except Exception as e:
                    print(f"Invalid comparison for {l}: {e}")
                    
        except Exception as e:
            print(f"Error processing comparison {comparison}: {e}")

    # Convert results to dataframe
    results_list = []
    for key, dge_result in results.items():
        comparison_parts = key.rsplit('_', 1)  # Split from right to separate label
        comparison_name = comparison_parts[0]
        group_label = comparison_parts[1] if len(comparison_parts) > 1 else 'No_cell_grouping'
        
        # Assuming dge_result has attributes or can be converted to dict
        result_df = dge_result.de_result.copy()  # or however your dgeResults stores the dataframe
        result_df['comparison'] = comparison_name
        result_df['cell_group'] = group_label
        
        results_list.append(result_df)

    # Combine all results into single dataframe
    results_df = pd.concat(results_list, ignore_index=True)
    print(results_df)
    return results_df

def detect_count_type(adata, tol=1e-6):
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    
    is_int = np.abs(X - np.round(X)) < tol
    gene_int = is_int.all(axis=0)

    return {
        "n_genes": X.shape[1],
        "n_integer_genes": gene_int.sum(),
        "n_fractional_genes": (~gene_int).sum(),
        "fraction_integer_genes": gene_int.mean()
    }

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def plot_volcano(
    result_df,
    gene_col="feature_name",       # column with gene labels
    logfc_col="log_fc",
    pval_col="adj_p_value",
    logcpm_col="logCPM",
    log_fc_thresh=2,
    fdr_thresh=0.01,
    logcpm_thresh=1.0,
    top_n_labels=15,            # label the top N most significant genes
    title=None,
    figsize=(8, 6),
):
    df = result_df.copy()

    # -log10 transform p-values, clip to avoid inf
    df["_neglog10p"] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Assign category for coloring
    def categorize(row):
        if row[logcpm_col] <= logcpm_thresh:
            return "low_expr"
        if row[pval_col] < fdr_thresh and row[logfc_col] > log_fc_thresh:
            return "up"
        if row[pval_col] < fdr_thresh and row[logfc_col] < -log_fc_thresh:
            return "down"
        return "ns"

    df["_cat"] = df.apply(categorize, axis=1)

    colors = {
        "up":       "#d62728",   # red
        "down":     "#1f77b4",   # blue
        "ns":       "#aaaaaa",   # grey
        "low_expr": "#dddddd",   # light grey
    }

    fig, ax = plt.subplots(figsize=figsize)

    for cat, grp in df.groupby("_cat"):
        ax.scatter(
            grp[logfc_col],
            grp["_neglog10p"],
            c=colors[cat],
            s=13,
            alpha=0.6,
            linewidths=0,
            label=cat,
            rasterized=True,    # important for large datasets
        )

    # Threshold lines
    ax.axhline(-np.log10(fdr_thresh), color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axvline( log_fc_thresh,        color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axvline(-log_fc_thresh,        color="black", lw=0.8, ls="--", alpha=0.5)

    # Label top N significant genes
    sig = df[df["_cat"].isin(["up", "down"])].nsmallest(top_n_labels, pval_col)

    from adjustText import adjust_text

    texts = []
    for _, row in sig.iterrows():
        texts.append(
            ax.text(
                row[logfc_col],
                row["_neglog10p"],
                row[gene_col],
                fontsize=6,
                ha="left",
                va="bottom",
            )
        )

    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5)
)

    # Counts in legend
    n_up   = (df["_cat"] == "up").sum()
    n_down = (df["_cat"] == "down").sum()

    legend_patches = [
        mpatches.Patch(color=colors["up"],       label=f"Up ({n_up})"),
        mpatches.Patch(color=colors["down"],      label=f"Down ({n_down})"),
        mpatches.Patch(color=colors["ns"],        label="NS"),
        mpatches.Patch(color=colors["low_expr"],  label="Low expr"),
    ]
    ax.legend(handles=legend_patches, fontsize=8, frameon=False)

    ax.set_xlabel("log2 Fold Change", fontsize=11)
    ax.set_ylabel("-log10 FDR", fontsize=11)
    ax.set_title(title or "Volcano plot", fontsize=12)
    sns.despine(ax=ax)

    plt.tight_layout()
    return fig

def plot_composite_umap_panel(adata, features, feature_type, figsize = None, cmap = "Reds", ncols = 3, vmin = None, vmax = None):
    
    # Validate UMAP exists
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates not found. Run sc.tl.umap(adata) first.")

    # Extract var_name from features:
    var_name = list()
    filtered_features = list()

    for feature in features:
        if feature_type == 'gene':
            mask = adata.var['gene_name'] == feature
        
        else:
            mask = adata.var['transcript_name'] == feature

        if not mask.any():
            print(f"Feature '{feature}' not found in var.")
            continue
        
        # Valid features:
        filtered_features.append(feature)

        # Get the corresponding var_name
        var_name.append(adata.var_names[mask][0])
            
    # Validate features exist:
    features = filtered_features
    all_valid_features = list(adata.var_names) + list(adata.obs.columns)
    invalid_features = [f for f in var_name if f not in all_valid_features]
    if invalid_features:
        raise ValueError(f"Features not found: {invalid_features}")
    
    # Auto-calculate figsize if not provided
    n_features = len(features)
    nrows = int(np.ceil(n_features / ncols))
    if figsize is None:
        figsize = (ncols * 4, nrows * 4)
    
    # Create subplot grid
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        constrained_layout=True
    )
    
    # Handle case of single subplot
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    else:
        axes = axes.flatten()
        
    # Plot each feature using sc.pl.umap
    for idx, (ind_var_name, feature) in enumerate(zip(var_name, features)):
        ax = axes[idx]
        
        # Use scanpy's umap plotting
        sc.pl.umap(
            adata,
            color=ind_var_name,
            ax=ax,
            show=False,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            title = feature, 
            use_raw=False,
        )
    
    # Hide unused subplots
    for idx in range(n_features, len(axes)):
        axes[idx].set_visible(False)
        
    return fig

# ---------------------------------------------------------------------------
# Step functions
# ---------------------------------------------------------------------------

def run_annotate(args: argparse.Namespace) -> None:
    """Annotate features from input matrices."""
    input_matrixes = args.input
    logger.info("=== ANNOTATE FEATURES ===")
    logger.info("  Input matrices : %s", [str(p) for p in input_matrixes])
    logger.info("  Feature type   : %s", args.feature_type)
    logger.info("  Platform       : %s", args.platform_type)

    # Create output directory:
    Path('{}/{}/unfiltered'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)
    Path('{}/{}/matrixes'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    # Use genomes from arguments:
    genomes = args.genomes

    # Build condition list (handle single value for all matrices or N values)
    n_matrices = len(input_matrixes)
    if args.condition:
        if len(args.condition) == 1:
            conditions = args.condition * n_matrices
            logger.info("Using single condition '%s' for all %d matrices", args.condition[0], n_matrices)
        elif len(args.condition) == n_matrices:
            conditions = args.condition
            logger.info("Using provided conditions for %d matrices", n_matrices)
        else:
            raise ValueError(f"Number of conditions ({len(args.condition)}) must be 1 or match number of matrices ({n_matrices})")
    else:
        conditions = []
        for matrix in input_matrixes:
            basename = os.path.splitext(os.path.basename(matrix))[0].split('_')
            conditions.append(basename[0])
        logger.info("Extracted conditions from filenames: %s", conditions)

    # Build replicate list (handle single value for all matrices or N values)
    if args.replicate:
        if len(args.replicate) == 1:
            replicates = args.replicate * n_matrices
            logger.info("Using single replicate '%s' for all %d matrices", args.replicate[0], n_matrices)
        elif len(args.replicate) == n_matrices:
            replicates = args.replicate
            logger.info("Using provided replicates for %d matrices", n_matrices)
        else:
            raise ValueError(f"Number of replicates ({len(args.replicate)}) must be 1 or match number of matrices ({n_matrices})")
    else:
        replicates = []
        for matrix in input_matrixes:
            basename = os.path.splitext(os.path.basename(matrix))[0].split('_')
            replicates.append(basename[1] if len(basename) > 1 else "1")
        logger.info("Extracted replicates from filenames: %s", replicates)

    annotated_matrix = annotate_merge_matrix(input_matrixes, args.feature_type, genomes, args.specie_threshold, conditions, replicates, args.platform_type, args.output)

    # Create condition name for output (join all unique conditions with hyphen)
    condition_name = '-'.join(dict.fromkeys(conditions))  # dict.fromkeys preserves order and removes duplicates
    logger.info("Output condition name: %s", condition_name)

    #Compute cell-specie assignment stats (if needed):
    if len(genomes)>=2:
        cell_specie_stats = specie_assignment_stats(annotated_matrix)
        cell_specie_stats.insert(0, 'condition', condition_name)
        cell_specie_stats.to_csv("{}/{}/unfiltered/{}_{}_unfiltered_cell_specie_stats.csv".format(args.output, args.feature_type, condition_name, args.feature_type), index=False)

    #Calculate QC metrics:
    sc.pp.calculate_qc_metrics(annotated_matrix, qc_vars=["mt", "ribosomal"], inplace=True, percent_top=[20], log1p=False)

    #Save merged matrix:
    annotated_matrix.write('{}/{}/matrixes/{}_{}_annotated_matrix.h5ad'.format(args.output, args.feature_type, condition_name, args.feature_type))

    # Use combined condition name for display
    sample_name = condition_name

    # Determine grouping: use 'sample' if multiple conditions, 'replicate' if single condition
    groupby_col = 'sample' if len(set(conditions)) > 1 else 'replicate'

    #Generate QC plots pre filtering by sample/replicate and by specie_assignment:
    plotting_QC(annotated_matrix, sample_name, args.feature_type, 'unfiltered', args.output, cluster_by=groupby_col)
    if len(genomes) >= 2 and 'specie_assignment' in annotated_matrix.obs.columns:
        plotting_QC(annotated_matrix, sample_name, args.feature_type, 'unfiltered', args.output, cluster_by='specie_assignment', scatter=False)

def run_filter(args: argparse.Namespace) -> None:
    """Filter cells and features."""
    logger.info("=== FILTERING ===")

    # Create output directory:
    Path('{}/{}/filtered'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    #Define genomes:
    genomes = ['hsapiens_gene_ensembl', 'mmusculus_gene_ensembl', 'clfamiliaris_gene_ensembl']

    #Extract condition:
    basename = os.path.splitext(os.path.basename(args.input[0]))[0].split('_')
    condition = basename[0]

    #Filtering:
    unfiltered_matrix = sc.read_h5ad(args.input[0])
    filtered_annotated_matrix = QC_filtering(unfiltered_matrix, args.feature_type, args.output)

    #Compute cell-specie assignment stats (if needed):
    if len(genomes)>=2:
        cell_specie_stats = specie_assignment_stats(filtered_annotated_matrix)
        cell_specie_stats.insert(0, 'condition', condition)
        cell_specie_stats.to_csv("{}/{}/filtered/{}_{}_filtered_cell_specie_stats.csv".format(args.output, args.feature_type, condition, args.feature_type), index=False)
    
    #Get statistics for the annotated matrix:
    number_cells, number_features, median_counts_per_cell, median_features_per_cell = feature_counts_stats(filtered_annotated_matrix)
    filtered_annotated_stats = list()
    filtered_annotated_stats.append({
        "condition": condition,
        "number_cells": number_cells,
        "number_features": number_features,
        "median_counts_per_cell": median_counts_per_cell,
        "median_feature_per_cell": median_features_per_cell,
    })

    #Generate the stats dataframe: 
    pd.DataFrame(filtered_annotated_stats).to_csv("{}/{}/filtered/{}_{}_filtered_stats.csv".format(args.output, args.feature_type, condition, args.feature_type), index=False)

    #Save the filtered matrix:
    filtered_annotated_matrix.write('{}/{}/matrixes/{}_{}_filtered_annotated_matrix.h5ad'.format(args.output, args.feature_type, condition, args.feature_type))

def run_cluster(args: argparse.Namespace) -> None:
    """Cluster cells (e.g. Leiden / Louvain)."""
    logger.info("=== CLUSTERING ===")

    # Create output directory:
    Path('{}/{}/clustering'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    #Extract condition:
    basename = os.path.splitext(os.path.basename(args.input[0]))[0].split('_')
    condition = basename[0]

    #Define normalization method:
    normalization_method = 'log-normalized'

    #Read input:
    data = sc.read_h5ad(args.input[0])

    #Filtering mixed cells:
    if 'specie_assignment' in data.obs.columns:
        annotated_matrix = data[data.obs['specie_assignment'] != 'mixed'].copy()
        print(f"Remaining cells after removing 'mixed': {data.n_obs}")
    else:
        annotated_matrix = data
    
    #Initial clustering:
    clustered_uncorrected_matrix = PCA_clustering(annotated_matrix, 'uncorrected', condition, normalization_method, args.feature_type, args.output, args.by_condition, number_pcs = 50, number_neighbours = 15, plotting = True)

    #Batch correction with harmonpy:
    corrected_matrix = harmonpy_correction(clustered_uncorrected_matrix, args.by_condition)

    #Re-cluster again: 
    clustered_corrected_matrix = PCA_clustering(corrected_matrix, 'corrected', condition, normalization_method, args.feature_type, args.output, args.by_condition, number_pcs = 50, number_neighbours = 15, plotting = True)

    #Save matrixes:
    clustered_uncorrected_matrix.write('{}/{}/matrixes/{}_{}_uncorrected_clustered_matrix.h5ad'.format(args.output, args.feature_type, condition, args.feature_type))
    clustered_corrected_matrix.write('{}/{}/matrixes/{}_{}_corrected_clustered_matrix.h5ad'.format(args.output, args.feature_type, condition, args.feature_type))

def run_marker_genes(args: argparse.Namespace) -> None:
    """Indentify marker genes per cluster analysis."""
    logger.info("=== MARKER GENES PER CLUSTER ===")

    # Create output directory:
    Path('{}/{}/diff_expression'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    #Extract condition:
    basename = os.path.splitext(os.path.basename(args.input[0]))[0].split('_')
    condition = basename[0]

    #Read input:
    data = sc.read_h5ad(args.input[0])
    
    # ── Guard: confirm raw counts layer exists ────────────────────────────────
    counts_layer = "raw_counts"
    if counts_layer not in data.layers:
        raise KeyError(
            f"Expected a '{counts_layer}' layer with raw integer counts. "
            f"Available layers: {list(data.layers.keys())}"
        )
 
    # ── Resolve resolutions ───────────────────────────────────────────────────
    # args.clustering_resolution can be a single string or a list
    if isinstance(args.clustering_resolution, str):
        resolutions = [args.clustering_resolution]
    else:
        resolutions = list(args.clustering_resolution)
 
    # ── Main loop ──────────────────────────────────────
    logger.info(f"  Cells: {data.n_obs}")

    for resolution in resolutions:
        logger.info(f"  Resolution: {resolution}")

        if resolution not in data.obs.columns:
            logger.warning(
                f"  '{resolution}' not found in obs "
            )
            continue

        clusters = data.obs[resolution].unique().tolist()
        cells_per_cluster = data.obs[resolution].value_counts()

        for cluster, n_cells in cells_per_cluster.items():
            logger.info(f"  Cluster {cluster}: {n_cells} cells")

        # Step 1: one-vs-rest
        logger.info("Comparing one cluster vs the rest to identify marker genes...")
        Global_DE_PerCluster(
            data                = data,
            condition           = condition,
            feature_type        = args.feature_type,
            selected_resolution = resolution,
            output              = args.output,
            number_top_DE_genes = 10,
        )

def run_merge(args: argparse.Namespace) -> None:
    """Merging N input h5ad matrixes into one."""
    logger.info("=== MERGING MATRIXES (*.h5ad) ===")

    # Create output directory:
    Path('{}/{}/matrixes'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    # Set up the dictonary:
    grouped = dict()

    # Read all matrixes:
    for input_matrix in args.input:
        single_matrix = sc.read_h5ad(input_matrix)
        
        condition = single_matrix.obs['condition'].iloc[0]

        #Add matrix to the dictionary for grouping later on:
        grouped[condition] = single_matrix
    
    #Merge matrixes:
    #Extracting metadata:
    all_var = [x.var for x in grouped.values()]

    #Concatenate them
    all_var = pd.concat(all_var, join="outer")

    #Remove duplicates
    all_var = all_var[~all_var.index.duplicated(keep='first')]

    merged_data = sc.concat(grouped, join="outer", label="sample_id", index_unique="-")
    merged_data.var = all_var.loc[merged_data.var_names]

    #Save matrixes:
    merged_data.write('{}/{}/matrixes/{}_matrix.h5ad'.format(args.output, args.feature_type, args.name))

def run_diffexp(args: argparse.Namespace) -> None:
    """Pseudobulk aggregation and differential expression (EdgeR)."""
    logger.info("=== PSEUDOBULK AND DIFFERENTIAL EXPRESSION ===")

    # Create output directory:
    Path('{}/{}/diff_expression'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    # Extract condition:
    basename = os.path.splitext(os.path.basename(args.input[0]))[0].split('_')
    condition = basename[0]

    # Read input:
    data = sc.read_h5ad(args.input[0])

    # ── Confirm raw counts layer exists ────────────────────────────────
    counts_layer = "raw_counts"
    if counts_layer not in data.layers:
        raise KeyError(
            f"Expected a '{counts_layer}' layer with raw integer counts. "
            f"Available layers: {list(data.layers.keys())}"
        )
    
    else:
        data.X = data.layers['raw_counts']

    # ── Generate display name ────────────────────────────────────────────────
    if args.feature_type == 'transcript':
        has_name = data.var['transcript_name'].notna()
        display_name = (
            data.var['organism'].astype(str) + '---' +
            data.var['transcript_id'].astype(str) + ' (' +
            data.var['transcript_name'].astype(str) + ')'
        ).where(
            has_name,
            data.var['organism'].astype(str) + '---' + data.var['features_id'].astype(str)
        )
    else:
        has_name = data.var['gene_name'].notna()
        display_name = (
            data.var['organism'].astype(str) + '---' +
            data.var['features_id'].astype(str) + ' (' +
            data.var['gene_name'].astype(str) + ')'
        ).where(
            has_name,
            data.var['organism'].astype(str) + '---' + data.var['features_id'].astype(str)
        )

    data.var['display_name'] = display_name.values
    data.var.index = display_name.values
    data.raw = data

    # Pseudobulk via decoupler aggregation:
    pseudobulk_data = decoupler_pseudobulk(data, args.grouping, args.sample_column, args.min_cells)
    
    ## NOTE: not good! Adapting due to isoquant quantification:
    pseudobulk_data.X = np.round(pseudobulk_data.X)
    
    # Run DE:
    conditions = data.obs['condition'].unique()

    # Define design:
    if len(conditions) == 1:

        if args.comparisons is None:
            result = pseudo_bulk_dge_engine(
                pseudobulk_data,
                label=args.grouping,
                sample_col=args.sample_column,
                pairwise=True,
            )

        else:
            result = pseudo_bulk_dge_engine(
                pseudobulk_data,
                None,
                args.sample_column,
                args.grouping,
                pairwise=False,
                comparisons=args.comparisons
            )

    else:

        if args.comparisons is None:
            result = pseudo_bulk_dge_engine(
                pseudobulk_data,
                label=args.grouping,
                sample_col=args.sample_column,
                pairwise=True,
                condition_col="condition",
            )

        else:
            result = pseudo_bulk_dge_engine(
                pseudobulk_data,
                label=args.grouping,
                sample_col=args.sample_column,
                pairwise=False,
                comparisons=args.comparisons,
                condition_col="condition"
            )

    # Filter results using provided thresholds:
    LOG_FC_THRESHOLD = args.log_fc_threshold
    FDR_THRESHOLD    = args.fdr_threshold
    LOGCPM_THRESHOLD = args.logcpm_threshold

    filtered = result[
        (result["adj_p_value"] < FDR_THRESHOLD) &
        (result["log_fc"].abs() > LOG_FC_THRESHOLD) &
        (result["logCPM"] > LOGCPM_THRESHOLD)
    ].copy()

    # Export results into a table (tsv):
    result.to_csv('{}/{}/diff_expression/DE_results.tsv'.format(args.output, args.feature_type), sep='\t')
    filtered.to_csv('{}/{}/diff_expression/DE_filtered_results.tsv'.format(args.output, args.feature_type), sep='\t')

    # Generate proper labels for visualization:
    inside_parens = result['variable'].str.extract(r'\((.*?)\)', expand=False)

    # split on one or more hyphens and take second part
    after_hyphen = result['variable'].str.split(r'-+', n=1).str[1]

    # combine: parentheses first, otherwise hyphen result
    result['feature_name'] = inside_parens.fillna(after_hyphen).str.strip()

    # Generate volcano plots:
    if len(conditions) == 1:
        for comparison, grp in result.groupby("comparison"):
            fig = plot_volcano(grp, title=comparison, log_fc_thresh = LOG_FC_THRESHOLD, fdr_thresh=FDR_THRESHOLD, logcpm_thresh=LOGCPM_THRESHOLD)
            fig.savefig(f"{args.output}/{args.feature_type}/diff_expression/volcano_{comparison}.pdf", dpi=150)
            plt.close()
    else:
        for (cluster, comparison), grp in result.groupby(["cluster", "comparison"]):
            fig = plot_volcano(grp, title=f"{cluster} | {comparison}", log_fc_thresh = LOG_FC_THRESHOLD, fdr_thresh=FDR_THRESHOLD, logcpm_thresh=LOGCPM_THRESHOLD)
            fig.savefig(f"{args.output}/{args.feature_type}/diff_expression/volcano_{cluster}_{comparison}.pdf", dpi=150)
            plt.close()

def run_umap_plotting(args: argparse.Namespace) -> None:
    """Plotting features' expression in UMAP."""
    logger.info("=== PLOTTING EXPRESSION UMAP ===")

    # Create output directory:
    Path('{}/{}/plotting'.format(args.output, args.feature_type)).mkdir(parents=True, exist_ok=True)

    # Extract condition:
    basename = os.path.splitext(os.path.basename(args.input[0]))[0].split('_')
    condition = basename[0]

    # Read input:
    data = sc.read_h5ad(args.input[0])
    
    # Plotting:
    fig = plot_composite_umap_panel(
        data,
        args.features,
        args.feature_type,
        ncols=3,
        cmap='Reds'
    )

    # Save manually
    fig.savefig('{}/{}/plotting/UMAP_{}_{}.pdf'.format(args.output, args.feature_type, condition, '_'.join(args.features)), dpi=300, bbox_inches='tight')

# ---------------------------------------------------------------------------
# Full pipeline (default when no subcommand is given)
# ---------------------------------------------------------------------------

def run_all(args: argparse.Namespace) -> None:
    """Run every step in the canonical order."""
    logger.info("No subcommand specified — running full pipeline.")
    run_annotate(args)
    run_filter(args)
    run_cluster(args)
    run_marker_genes(args)


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    # ── Root parser ──────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        prog="scrna_pipeline",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )

    # Global arguments (available in all subcommands and in the default full-pipeline run)
    global_group = parser.add_argument_group("global options")

    global_group.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Path to input matrixes (*.h5ad)"
    )

    global_group.add_argument(
        "-o", "--output",
        required=True,
        help="Output name"
    )

    global_group.add_argument(
        "-feature_type",
        default="gene",
        choices=["gene", "transcript"],
        help="Feature type to annotate (default: gene).",
    )

    # ── Subparsers (mutually exclusive by design) ─────────────────────────
    subparsers = parser.add_subparsers(
        dest="subcommand",
        metavar="SUBCOMMAND",
        help="Step to execute. Omit to run all steps in order.",
    )

    # ── annotate ──────────────────────────────────────────────────────────
    p_annotate = subparsers.add_parser(
        "annotate_features",
        help="Annotate features from input matrices.",
        description="Annotate features from one or more count matrices.",
    )

    p_annotate.add_argument(
        "--platform_type",
        default="10X",
        choices=["10X", "Argentag", "Parse"],
        help="Platform type (default: 10X).",
    )

    p_annotate.add_argument(
        "-st", "--specie_threshold",
        default=0.9,
        help="Threshold to assign cell to specie"
    )

    p_annotate.add_argument(
        "--genomes",
        nargs="+",
        default=["hsapiens_gene_ensembl"],
        help="Genomes to use for feature annotation. Available: hsapiens_gene_ensembl (human), mmusculus_gene_ensembl (mouse), clfamiliaris_gene_ensembl (dog). Default: hsapiens_gene_ensembl",
    )

    p_annotate.add_argument(
        "--condition",
        nargs="+",
        default=None,
        help="Condition name(s). Can be: single value (applied to all matrices), or N values (one per matrix). Default: extracted from filenames.",
    )

    p_annotate.add_argument(
        "--replicate",
        nargs="+",
        default=None,
        help="Replicate identifier(s). Can be: single value (applied to all matrices), or N values (one per matrix). Default: extracted from filenames.",
    )

    p_annotate.set_defaults(func=run_annotate)

    # ── filter ────────────────────────────────────────────────────────────
    p_filter = subparsers.add_parser(
        "filter",
        help="Filter cells and features.",
        description="Apply quality-control filters to an annotated dataset.",
    )

    p_filter.add_argument(
        "--min_genes",
        type=int,
        default=1000,
        help="Minimum number of genes per cell (default: 1000).",
    )

    p_filter.add_argument(
        "--max_gene",
        type=int,
        default=5000,
        help="Maximum number of genes per cell (default: 5000).",
    )

    p_filter.add_argument(
        "--min_cells_gene",
        type=int,
        default=15,
        help="Minimum number of cells per gene (default: 15).",
    )

    p_filter.add_argument(
        "--max_total_counts",
        type=int,
        default=20000,
        help="Maximum counts per cell (default:20000).",
    )

    p_filter.add_argument(
        "--max_mt_pct",
        type=float,
        default=15.0,
        help="Maximum mitochondrial gene percentage per cell (default: 15).",
    )

    p_filter.set_defaults(func=run_filter)

    # ── cluster ───────────────────────────────────────────────────────────
    p_cluster = subparsers.add_parser(
        "cluster",
        help="Cluster cells.",
        description="Dimensionality reduction and community detection.",
    )
    p_cluster.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        metavar="N",
        help="Number of principal components (default: 30).",
    )
    p_cluster.add_argument(
        "--algorithm",
        default="leiden",
        choices=["leiden", "louvain"],
        help="Community detection algorithm (default: leiden).",
    )

    p_cluster.add_argument(
        "--by_condition",
        action="store_true",
        help="Clustering at condition level instead of sample.",
    )

    p_cluster.set_defaults(func=run_cluster)

    # ── marker genes ───────────────────────────────────────────────────────────
    p_marker = subparsers.add_parser(
        "marker",
        help="Differential expression analysis.",
        description="Find marker genes between groups of cells.",
    )
    p_marker.add_argument(
        "--clustering_resolution",
        type=str,
        help="Clustering resolution to take into account for marker detection.",
    )
    p_marker.add_argument(
        "--method",
        default="wilcoxon",
        choices=["wilcoxon", "t-test", "logreg"],
        help="Statistical test method (default: wilcoxon).",
    )
    p_marker.add_argument(
        "--groupby",
        default="leiden",
        metavar="OBS_KEY",
        help="Observation key to group cells by (default: leiden).",
    )
    p_marker.add_argument(
        "--n-top-genes",
        type=int,
        default=100,
        metavar="N",
        help="Number of top marker genes to report per group (default: 100).",
    )
    p_marker.set_defaults(func=run_marker_genes)

    # ── diffexp ───────────────────────────────────────────────────────────
    #grouping_variable, sample_column, min_cells
    p_diffexp = subparsers.add_parser(
        "diffexp",
        help="Differential expression analysis.",
        description="Run pseudobulk and differential expression analysis (EdgeR).",
    )
    p_diffexp.add_argument(
        "--grouping",
        type=str,
        help="Variable to based the grouping on the pseudobulk on.",
    )
    p_diffexp.add_argument(
        "--comparisons",
        nargs=2,
        metavar=("CLUSTER_A", "CLUSTER_B"),
        action="append",
        default=None,
        help="Pairwise comparisons that must be run.",
    )
    p_diffexp.add_argument(
        "--sample_column",
        default="sample",
        help="Variable where the sample information is stored (default: sample).",
    )
    p_diffexp.add_argument(
        "--min_cells",
        type=int,
        default=20,
        help="Minimal number of cells required in a group to do pseudobulk aggregation (default: 20).",
    )
    p_diffexp.add_argument(
        "--log_fc_threshold",
        type=float,
        default=1,
        help="Log2 fold-change threshold for filtering DE genes (default: 1).",
    )
    p_diffexp.add_argument(
        "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold for filtering DE genes (default: 0.05).",
    )
    p_diffexp.add_argument(
        "--logcpm_threshold",
        type=float,
        default=1.0,
        help="logCPM threshold to filter out very lowly expressed genes (default: 1.0).",
    )
    p_diffexp.set_defaults(func=run_diffexp)

    # ── merge ───────────────────────────────────────────────────────────

    p_merge = subparsers.add_parser(
        "merge",
        help="Merge multiple input h5ad matrixes.",
        description="Merge multiple input h5ad matrixes.",
    )

    p_merge.add_argument(
        "--name",
        type=str,
        help="Name to give to the merged matrix.",
    )

    p_merge.set_defaults(func=run_merge)

    # ── umap plotting ───────────────────────────────────────────────────

    p_umap_plotting = subparsers.add_parser(
        "umap_plotting",
        help="Plotting the expression of specific features in the generated UMAP.",
        description="Plotting feature's expression.",
    )

    p_umap_plotting.add_argument(
        "--features",
        type=str,
        nargs='+',
        help="Features to plot their expression in the UMAP.",
    )

    p_umap_plotting.set_defaults(func=run_umap_plotting)

    # ── full-pipeline defaults (used when no subcommand is given) ─────────
    # These mirror the annotate arguments so that run_all() can forward them.
    full_group = parser.add_argument_group(
        "full-pipeline options (used when no subcommand is given)"
    )
    '''
    full_group.add_argument(
        "--feature-type",
        default="gene",
        choices=["gene", "peak", "protein"],
        help="Feature type to annotate (default: gene).",
    )
    full_group.add_argument(
        "--platform",
        default="illumina",
        choices=["illumina", "pacbio", "nanopore"],
        help="Sequencing platform (default: illumina).",
    )
    '''
    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    #Define input matrixes:
    input_matrixes = args.input

    # Create output directory
    Path(args.output).mkdir(parents=True, exist_ok=True)

    if args.subcommand is None:
        # No subcommand → full pipeline
        if not args.input:
            parser.error(
                "At least one --input matrix is required when running the full pipeline."
            )
        run_all(args)
    else:
        # Dispatch to the subcommand's registered function
        args.func(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
