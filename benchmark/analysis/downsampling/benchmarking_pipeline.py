"""
Platform benchmarking: Downsampling from raw counts and measuring proliferation resolution.
Uses actual column names from the data.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from scipy.stats import mannwhitneyu
import harmonypy


def create_proliferation_binary(adata_mouse):
    """
    Create binary proliferation label (once, on full depth).
    Uses existing phase labels: high = S+G2M, low = G1
    """
    adata_mouse.obs['proliferation_binary'] = (
        adata_mouse.obs['phase'].isin(['S', 'G2M'])
    ).astype(str)  # 'True' = S/G2M, 'False' = G1

    return adata_mouse


def downsample_adata(adata, target_umis, seed=None):
    """
    Downsample raw counts to target UMI depth per cell using scanpy's downsample_counts.

    Parameters:
    -----------
    adata : AnnData
        Input with raw counts in X
    target_umis : int
        Target total count per cell
    seed : int
        Random seed for reproducibility

    Returns:
    --------
    adata_ds : AnnData
        Downsampled copy
    """
    adata_ds = adata.copy()

    # Use scanpy's downsample_counts: samples each cell to target count
    sc.pp.downsample_counts(adata_ds, counts_per_cell=target_umis, random_state=seed)

    # Update n_counts metadata
    adata_ds.obs['n_counts'] = np.asarray(adata_ds.X.sum(axis=1)).ravel()

    return adata_ds


def preprocess_pipeline(adata, seed=None, batch_col='sample', norm_method='log1p', target_sum=1e4, use_harmony=True):
    """
    Shared preprocessing: QC → Normalize → Score genes → HVG → PCA → Neighbors → UMAP → [optional Harmony].
    Returns adata ready for clustering.

    Parameters:
    -----------
    norm_method : str
        Normalization method: 'log1p' (default), 'clr', or 'none'
    target_sum : float
        Target sum for normalize_total (default: 1e4 = 10K UMIs)
    use_harmony : bool
        Whether to apply harmony batch correction (default: True)
    """
    adata = adata.copy()

    # **QC**: Remove cells with no counts after downsampling
    adata = adata[adata.obs['n_counts'] > 0]

    if adata.n_obs == 0:
        raise ValueError("No cells left after QC!")

    # **Filter genes**: Remove genes with zero counts (important for downsampled data)
    sc.pp.filter_genes(adata, min_counts=1)

    # **Normalize** to target count per cell, then apply transform
    sc.pp.normalize_total(adata, target_sum=target_sum)

    if norm_method == 'log1p':
        sc.pp.log1p(adata)
    elif norm_method == 'clr':
        # CLR (centered log-ratio) normalization
        from scipy.stats import gmean
        log_counts = np.log(adata.X.data)
        geom_mean = gmean(adata.X.data)
        adata.X = adata.X / geom_mean
        sc.pp.log1p(adata)
    elif norm_method == 'none':
        pass  # No additional transform, just total count normalization

    # **Cell cycle scoring** (on normalized data)
    if "gene_id" in adata.var.columns:
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

        gene_names_set = set(adata.var["gene_id"].dropna())
        s_genes_found = [g for g in s_genes if g in gene_names_set]
        g2m_genes_found = [g for g in g2m_genes if g in gene_names_set]

        print(f"    Cell cycle genes available: S={len(s_genes_found)}, G2M={len(g2m_genes_found)}")

        if len(s_genes_found) > 0 and len(g2m_genes_found) > 0:
            # Create temporary adata with gene_id as var_names for scoring
            adata_temp = sc.AnnData(
                X=adata.X,
                obs=adata.obs.copy(),
                var=pd.DataFrame(index=adata.var["gene_id"].astype(str).values)
            )
            adata_temp.var_names_make_unique()

            # Score (no need to re-normalize, already done)
            sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_genes_found, g2m_genes=g2m_genes_found)

            # Copy scores back
            adata.obs["S_score"] = adata_temp.obs["S_score"].values
            adata.obs["G2M_score"] = adata_temp.obs["G2M_score"].values

    # **HVG selection** (10% of non-MT features, flavor='seurat')
    adata_without_mt = adata[:, ~adata.var["mt"]].copy() if "mt" in adata.var.columns else adata.copy()
    n_hvg = max(1, int(adata_without_mt.n_vars * 0.1))  # 10% of non-MT features, minimum 1
    sc.pp.highly_variable_genes(adata_without_mt, n_top_genes=n_hvg, flavor='seurat', batch_key=batch_col)

    # Transfer HVG mask back to the full object
    adata.var["highly_variable"] = False
    adata.var.loc[adata_without_mt.var_names, "highly_variable"] = adata_without_mt.var["highly_variable"].values
    adata = adata[:, adata.var['highly_variable']]

    print(f"    HVGs selected: {adata.n_vars} genes (10% of {adata_without_mt.n_vars} non-MT genes = {n_hvg} target)")

    # **PCA** (before harmony, on full feature set)
    sc.tl.pca(adata, n_comps=100)

    # **Neighbors & UMAP** (before harmony, on standard PCA)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, random_state=0)

    # **Harmony batch correction** (if enabled and batch column exists)
    if use_harmony and batch_col in adata.obs.columns:
        try:
            sc.external.pp.harmony_integrate(
                adata,
                key=batch_col,
                basis='X_pca',
                adjusted_basis='X_pca_harmony',

            )
            # Recalculate neighbors and UMAP on harmony-corrected PCs
            sc.pp.neighbors(adata, use_rep='X_pca_harmony', key_added='harmony', n_neighbors=15, n_pcs=50)
            sc.tl.umap(adata, neighbors_key='harmony', random_state=0)
        except Exception as e:
            print(f"Warning: Harmony failed ({e}), using standard PCA")

    return adata


def apply_genome_pipeline(adata, seed=None, batch_col='sample', norm_method='log1p', target_sum=1e4, use_harmony=True):
    """
    Genome analysis: All cells, Leiden clustering (resolution=0.05 for ~3 clusters).
    Expected clusters: 1 human, 1 dog, 1 mouse
    """
    adata = preprocess_pipeline(adata, seed=seed, batch_col=batch_col, norm_method=norm_method, target_sum=target_sum, use_harmony=use_harmony)

    # Leiden clustering on neighborhood graph (respects UMAP structure)
    sc.tl.leiden(adata, resolution=0.05, random_state=seed)

    return adata


def count_detected_genes_and_isoforms(adata):
    """
    Count detected genes and isoforms based on gene_id and transcript_id columns.
    A gene/isoform is detected if the sum of its counts across all cells > 0.
    """
    metrics = {}

    # Get raw counts
    X = adata.X

    # Sum counts for each variable
    var_counts = np.asarray(X.sum(axis=0)).ravel()

    # Count detected genes (unique gene_ids with counts > 0)
    if 'gene_id' in adata.var.columns:
        gene_id_counts = {}
        for idx, gene_id in enumerate(adata.var['gene_id']):
            if gene_id not in gene_id_counts:
                gene_id_counts[gene_id] = 0
            gene_id_counts[gene_id] += var_counts[idx]

        detected_genes = sum(1 for count in gene_id_counts.values() if count > 0)
        metrics['n_genes_detected'] = detected_genes

    # Count detected isoforms (unique transcript_ids with counts > 0)
    if 'transcript_id' in adata.var.columns:
        transcript_id_counts = {}
        for idx, transcript_id in enumerate(adata.var['transcript_id']):
            if transcript_id not in transcript_id_counts:
                transcript_id_counts[transcript_id] = 0
            transcript_id_counts[transcript_id] += var_counts[idx]

        detected_isoforms = sum(1 for count in transcript_id_counts.values() if count > 0)
        metrics['n_isoforms_detected'] = detected_isoforms

    return metrics


def calculate_knn_purity(adata, label_col='leiden', k=15):
    """
    Calculate KNN purity: mean fraction of k-nearest neighbors with the same label.
    Uses harmony-corrected PCA if available, otherwise standard PCA.
    """
    # Get the representation to use
    pca_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    X = adata.obsm[pca_rep]

    # Fit KNN model
    knn = NearestNeighbors(n_neighbors=k+1)  # +1 because it includes the cell itself
    knn.fit(X)
    distances, indices = knn.kneighbors(X)

    # Get labels
    labels = adata.obs[label_col].values.astype(str)

    # Calculate purity for each cell
    purities = []
    for i in range(adata.n_obs):
        neighbor_indices = indices[i][1:]  # Exclude the cell itself (first neighbor)
        neighbor_labels = labels[neighbor_indices]
        cell_label = labels[i]

        # Purity: fraction of neighbors with same label
        purity = (neighbor_labels == cell_label).sum() / len(neighbor_labels)
        purities.append(purity)

    return np.mean(purities)


def evaluate_genome_clustering(adata):
    """
    Evaluate genome analysis clustering (k=3, all cells).
    Metrics: ARI_genome, Silhouette, KNN_purity
    """
    metrics = {}

    ari_genome = adjusted_rand_score(
        adata.obs['specie_assignment'],
        adata.obs['leiden']
    )
    metrics['ari_genome'] = ari_genome

    # Use harmony-corrected PCA if available for silhouette, otherwise standard PCA
    pca_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'

    # Silhouette requires at least 2 clusters
    n_clusters = len(adata.obs['leiden'].unique())
    if n_clusters >= 2:
        sil = silhouette_score(adata.obsm[pca_rep], adata.obs['leiden'].astype(int))
    else:
        sil = np.nan  # Only 1 cluster found (clustering failed)
    metrics['silhouette'] = sil

    # Calculate KNN purity (also requires multiple clusters to be meaningful)
    if n_clusters >= 2:
        knn_purity = calculate_knn_purity(adata, label_col='leiden', k=15)
    else:
        knn_purity = np.nan
    metrics['knn_purity'] = knn_purity

    metrics['median_umis'] = adata.obs['n_counts'].median()
    metrics['n_features'] = adata.n_vars

    return metrics


def evaluate_proliferation_clustering(adata):
    """
    Evaluate proliferation signal in mouse cells (True vs False phase labels).
    Compares S_score + G2M_score between proliferating and non-proliferating cells.

    Metrics:
    - prolif_pvalue: Mann-Whitney U p-value (True vs False)
    - prolif_effect_size: Rank-Biserial correlation
    - prolif_true_median/mean: Proliferation score in True cells
    - prolif_false_median/mean: Proliferation score in False cells
    - prolif_fold_change: Mean True / Mean False
    - prolif_n_true/n_false: Cell counts
    - silhouette: Cluster compactness
    """
    metrics = {}

    adata.obs['prolif_score'] = adata.obs['S_score'] + adata.obs['G2M_score']

    true_mask = adata.obs['proliferation_binary'] == 'True'
    false_mask = adata.obs['proliferation_binary'] == 'False'

    true_scores = adata.obs.loc[true_mask, 'prolif_score'].values
    false_scores = adata.obs.loc[false_mask, 'prolif_score'].values

    if len(true_scores) > 0 and len(false_scores) > 0:
        stat, pval = mannwhitneyu(true_scores, false_scores)
        metrics['prolif_pvalue'] = pval

        n1, n2 = len(true_scores), len(false_scores)
        metrics['prolif_effect_size'] = 1 - (2 * stat / (n1 * n2))

        metrics['prolif_true_median'] = np.median(true_scores)
        metrics['prolif_false_median'] = np.median(false_scores)
        metrics['prolif_true_mean'] = np.mean(true_scores)
        metrics['prolif_false_mean'] = np.mean(false_scores)

        if metrics['prolif_false_mean'] > 0:
            metrics['prolif_fold_change'] = metrics['prolif_true_mean'] / metrics['prolif_false_mean']
        else:
            metrics['prolif_fold_change'] = np.nan

        metrics['prolif_n_true'] = n1
        metrics['prolif_n_false'] = n2
    else:
        metrics['prolif_pvalue'] = np.nan
        metrics['prolif_effect_size'] = np.nan
        metrics['prolif_true_median'] = np.nan
        metrics['prolif_false_median'] = np.nan
        metrics['prolif_true_mean'] = np.nan
        metrics['prolif_false_mean'] = np.nan
        metrics['prolif_fold_change'] = np.nan
        metrics['prolif_n_true'] = len(true_scores)
        metrics['prolif_n_false'] = len(false_scores)

    metrics['median_umis'] = adata.obs['n_counts'].median()
    metrics['n_features'] = adata.n_vars

    return metrics


def generate_depth_targets(min_depth, max_depth, n_points=4):
    """
    Generate evenly-spaced target depths between min and max.

    Parameters:
    -----------
    min_depth : int
        Minimum target UMI depth
    max_depth : int
        Maximum target UMI depth
    n_points : int
        Number of depth points to generate (default: 4)

    Returns:
    --------
    target_depths : list
        Evenly-spaced depth levels
    """
    return sorted([int(d) for d in np.linspace(min_depth, max_depth, n_points)])


def compute_downsample_targets(adata_dict, custom_depths=None, depth_fractions=None):
    """
    Find minimum median UMI depth across all platforms.
    Define target depths as fractions of that minimum, or use custom depths.

    Parameters:
    -----------
    adata_dict : dict
        Dictionary of platform names to adata objects
    custom_depths : list, optional
        Manually specified depth levels. If provided, uses these instead of computing from fractions.
    depth_fractions : list, optional
        Custom fractions of minimum median to use (default: [0.25, 0.50, 0.75, 1.0])

    Returns:
    --------
    target_depths : list
        Depth levels to test
    min_median : float
        Minimum median UMIs across platforms
    original_medians : dict
        Original median UMIs per platform
    """
    medians = {}
    for platform, adata in adata_dict.items():
        # Raw counts in layers
        X_raw = adata.layers['raw_counts']
        median_umis = np.median(np.asarray(X_raw.sum(axis=1)).ravel())
        medians[platform] = median_umis

    min_median = min(medians.values())

    # Use custom depths if provided, otherwise compute from fractions
    if custom_depths is not None:
        target_depths = sorted(custom_depths)
        print(f"\nOriginal Sequencing Depth (Median UMIs/cell):")
        for platform in sorted(medians.keys()):
            print(f"  {platform:15s}: {medians[platform]:7.0f}")
        print(f"\nDownsampling Design:")
        print(f"  Using custom target depths: {target_depths}")
    else:
        if depth_fractions is None:
            depth_fractions = [0.25, 0.50, 0.75, 1.0]
        target_depths = [int(min_median * frac) for frac in depth_fractions]

        print(f"\nOriginal Sequencing Depth (Median UMIs/cell):")
        for platform in sorted(medians.keys()):
            print(f"  {platform:15s}: {medians[platform]:7.0f}")

        print(f"\nDownsampling Design:")
        print(f"  Reference (min median): {min_median:.0f} UMIs/cell")
        print(f"  Target depths:")
        for frac, depth in zip(depth_fractions, target_depths):
            print(f"    {frac*100:5.0f}% of ref: {depth:6.0f} UMIs/cell")

    return target_depths, min_median, medians


def benchmark_platforms(matrix_paths_dict, target_depths, n_seeds=3, batch_col='batch', norm_method='log1p', target_sum=1e4, save_umap_data=True, umap_data_dir='.umap_temp', use_harmony=True):
    """
    Full benchmarking loop: all platforms × all matrix replicates × all depths × all seeds.

    Parameters:
    -----------
    matrix_paths_dict : dict
        {'platform1': [path1, path2, ...], 'platform2': [...], ...}
        Each list contains paths to replicate matrices for that platform.
        Each matrix should have:
        - X: empty (will be filled from layers['raw_counts'])
        - layers['raw_counts']: raw UMI counts
        - obs['specie_assignment']: genome/specie label
        - obs['S_score'], obs['G2M_score']: cell cycle scores
        - obs['phase']: phase label (G1, S, G2M)
        - obs['batch']: (optional) batch/platform identifier for harmony

    target_depths : list
        Depth levels to test (UMIs per cell)

    n_seeds : int
        Number of random seeds per depth level (per matrix)

    batch_col : str
        Column name for batch correction in harmony (default: 'batch')

    Returns:
    --------
    results_df : DataFrame
        All results (platform × matrix_seed × depth × downsample_seed × metrics)

    summary : DataFrame
        Aggregated mean ± IQR per platform × depth (across all matrix replicates and downsampling seeds)
    """
    all_results = []

    # Create directory for saving preprocessed adata objects for UMAP generation
    if save_umap_data:
        import os
        os.makedirs(umap_data_dir, exist_ok=True)

    for platform_name, matrix_paths in matrix_paths_dict.items():
        print(f"\n{'='*60}")
        print(f"Processing: {platform_name}")
        print(f"{'='*60}")

        # Loop over matrix replicates (seed0, seed1, seed2)
        for matrix_seed, matrix_path in enumerate(matrix_paths):
            print(f"\n  Matrix replicate {matrix_seed}...")

            # Load matrix
            adata_full = sc.read_h5ad(matrix_path)

            # **IMPORTANT: Set X to raw counts (original h5ad has log-transformed data in X)**
            adata_full.X = adata_full.layers['raw_counts'].copy()

            # **Create proliferation binary label ONCE on full depth (mouse cells only)**
            mouse_mask = adata_full.obs['specie_assignment'] == 'mouse'
            adata_full.obs['proliferation_binary'] = 'False'  # Default to False
            adata_full.obs.loc[mouse_mask, 'proliferation_binary'] = (
                adata_full.obs.loc[mouse_mask, 'phase'].isin(['S', 'G2M'])
            ).astype(str)

            n_mouse = mouse_mask.sum()
            print(f"  Mouse cells: {n_mouse}")

            # Get full depth value from all cells
            actual_full_depth = int(np.median(np.asarray(adata_full.layers['raw_counts'].sum(axis=1)).ravel()))

            print(f"  Actual full depth (median UMIs): {actual_full_depth}")
            print(f"  Target depths to downsample: {sorted(target_depths)}")

            # First, run all target depths with downsampling
            for depth in sorted(target_depths):
                for seed in range(n_seeds):
                    print(f"  Depth {depth:5.0f} UMIs, seed {seed}... ", end='', flush=True)

                    try:
                        # Always downsample to target depth
                        adata_full_ds = downsample_adata(adata_full, target_umis=depth, seed=seed)

                        # Count detected genes and isoforms before preprocessing
                        detection_metrics = count_detected_genes_and_isoforms(adata_full_ds)

                        # Preprocess once for both analyses
                        adata_preprocessed = preprocess_pipeline(
                            adata_full_ds,
                            seed=seed,
                            batch_col=batch_col,
                            norm_method=norm_method,
                            target_sum=target_sum,
                            use_harmony=use_harmony
                        )

                        # -------- GENOME ANALYSIS: All cells, Leiden clustering (resolution=0.05) --------
                        adata_genome_analyzed = adata_preprocessed.copy()
                        sc.tl.leiden(adata_genome_analyzed, resolution=0.05, random_state=seed)

                        # Save preprocessed adata with clustering for UMAP generation
                        if save_umap_data:
                            import os
                            safe_platform = platform_name.replace('-', '_').replace('+', 'plus').lower()
                            # Store metadata for UMAP display
                            adata_genome_analyzed.uns['depth_umis'] = int(depth)
                            adata_genome_analyzed.uns['n_hvgs'] = adata_genome_analyzed.n_vars
                            save_path = os.path.join(umap_data_dir, f'{safe_platform}_{int(depth)}_m{matrix_seed}_s{seed}.h5ad')
                            adata_genome_analyzed.write(save_path)

                        metrics_genome = evaluate_genome_clustering(adata_genome_analyzed)
                        metrics_genome.update(detection_metrics)  # Add detected genes and isoforms
                        metrics_genome['analysis_type'] = 'genome'
                        metrics_genome['platform'] = platform_name
                        metrics_genome['matrix_seed'] = matrix_seed
                        metrics_genome['depth'] = depth
                        metrics_genome['downsample_seed'] = seed
                        metrics_genome['is_full_depth'] = False  # These are downsampled depths
                        all_results.append(metrics_genome)

                        # -------- PROLIFERATION ANALYSIS: Mouse only --------
                        adata_prolif_analyzed = adata_preprocessed[adata_preprocessed.obs['specie_assignment'] == 'mouse'].copy()

                        metrics_prolif = evaluate_proliferation_clustering(adata_prolif_analyzed)
                        metrics_prolif['analysis_type'] = 'proliferation'
                        metrics_prolif['platform'] = platform_name
                        metrics_prolif['matrix_seed'] = matrix_seed
                        metrics_prolif['depth'] = depth
                        metrics_prolif['downsample_seed'] = seed
                        metrics_prolif['is_full_depth'] = False  # These are downsampled depths
                        all_results.append(metrics_prolif)

                        pval_str = f"{metrics_prolif['prolif_pvalue']:.2e}" if metrics_prolif['prolif_pvalue'] > 0 else "p<1e-15"
                        print(f"✓ p={pval_str}, fold={metrics_prolif['prolif_fold_change']:.1f}x, ARI_gen={metrics_genome['ari_genome']:.3f}")

                    except Exception as e:
                        print(f"✗ Error: {e}")
                        continue

            # Second, run full depth (no downsampling) for all seeds
            print(f"\n  Running FULL DEPTH (no downsampling)...")
            for seed in range(n_seeds):
                print(f"  Full depth, seed {seed}... ", end='', flush=True)

                try:
                    # Use full dataset without downsampling
                    adata_full_ds = adata_full.copy()
                    adata_full_ds.obs['n_counts'] = np.asarray(adata_full_ds.X.sum(axis=1)).ravel()

                    # Count detected genes and isoforms before preprocessing
                    detection_metrics = count_detected_genes_and_isoforms(adata_full_ds)

                    # Preprocess once for both analyses
                    adata_preprocessed = preprocess_pipeline(
                        adata_full_ds,
                        seed=seed,
                        batch_col=batch_col,
                        norm_method=norm_method,
                        target_sum=target_sum,
                        use_harmony=use_harmony
                    )

                    # -------- GENOME ANALYSIS: All cells, Leiden clustering (resolution=0.05) --------
                    adata_genome_analyzed = adata_preprocessed.copy()
                    sc.tl.leiden(adata_genome_analyzed, resolution=0.05, random_state=seed)

                    # Save preprocessed adata with clustering for UMAP generation
                    if save_umap_data:
                        import os
                        safe_platform = platform_name.replace('-', '_').replace('+', 'plus').lower()
                        # Store metadata for UMAP display
                        adata_genome_analyzed.uns['depth_umis'] = actual_full_depth
                        adata_genome_analyzed.uns['n_hvgs'] = adata_genome_analyzed.n_vars
                        save_path = os.path.join(umap_data_dir, f'{safe_platform}_full_depth_m{matrix_seed}_s{seed}.h5ad')
                        adata_genome_analyzed.write(save_path)

                    metrics_genome = evaluate_genome_clustering(adata_genome_analyzed)
                    metrics_genome.update(detection_metrics)  # Add detected genes and isoforms
                    metrics_genome['analysis_type'] = 'genome'
                    metrics_genome['platform'] = platform_name
                    metrics_genome['matrix_seed'] = matrix_seed
                    metrics_genome['depth'] = actual_full_depth
                    metrics_genome['downsample_seed'] = seed
                    metrics_genome['is_full_depth'] = True
                    all_results.append(metrics_genome)

                    # -------- PROLIFERATION ANALYSIS: Mouse only --------
                    adata_prolif_analyzed = adata_preprocessed[adata_preprocessed.obs['specie_assignment'] == 'mouse'].copy()

                    metrics_prolif = evaluate_proliferation_clustering(adata_prolif_analyzed)
                    metrics_prolif['analysis_type'] = 'proliferation'
                    metrics_prolif['platform'] = platform_name
                    metrics_prolif['matrix_seed'] = matrix_seed
                    metrics_prolif['depth'] = actual_full_depth
                    metrics_prolif['downsample_seed'] = seed
                    metrics_prolif['is_full_depth'] = True
                    all_results.append(metrics_prolif)

                    pval_str = f"{metrics_prolif['prolif_pvalue']:.2e}" if metrics_prolif['prolif_pvalue'] > 0 else "p<1e-15"
                    print(f"✓ p={pval_str}, fold={metrics_prolif['prolif_fold_change']:.1f}x, ARI_gen={metrics_genome['ari_genome']:.3f}")

                except Exception as e:
                    print(f"✗ Error: {e}")
                    continue

    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)

    if len(results_df) == 0:
        print("ERROR: No results collected. Check for exceptions in the main loop.")
        return None, None

    # Aggregate: median ± IQR per analysis_type × platform × depth
    summaries = []

    genome_results = results_df[results_df['analysis_type'] == 'genome']
    if len(genome_results) > 0:
        summary_genome = genome_results.groupby(['platform', 'depth']).agg({
            'ari_genome': ['median', ('q1', lambda x: x.quantile(0.25)), ('q3', lambda x: x.quantile(0.75))],
            'silhouette': ['median', ('q1', lambda x: x.quantile(0.25)), ('q3', lambda x: x.quantile(0.75))],
            'median_umis': ['median'],
        }).round(3)
        summaries.append(summary_genome.assign(analysis_type='genome'))

    prolif_results = results_df[results_df['analysis_type'] == 'proliferation']
    if len(prolif_results) > 0:
        summary_prolif = prolif_results.groupby(['platform', 'depth']).agg({
            'prolif_pvalue': ['median', ('q1', lambda x: x.quantile(0.25)), ('q3', lambda x: x.quantile(0.75))],
            'prolif_effect_size': ['median', ('q1', lambda x: x.quantile(0.25)), ('q3', lambda x: x.quantile(0.75))],
            'prolif_fold_change': ['median', ('q1', lambda x: x.quantile(0.25)), ('q3', lambda x: x.quantile(0.75))],
            'prolif_true_median': ['median'],
            'prolif_false_median': ['median'],
            'median_umis': ['median'],
        }).round(3)
        summaries.append(summary_prolif.assign(analysis_type='proliferation'))

    # Combine summaries
    summary = pd.concat(summaries) if summaries else pd.DataFrame()

    return results_df, summary


def generate_genome_umaps_from_saved(platform_name, target_depths, n_seeds=3, n_matrix_seeds=3, umap_data_dir='.umap_temp', cleanup=True, include_full_depth=True):
    """
    Generate UMAP grid for genome analysis from pre-saved adata objects.
    Loads preprocessed adata files saved during benchmarking.

    Parameters:
    -----------
    platform_name : str
        Platform name (for figure title)
    target_depths : list or array
        Depth levels to visualize (downsampled depths)
    n_seeds : int
        Number of downsample seeds
    n_matrix_seeds : int
        Number of input matrix replicates
    umap_data_dir : str
        Directory containing saved adata objects
    cleanup : bool
        Whether to delete saved files after generating figures
    include_full_depth : bool
        Whether to include full depth (no downsampling) as the last row

    Returns:
    --------
    fig : matplotlib figure
        Grid of UMAP plots (depths × [matrix_seeds × downsample_seeds])
    """
    import matplotlib.pyplot as plt
    import os
    import glob

    safe_platform = platform_name.replace('-', '_').replace('+', 'plus').lower()
    depths_sorted = sorted(set(target_depths))

    # Add full depth as the last row if requested and file exists
    if include_full_depth:
        full_depth_file = os.path.join(umap_data_dir, f'{safe_platform}_full_depth_m0_s0.h5ad')
        if os.path.exists(full_depth_file):
            depths_sorted = list(depths_sorted) + ['full_depth']

    n_depths = len(depths_sorted)
    n_cols = n_matrix_seeds * n_seeds  # Total columns: all matrix and downsample seeds

    # Create grid figure: rows = depths, columns = matrix_seeds × downsample_seeds
    fig, axes = plt.subplots(n_depths, n_cols, figsize=(4.5 * n_cols, 4.5 * n_depths))

    # Ensure 2D axes array
    if n_depths == 1:
        axes = axes.reshape(1, -1)
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for depth_idx, depth in enumerate(depths_sorted):
        # Add row label to the left
        if depth == 'full_depth':
            row_label = 'FULL DEPTH'
        else:
            row_label = f'{int(depth)} UMIs'

        fig.text(0.01, 0.95 - (depth_idx + 0.5) * (1.0 / n_depths), row_label,
                 fontsize=11, fontweight='bold', va='center',
                 bbox=dict(boxstyle='round', facecolor='#ffffcc' if depth == 'full_depth' else '#f0f0f0', alpha=0.7))
        for matrix_idx in range(n_matrix_seeds):
            for seed_idx in range(n_seeds):
                col_idx = matrix_idx * n_seeds + seed_idx
                ax = axes[depth_idx, col_idx]

                # Construct appropriate filename based on depth type
                if depth == 'full_depth':
                    depth_label = 'Full Depth'
                    adata_path = os.path.join(umap_data_dir, f'{safe_platform}_full_depth_m{matrix_idx}_s{seed_idx}.h5ad')
                else:
                    depth_label = f'{int(depth)}'
                    adata_path = os.path.join(umap_data_dir, f'{safe_platform}_{int(depth)}_m{matrix_idx}_s{seed_idx}.h5ad')

                print(f"    UMAP: {depth_label}, Matrix {matrix_idx}, Downsample {seed_idx}...", end=' ', flush=True)

                try:
                    # Load saved adata
                    if not os.path.exists(adata_path):
                        raise FileNotFoundError(f"Saved adata not found: {adata_path}")

                    adata_processed = sc.read_h5ad(adata_path)

                    # Ensure we're using harmony-corrected UMAP if available
                    if 'X_pca_harmony' in adata_processed.obsm and 'harmony' not in adata_processed.obsp:
                        # Recompute neighbors from harmony PCA if needed
                        sc.pp.neighbors(adata_processed, use_rep='X_pca_harmony', key_added='harmony', n_neighbors=15, n_pcs=50)
                        sc.tl.umap(adata_processed, neighbors_key='harmony')

                    # Plot UMAP colored by specie_assignment
                    sc.pl.umap(
                        adata_processed,
                        color='specie_assignment',
                        ax=ax,
                        show=False,
                        legend_fontsize=7,
                        title=''
                    )

                    # Compute metrics from loaded data
                    ari = adjusted_rand_score(
                        adata_processed.obs['specie_assignment'],
                        adata_processed.obs['leiden']
                    )
                    n_clusters = len(adata_processed.obs['leiden'].unique())

                    # Silhouette requires at least 2 clusters
                    if n_clusters >= 2:
                        sil = silhouette_score(
                            adata_processed.obsm['X_pca_harmony'] if 'X_pca_harmony' in adata_processed.obsm else adata_processed.obsm['X_pca'],
                            adata_processed.obs['leiden'].astype(int)
                        )
                    else:
                        sil = np.nan

                    # KNN purity
                    if n_clusters >= 2:
                        knn_pur = calculate_knn_purity(adata_processed, label_col='leiden', k=15)
                    else:
                        knn_pur = np.nan

                    # Get metadata from saved adata
                    depth = adata_processed.uns.get('depth_umis', 'N/A')
                    n_hvgs = adata_processed.uns.get('n_hvgs', 'N/A')

                    # Set title with metadata and metrics
                    title_str = f"Depth: {depth} UMIs | HVGs: {n_hvgs} | M{matrix_idx}S{seed_idx}\nARI: {ari:.3f}"
                    if not np.isnan(sil):
                        title_str += f" | SIL: {sil:.3f}"
                    if not np.isnan(knn_pur):
                        title_str += f" | KNN: {knn_pur:.3f}"

                    ax.set_title(title_str, fontsize=8, fontweight='bold')

                    # Clean up axes labels for interior subplots
                    if col_idx > 0:
                        ax.set_ylabel('')
                    if depth_idx < n_depths - 1:
                        ax.set_xlabel('')

                    print("✓")

                except Exception as e:
                    ax.text(0.5, 0.5, f"Error: {str(e)[:25]}", ha='center', va='center', transform=ax.transAxes, fontsize=8)
                    ax.set_title(f"M{matrix_idx}S{seed_idx}\n✗ Failed", fontsize=9, color='red')
                    print(f"✗ {e}")

    # Overall title
    fig.suptitle(
        f'Genome Analysis: {platform_name}\nUMAP colored by specie_assignment (M=Matrix seed, S=Downsample seed)',
        fontsize=13,
        fontweight='bold',
        y=0.995
    )

    plt.tight_layout()

    # Cleanup saved files if requested
    if cleanup:
        for saved_file in glob.glob(os.path.join(umap_data_dir, f'{safe_platform}_*.h5ad')):
            try:
                os.remove(saved_file)
            except Exception as e:
                print(f"Warning: Could not delete {saved_file}: {e}")

    return fig


if __name__ == "__main__":

    print("Benchmarking pipeline loaded. Use benchmark_platforms() to run analysis.")