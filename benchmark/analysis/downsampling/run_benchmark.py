"""
Run the benchmarking analysis and generate visualizations.

Can be re-run for plot generation only if benchmark_results.csv and
benchmark_summary.csv already exist.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from benchmarking_pipeline import (
    benchmark_platforms,
    compute_downsample_targets,
    generate_depth_targets,
    generate_genome_umaps_from_saved
)


def main():
    """Load data, run benchmark, save results, generate plots."""

    # ========================================
    # CONFIGURATION
    # ========================================
    USE_HARMONY = True  # Set to False to test without batch correction

    # PLATFORM COLORS (customize here)
    platform_colors = {
        '10X-3PRIME': '#658b66',
        '10X-5PRIME': '#8aa783',
        'ARGENTAG': '#2680a6',
        'PARSE': '#875670',
    }

    # Check if results already exist
    results_exist = os.path.exists("benchmark_results.csv") and os.path.exists("benchmark_summary.csv")

    if results_exist:
        print("Found existing benchmark results. Loading from disk...")
        results_df = pd.read_csv("benchmark_results.csv")
        summary = pd.read_csv("benchmark_summary.csv")
        print("✓ Loaded benchmark_results.csv and benchmark_summary.csv")

        # Extract target depths from loaded results (exclude full depth)
        genome_results_loaded = results_df[results_df['analysis_type'] == 'genome'].copy()
        if 'is_full_depth' in genome_results_loaded.columns:
            target_depths = sorted(genome_results_loaded[genome_results_loaded['is_full_depth'] == False]['depth'].unique())
        else:
            target_depths = sorted(genome_results_loaded['depth'].unique())
    else:
        # ========================================
        # 1. LOAD DATA
        # ========================================
        print("Loading data...")

        # Load multiple matrices per platform (seed0, seed1, seed2 as replicates)
        base_path = "<path_to_downsampled_matrices>"
        matrix_paths = {
            '10X-3PRIME': [
                f"{base_path}/10X-3PRIME_transcript_corrected_clustered_matrix_cell_cycle_scored_3000_cells_seed{i}.h5ad"
                for i in range(3)
            ],
            '10X-5PRIME': [
                f"{base_path}/10X-5PRIME_transcript_corrected_clustered_matrix_cell_cycle_scored_3000_cells_seed{i}.h5ad"
                for i in range(3)
            ],
            'ARGENTAG': [
                f"{base_path}/ARGENTAG_transcript_corrected_clustered_matrix_cell_cycle_scored_3000_cells_seed{i}.h5ad"
                for i in range(3)
            ],
            'PARSE': [
                f"{base_path}/PARSE_transcript_corrected_clustered_matrix_cell_cycle_scored_3000_cells_seed{i}.h5ad"
                for i in range(3)
            ],
        }

        # ========================================
        # 2. COMPUTE TARGET DEPTHS & ORIGINAL MEDIANS
        # ========================================
        print("\nComputing target depths...")

        # -------- CUSTOMIZE DOWNSAMPLING DEPTHS HERE --------

        # Load first matrix of each platform to compute target depths
        adata_dict_sample = {platform: sc.read_h5ad(paths[0]) for platform, paths in matrix_paths.items()}
        target_depths, min_median, original_medians = compute_downsample_targets(
            adata_dict_sample,
            custom_depths=[100, 500, 1000, 2000, 3000]
        )

        # ========================================
        # 3. RUN BENCHMARKING
        # ========================================
        print("\nRunning benchmarking analysis...")
        print("(This may take 10-30 minutes depending on data size)")

        results_df, summary = benchmark_platforms(
            matrix_paths,
            target_depths=target_depths,
            n_seeds=3,
            batch_col='replicate',  # Correct for within-platform replicates
            norm_method='log1p',    # Normalization: 'log1p' (default), 'clr', or 'none'
            target_sum=1e4,         # Normalize to 10K UMIs per cell
            save_umap_data=True,    # Save preprocessed data for UMAP generation
            umap_data_dir='.umap_temp',  # Directory for temporary adata files
            use_harmony=USE_HARMONY  # Batch correction: True (default) or False to test
        )

        # ========================================
        # 4. SAVE RESULTS
        # ========================================
        print("\nSaving results...")
        results_df.to_csv("benchmark_results.csv", index=False)

        # Flatten summary columns before saving for easier reloading
        summary_to_save = summary.reset_index()
        summary_flat_cols = []
        for col in summary_to_save.columns:
            if isinstance(col, tuple):
                metric, func = col
                if func == '':
                    summary_flat_cols.append(metric)
                else:
                    summary_flat_cols.append(f"{metric}_{func}")
            else:
                summary_flat_cols.append(col)
        summary_to_save.columns = summary_flat_cols
        summary_to_save.to_csv("benchmark_summary.csv", index=False)

        # Save genes and isoforms detected per depth and matrix
        genes_isoforms = results_df[
            (results_df['analysis_type'] == 'genome') &
            (results_df['n_genes_detected'].notna())
        ][['platform', 'matrix_seed', 'depth', 'n_genes_detected', 'n_isoforms_detected']].drop_duplicates()
        genes_isoforms = genes_isoforms.sort_values(['platform', 'depth', 'matrix_seed'])
        genes_isoforms.to_csv("genes_isoforms_detected.csv", index=False)
        print("✓ Saved: genes_isoforms_detected.csv")

        # Split and save genome and proliferation summaries separately
        genome_summary = summary[summary['analysis_type'] == 'genome'].reset_index()
        prolif_summary = summary[summary['analysis_type'] == 'proliferation'].reset_index()

        genome_summary.to_csv("genome_ari_summary.csv", index=False)
        prolif_summary.to_csv("proliferation_summary.csv", index=False)

        # Save original sequencing depths
        with open("original_sequencing_depths.txt", "w") as f:
            f.write("="*70 + "\n")
            f.write("ORIGINAL SEQUENCING DEPTH (Median UMIs/cell)\n")
            f.write("="*70 + "\n\n")
            for platform in sorted(original_medians.keys()):
                f.write(f"{platform:20s}: {original_medians[platform]:7.0f} UMIs/cell\n")
            f.write(f"\n{'='*70}\n")
            f.write("DOWNSAMPLING LEVELS (as % of minimum)\n")
            f.write(f"{'='*70}\n\n")
            f.write(f"Reference depth (min median): {min_median:.0f} UMIs/cell\n\n")
            for idx, depth in enumerate(sorted(set(target_depths))):
                pct = (depth / min_median) * 100
                f.write(f"{pct:5.0f}% →  {depth:6.0f} UMIs/cell\n")

        print(f"\n{'='*70}")
        print("ORIGINAL SEQUENCING DEPTH (Median UMIs/cell)")
        print(f"{'='*70}")
        for platform in sorted(original_medians.keys()):
            print(f"  {platform:15s}: {original_medians[platform]:7.0f}")

        print(f"\n{'='*70}")
        print("DOWNSAMPLING REFERENCE")
        print(f"{'='*70}")
        print(f"Minimum median (reference): {min_median:.0f} UMIs/cell\n")
        print("Target depths as % of reference:")
        for depth in sorted(set(target_depths)):
            pct = (depth / min_median) * 100
            print(f"  {pct:5.0f}% → {depth:6.0f} UMIs/cell")

        print(f"\n{'='*70}")
        print("BENCHMARK SUMMARY (Aggregated across seeds)")
        print(f"{'='*70}")
        print(summary)

    # ========================================
    # 5. GENERATE VISUALIZATIONS
    # ========================================
    print("\nGenerating visualizations...")

    # If summary is not yet flattened (freshly computed), flatten it
    if any(isinstance(col, tuple) for col in summary.columns):
        summary_flat = summary.reset_index()
        # Flatten multi-level columns
        new_cols = []
        for col in summary_flat.columns:
            if isinstance(col, tuple):
                metric, func = col
                if func == '':
                    new_cols.append(metric)
                else:
                    new_cols.append(f"{metric}_{func}")
            else:
                new_cols.append(col)
        summary_flat.columns = new_cols
    else:
        # Already flattened (loaded from CSV)
        summary_flat = summary

    # Exclude full_depth samples for consistent comparison across platforms
    if 'is_full_depth' in summary_flat.columns:
        summary_flat = summary_flat[summary_flat['is_full_depth'] == False].copy()
    elif 'is_full_depth' in results_df.columns:
        # Get depths to exclude (the full_depth for each platform)
        full_depths = results_df[results_df['is_full_depth'] == True]['depth'].unique()
        summary_flat = summary_flat[~summary_flat['depth'].isin(full_depths)].copy()

    platforms = sorted(results_df['platform'].unique())
    

    # --- Figure 4b: Proliferation Score Trajectories (True vs False) ---
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    for idx, platform in enumerate(platforms):
        ax = axes[idx]
        subset = summary_flat[
            (summary_flat['platform'] == platform) &
            (summary_flat['analysis_type'] == 'proliferation')
        ].sort_values('depth')

        depths = subset['depth'].values
        true_medians = subset['prolif_true_median_median'].values
        false_medians = subset['prolif_false_median_median'].values

        ax.plot(depths, true_medians, marker='o', label='True (S/G2M)', linewidth=2.5, markersize=8,
                color=platform_colors.get(platform, '#cccccc'), alpha=0.8)
        ax.plot(depths, false_medians, marker='s', label='False (G1)', linewidth=2.5, markersize=8,
                color=platform_colors.get(platform, '#cccccc'), alpha=0.4, linestyle='--')

        ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Median Proliferation Score', fontsize=11, fontweight='bold')
        ax.set_title(f'{platform}', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xscale('log')

    plt.suptitle('Proliferation Score Trajectory: True vs False Cells', fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig('fig4b_prolif_score_trajectory.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig4b_prolif_score_trajectory.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig4b_prolif_score_trajectory.{pdf,png}")

    # --- Figure 5b: Score Distributions Across Depths (Box Plots) ---
    print("\nGenerating proliferation score distribution plots...")

    prolif_full = results_df[results_df['analysis_type'] == 'proliferation'].copy()

    # Show all analyzed depths
    depths_to_show = sorted(set(target_depths))

    for platform in platforms:
        fig, axes = plt.subplots(1, len(depths_to_show), figsize=(5*len(depths_to_show), 6), sharey=True)
        if len(depths_to_show) == 1:
            axes = [axes]

        platform_data = prolif_full[prolif_full['platform'] == platform].sort_values('depth')

        for ax_idx, depth in enumerate(depths_to_show):
            ax = axes[ax_idx]
            depth_data = platform_data[platform_data['depth'] == depth]

            if len(depth_data) == 0:
                continue

            # For this depth, show the true/false values from the seeds
            true_vals = depth_data['prolif_true_median'].values
            false_vals = depth_data['prolif_false_median'].values

            # Prepare data for violin plot
            data_to_plot = [
                {'Score': val, 'Group': 'False (G1)'} for val in false_vals
            ] + [
                {'Score': val, 'Group': 'True (S/G2M)'} for val in true_vals
            ]
            plot_df = pd.DataFrame(data_to_plot)

            # Create violin plot
            sns.violinplot(data=plot_df, x='Group', y='Score', hue='Group', ax=ax, palette=['#aaaaaa', platform_colors.get(platform, '#cccccc')], inner='box', legend=False)

            ax.set_ylabel('Proliferation Score', fontsize=10, fontweight='bold')
            ax.set_xlabel('')
            ax.set_title(f'Depth: {int(depth)} UMIs\n(n={len(depth_data)} seeds)', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')

        fig.suptitle(f'{platform}: Score Distributions Across Depths', fontsize=13, fontweight='bold')
        plt.tight_layout()

        safe_name = platform.replace('-', '_').replace('+', 'plus').lower()
        plt.savefig(f'fig5b_prolif_distributions_{safe_name}.pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f'fig5b_prolif_distributions_{safe_name}.png', dpi=150, bbox_inches='tight')
        print(f"✓ Saved: fig5b_prolif_distributions_{safe_name}.{{pdf,png}}")
        plt.close(fig)

    # --- Figure 5c: Genome ARI (Barplot) ---
    print("\nGenerating ARI (genome) barplot...")

    genome_results = results_df[results_df['analysis_type'] == 'genome'].copy()

    # Exclude full depth samples
    if 'is_full_depth' in genome_results.columns:
        genome_results = genome_results[genome_results['is_full_depth'] == False]

    genome_agg = genome_results.groupby(['platform', 'depth']).agg({
        'ari_genome': ['mean', 'std']
    }).reset_index()
    genome_agg.columns = ['platform', 'depth', 'ari_mean', 'ari_std']

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(genome_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = genome_agg[genome_agg['platform'] == platform].sort_values('depth')
        means = subset['ari_mean'].values
        stds = subset['ari_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Mean ARI (Genome Separation)', fontsize=13, fontweight='bold')
    ax.set_title('Genome/Specie Separation vs Depth (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax.set_ylim([-0.05, 1.05])
    plt.tight_layout()
    plt.savefig('fig5c_ari_genome_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig5c_ari_genome_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig5c_ari_genome_mean_sd.{pdf,png}")

    # --- Figure 5d: Proliferation Effect Size (Barplot) ---
    print("\nGenerating absolute effect size (proliferation) barplot...")

    prolif_results = results_df[results_df['analysis_type'] == 'proliferation'].copy()

    # Exclude full depth samples
    if 'is_full_depth' in prolif_results.columns:
        prolif_results = prolif_results[prolif_results['is_full_depth'] == False]

    prolif_agg = prolif_results.groupby(['platform', 'depth']).agg({
        'prolif_effect_size': ['mean', 'std']
    }).reset_index()
    prolif_agg.columns = ['platform', 'depth', 'effect_size_mean', 'effect_size_std']

    # Take absolute values
    prolif_agg['effect_size_mean'] = prolif_agg['effect_size_mean'].abs()

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(prolif_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = prolif_agg[prolif_agg['platform'] == platform].sort_values('depth')
        means = subset['effect_size_mean'].values
        stds = subset['effect_size_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Mean Absolute Effect Size (Rank-Biserial)', fontsize=13, fontweight='bold')
    ax.set_title('Proliferation Signal vs Depth - Absolute Effect Size (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    plt.savefig('fig5d_effect_size_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig5d_effect_size_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig5d_effect_size_mean_sd.{pdf,png}")

    # --- Figure 5e: KNN Purity (Barplot) ---
    print("\nGenerating KNN purity (genome) barplot...")

    genome_results_no_full = genome_results.copy()

    # Exclude full depth samples
    if 'is_full_depth' in genome_results_no_full.columns:
        genome_results_no_full = genome_results_no_full[genome_results_no_full['is_full_depth'] == False]

    knn_agg = genome_results_no_full.groupby(['platform', 'depth']).agg({
        'knn_purity': ['mean', 'std']
    }).reset_index()
    knn_agg.columns = ['platform', 'depth', 'knn_mean', 'knn_std']

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(knn_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = knn_agg[knn_agg['platform'] == platform].sort_values('depth')
        means = subset['knn_mean'].values
        stds = subset['knn_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Mean KNN Purity', fontsize=13, fontweight='bold')
    ax.set_title('KNN Purity vs Depth (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    plt.savefig('fig5e_knn_purity_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig5e_knn_purity_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig5e_knn_purity_mean_sd.{pdf,png}")

    # --- Figure 6a: Genes Detected vs Depth (Barplot) ---
    print("\nGenerating genes detected barplot...")

    genome_results_no_full_6a = genome_results.copy()

    # Exclude full depth samples
    if 'is_full_depth' in genome_results_no_full_6a.columns:
        genome_results_no_full_6a = genome_results_no_full_6a[genome_results_no_full_6a['is_full_depth'] == False]

    genes_agg = genome_results_no_full_6a.groupby(['platform', 'depth']).agg({
        'n_genes_detected': ['mean', 'std']
    }).reset_index()
    genes_agg.columns = ['platform', 'depth', 'n_genes_mean', 'n_genes_std']

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(genes_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = genes_agg[genes_agg['platform'] == platform].sort_values('depth')
        means = subset['n_genes_mean'].values
        stds = subset['n_genes_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of Genes Detected', fontsize=13, fontweight='bold')
    ax.set_title('Genes Detected vs Depth (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    plt.savefig('fig6a_genes_detected_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig6a_genes_detected_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig6a_genes_detected_mean_sd.{pdf,png}")

    # --- Figure 6a2: Isoforms Detected vs Depth (Barplot) ---
    print("\nGenerating isoforms detected barplot...")

    isoforms_agg = genome_results_no_full_6a.groupby(['platform', 'depth']).agg({
        'n_isoforms_detected': ['mean', 'std']
    }).reset_index()
    isoforms_agg.columns = ['platform', 'depth', 'n_isoforms_mean', 'n_isoforms_std']

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(isoforms_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = isoforms_agg[isoforms_agg['platform'] == platform].sort_values('depth')
        means = subset['n_isoforms_mean'].values
        stds = subset['n_isoforms_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of Isoforms Detected', fontsize=13, fontweight='bold')
    ax.set_title('Isoforms Detected vs Depth (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    plt.savefig('fig6a2_isoforms_detected_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig6a2_isoforms_detected_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig6a2_isoforms_detected_mean_sd.{pdf,png}")

    # --- Figure 6b: Median UMI Counts per Cell vs Depth (Barplot) ---
    print("\nGenerating median UMI counts per cell barplot...")

    genome_results_no_full_6b = genome_results.copy()

    # Exclude full depth samples
    if 'is_full_depth' in genome_results_no_full_6b.columns:
        genome_results_no_full_6b = genome_results_no_full_6b[genome_results_no_full_6b['is_full_depth'] == False]

    umis_agg = genome_results_no_full_6b.groupby(['platform', 'depth']).agg({
        'median_umis': ['mean', 'std']
    }).reset_index()
    umis_agg.columns = ['platform', 'depth', 'median_umis_mean', 'median_umis_std']

    fig, ax = plt.subplots(figsize=(13, 6))

    depths_sorted = sorted(umis_agg['depth'].unique())
    n_depths = len(depths_sorted)
    n_platforms = len(platforms)
    bar_width = 0.2
    x_offset = np.arange(n_depths)

    for platform_idx, platform in enumerate(platforms):
        subset = umis_agg[umis_agg['platform'] == platform].sort_values('depth')
        means = subset['median_umis_mean'].values
        stds = subset['median_umis_std'].values

        x_pos = x_offset[:len(means)] + (platform_idx - n_platforms/2 + 0.5) * bar_width

        ax.bar(
            x_pos, means,
            width=bar_width,
            label=platform,
            color=platform_colors.get(platform, '#cccccc'),
            alpha=0.8,
            edgecolor='black',
            linewidth=1.5
        )
        ax.errorbar(
            x_pos, means,
            yerr=stds,
            fmt='none',
            ecolor='black',
            capsize=5,
            linewidth=1.5,
            alpha=0.7
        )

    ax.set_xlabel('Sequencing Depth (UMIs/cell)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Median UMI Counts per Cell', fontsize=13, fontweight='bold')
    ax.set_title('Median Transcripts per Cell vs Depth (Mean ± SD)', fontsize=15, fontweight='bold')
    ax.set_xticks(x_offset)
    ax.set_xticklabels([f'{int(d)}' for d in depths_sorted], rotation=45)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    plt.tight_layout()
    plt.savefig('fig6b_median_umis_mean_sd.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('fig6b_median_umis_mean_sd.png', dpi=150, bbox_inches='tight')
    print("✓ Saved: fig6b_median_umis_mean_sd.{pdf,png}")

    # --- Figure 6+: Genome Analysis UMAPs per Platform ---
    print("\nGenerating genome UMAP grids by platform...")

    depths_to_plot = sorted(set(target_depths))

    for platform in platforms:
        print(f"  {platform}...")

        fig = generate_genome_umaps_from_saved(
            platform_name=platform,
            target_depths=depths_to_plot,
            n_seeds=3,
            n_matrix_seeds=3,
            umap_data_dir='.umap_temp',
            cleanup=False,  # Set to True to delete saved files after generating figures
            include_full_depth=True  # Automatically includes full depth if available
        )

        # Save figure
        safe_name = platform.replace('-', '_').replace('+', 'plus').lower()
        fig.savefig(f'fig6_genome_umaps_{safe_name}.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(f'fig6_genome_umaps_{safe_name}.png', dpi=150, bbox_inches='tight')
        print(f"  ✓ Saved: fig6_genome_umaps_{safe_name}.{{pdf,png}}")
        plt.close(fig)

    # Only print additional summaries if we just ran benchmarking (not loaded from CSV)
    if not results_exist:
        # Create a metadata summary combining original depths and targets
        metadata_summary = []
        for platform in sorted(original_medians.keys()):
            orig_depth = original_medians[platform]
            for target_depth in sorted(set(target_depths)):
                downsampling_pct = (target_depth / orig_depth) * 100
                metadata_summary.append({
                    'platform': platform,
                    'original_median_umis': orig_depth,
                    'target_depth_umis': target_depth,
                    'downsampling_percent': downsampling_pct
                })

        metadata_df = pd.DataFrame(metadata_summary)
        metadata_df.to_csv("downsampling_metadata.csv", index=False)

        print(f"\n{'='*70}")
        print("Downsampling Metadata (Original vs Target)")
        print(f"{'='*70}")
        print(metadata_df.to_string(index=False))

    print(f"\n{'='*70}")
    print("✓ Analysis complete! Figures saved as PNG and PDF.")
    print("Note: .umap_temp directory retained for future plot regeneration.")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
