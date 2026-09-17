import scanpy as sc
import argparse
from pathlib import Path


def downsample_h5ad(input_file, n_cells, n_replicates=3, cluster_key=None, clusters=None, dry_run=False):
    """
    Downsample an h5ad file to a specific number of cells with multiple replicates.

    Parameters
    ----------
    input_file : str
        Path to input h5ad file
    n_cells : int
        Number of cells to downsample to
    n_replicates : int
        Number of replicates to create (default: 3)
    cluster_key : str or None
        obs column name containing cluster labels (e.g. 'leiden_0.5').
        Required when `clusters` is provided.
    clusters : list of str or None
        Cluster labels to restrict downsampling to. Requires `cluster_key`.
    dry_run : bool
        If True, only report cell counts without writing any files.
    """
    print(f"Reading {input_file}...")
    adata = sc.read_h5ad(input_file)
    print(f"Total cells in dataset: {adata.n_obs}")

    if clusters is not None:
        if cluster_key is None:
            raise ValueError("--cluster-key is required when --clusters is provided.")
        if cluster_key not in adata.obs.columns:
            raise KeyError(
                f"Column '{cluster_key}' not found in obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        obs_labels = adata.obs[cluster_key].astype(str)
        mask = obs_labels.isin([str(c) for c in clusters])
        adata_subset = adata[mask]

        print(f"Cluster key : {cluster_key}")
        print(f"Cluster(s)  : {clusters}")
        print(f"Cells in selected cluster(s): {adata_subset.n_obs}")

        if adata_subset.n_obs == 0:
            raise ValueError(
                f"No cells found for cluster(s) {clusters} in column '{cluster_key}'. "
                f"Available values: {sorted(obs_labels.unique().tolist())}"
            )

        if dry_run:
            return

        adata = adata_subset.copy()
    elif dry_run:
        # dry-run without cluster filtering: just report total count
        return

    if n_cells > adata.n_obs:
        raise ValueError(
            f"Cannot downsample to {n_cells} cells. "
            f"{'Selected subset' if clusters else 'Input'} has only {adata.n_obs} cells."
        )

    input_path = Path(input_file)
    sample_name = input_path.stem

    cluster_suffix = ""
    if clusters is not None:
        cluster_tag = "_".join(str(c) for c in clusters)
        cluster_suffix = f"_{cluster_key}_clusters{cluster_tag}"

    for seed in range(n_replicates):
        print(f"\nCreating replicate {seed + 1}/{n_replicates} (seed={seed})...")

        adata_sampled = adata.copy()
        sc.pp.subsample(adata_sampled, n_obs=n_cells, random_state=seed)

        output_file = f"{sample_name}_{n_cells}_cells{cluster_suffix}_seed{seed}.h5ad"
        adata_sampled.write_h5ad(output_file)
        print(f"Saved: {output_file} ({adata_sampled.n_obs} cells)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Downsample an h5ad file to a fixed number of cells."
    )
    parser.add_argument("input_file", help="Path to input h5ad file")
    parser.add_argument("n_cells", type=int, help="Number of cells to downsample to")
    parser.add_argument(
        "n_replicates", type=int, nargs="?", default=3,
        help="Number of replicates (default: 3)"
    )
    parser.add_argument(
        "--cluster-key",
        help="obs column name with cluster labels (e.g. 'leiden_0.5')"
    )
    parser.add_argument(
        "--clusters",
        help="Comma-separated cluster labels to restrict downsampling to (requires --cluster-key)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report cell counts only; do not write any output files"
    )

    args = parser.parse_args()

    clusters = [c.strip() for c in args.clusters.split(",")] if args.clusters else None

    downsample_h5ad(
        args.input_file,
        args.n_cells,
        args.n_replicates,
        cluster_key=args.cluster_key,
        clusters=clusters,
        dry_run=args.dry_run,
    )
