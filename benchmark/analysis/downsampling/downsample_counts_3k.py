"""
Downsample all h5ad matrices in INPUT_DIR at counts level using sc.pp.downsample_counts.

Scans INPUT_DIR for every *.h5ad file and applies count downsampling for each
combination of target depth × downsample seed, preserving the original filename
as the output prefix so cluster/resolution/platform labels are carried through.

Output directory: INPUT_DIR/downsampled_counts/
Output naming:    {original_stem}_{depth}_s{downsample_seed}.h5ad
"""

import os
import glob
import scanpy as sc
import numpy as np

# ──────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR  = BASE_DIR
OUTPUT_DIR = os.path.join(BASE_DIR, "downsampled_counts")

DOWNSAMPLE_SEEDS = [0, 1, 2]
TARGET_DEPTHS    = [3000]

os.makedirs(OUTPUT_DIR, exist_ok=True)


def downsample_adata(adata, target_umis, seed):
    adata_ds = adata.copy()
    sc.pp.downsample_counts(adata_ds, counts_per_cell=target_umis, random_state=seed)
    adata_ds.obs["n_counts"] = np.asarray(adata_ds.X.sum(axis=1)).ravel()
    return adata_ds


def main():
    input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))

    if not input_files:
        print(f"No .h5ad files found in {INPUT_DIR}")
        return

    total = len(input_files) * len(TARGET_DEPTHS) * len(DOWNSAMPLE_SEEDS)
    done  = 0

    for input_path in input_files:
        stem = os.path.splitext(os.path.basename(input_path))[0]
        print(f"\nLoading {stem}...", flush=True)
        adata = sc.read_h5ad(input_path)
        print(f"  {adata.n_obs} cells", flush=True)

        for depth in TARGET_DEPTHS:
            for ds_seed in DOWNSAMPLE_SEEDS:
                out_name = f"{stem}_{depth}_s{ds_seed}.h5ad"
                out_path = os.path.join(OUTPUT_DIR, out_name)

                if os.path.exists(out_path):
                    print(f"  [SKIP] already exists: {out_name}", flush=True)
                    done += 1
                    continue

                print(
                    f"  depth={depth}  seed={ds_seed}  →  {out_name}",
                    end="  ", flush=True
                )
                adata_ds = downsample_adata(adata, target_umis=depth, seed=ds_seed)
                adata_ds.write_h5ad(out_path)
                done += 1
                print(f"✓  ({done}/{total})", flush=True)

    print(f"\nDone. {done}/{total} files written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
