#!/usr/bin/env python3
"""
sqanti_metrics_multi.py

Extract key SQANTI3 structural-category metrics from one or more
classification.txt files, broken down by sample AND organism.

Organism is inferred from the 'chrom' column.
Handles mixed separators automatically:
  - triple underscore:  dog___1   → dog
  - single underscore:  human_1   → human, mouse_MT → mouse

Output: one row per sample + organism combination.

Examples
"""

__version__ = "1.0.0"
--------
# Infer sample names from file stems:
python sqanti_metrics_multi.py a/classification.txt b/classification.txt

# Provide explicit sample names:
python sqanti_metrics_multi.py a.txt b.txt --samples sampleA sampleB

# Write TSV:
python sqanti_metrics_multi.py a.txt b.txt --out metrics.tsv

# Filter by minimum isoform expression:
python sqanti_metrics_multi.py a.txt b.txt --min-iso-exp 1.0 --out metrics.tsv
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import sys
import pandas as pd


STRUCTURAL_LABELS = {
    "FSM": "full-splice_match",
    "ISM": "incomplete-splice_match",
    "NIC": "novel_in_catalog",
    "NNC": "novel_not_in_catalog",
}

CHROM_SEP = "___"  # fallback separator (overridable via --chrom-sep)


def pct(x: int, d: int) -> float:
    return 100.0 * x / d if d and d > 0 else float("nan")


def read_classification(path: Path) -> pd.DataFrame:
    """
    Read a SQANTI3 classification.txt file robustly.

    SQANTI3 output often has an implicit row-index as the first column
    (the isoform ID), with the header starting one field short — this
    causes pandas to misalign all columns by one, putting numeric data
    into the structural_category column.

    We detect this by comparing header field count vs data field count
    and set index_col=0 when needed.
    """
    # Peek at first two lines to detect column shift
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as fh:
        header_fields = len(fh.readline().rstrip("\n").split("\t"))
        data_fields   = len(fh.readline().rstrip("\n").split("\t"))

    # If data has one more field than header → unnamed index column present
    index_col = 0 if data_fields > header_fields else None

    try:
        df = pd.read_csv(path, sep="\t", low_memory=False, index_col=index_col)
    except Exception:
        df = pd.read_csv(path, sep=None, engine="python", low_memory=False,
                         index_col=index_col)

    # Strip accidental whitespace from column names
    df.columns = df.columns.str.strip()
    return df


def require_cols(df: pd.DataFrame, path: Path, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: Missing required columns in {path}:\n"
            f"  Missing: {missing}\n"
            f"  Available: {', '.join(df.columns)}"
        )


def extract_organism(chrom_series: pd.Series) -> pd.Series:
    """
    Robustly extract organism prefix from chrom values.

    Handles mixed separator styles in the same file:
      dog___1   -> dog      (triple underscore)
      human_1   -> human    (single underscore)
      mouse_MT  -> mouse    (single underscore)

    Strategy:
      1. If triple underscore '___' is present → split on '___', take [0].
      2. Otherwise auto-detect: collect all unique single-underscore prefixes
         across the series, pick the longest matching prefix per value.
         This correctly handles chromosomes like 'mouse_1' vs 'mouse_MT'.
      3. Any chrom with no underscore at all → 'unknown'.
    """
    chrom = chrom_series.astype(str)

    # Split on triple underscore where present
    has_triple = chrom.str.contains("___", na=False)

    # For single-underscore chroms, auto-detect all prefixes present
    single_uc = chrom[~has_triple]
    # Collect candidate prefixes (everything before the last underscore segment
    # that looks like a chrom number/name: digit, MT, X, Y)
    import re as _re
    prefix_pattern = _re.compile(r'^([A-Za-z][A-Za-z0-9]*)_+(?:\d+|MT|X|Y)$')
    detected_prefixes = set()
    for val in single_uc.unique():
        m = prefix_pattern.match(val)
        if m:
            detected_prefixes.add(m.group(1))

    def _get_organism(val: str) -> str:
        if "___" in val:
            return val.split("___")[0]
        # Try matching against detected prefixes (longest match wins)
        for prefix in sorted(detected_prefixes, key=len, reverse=True):
            if val.startswith(prefix + "_"):
                return prefix
        # Generic fallback: take everything before the first underscore
        if "_" in val:
            return val.split("_")[0]
        return "unknown"

    return chrom.map(_get_organism)


def compute_metrics_for_group(df: pd.DataFrame) -> dict:
    """Compute all metrics for a single (sample, organism) subset."""
    total_isoforms = len(df)
    mono_exonic    = int((df["exons"] == 1).sum())
    multi_exonic   = int((df["exons"] >= 2).sum())

    # Convert to plain dict — Series.get() can misbehave when the index
    # contains substrings of the lookup key, returning a slice instead of 0.
    sc_counts = df["structural_category"].value_counts().to_dict()

    fsm = int(sc_counts.get(STRUCTURAL_LABELS["FSM"], 0))
    ism = int(sc_counts.get(STRUCTURAL_LABELS["ISM"], 0))
    nic = int(sc_counts.get(STRUCTURAL_LABELS["NIC"], 0))
    nnc = int(sc_counts.get(STRUCTURAL_LABELS["NNC"], 0))
    splice_total = fsm + ism + nic + nnc

    unique_genes = (
        int(df["associated_gene"].dropna().nunique())
        if "associated_gene" in df.columns else None
    )

    canon_true = None
    canon_pct  = None
    if "all_canonical" in df.columns:
        canon      = df["all_canonical"].astype(str).str.upper()
        canon_true = int((canon == "TRUE").sum() + (canon == "1").sum())
        canon_pct  = pct(canon_true, total_isoforms)

    out = {
        "n_isoforms":       total_isoforms,
        "n_genes":          unique_genes,
        "n_mono_exonic":    mono_exonic,
        "n_multi_exonic":   multi_exonic,
        "FSM":              fsm,
        "ISM":              ism,
        "NIC":              nic,
        "NNC":              nnc,
        "splice_total":     splice_total,
        "pct_FSM_of_splice": pct(fsm, splice_total),
        "pct_ISM_of_splice": pct(ism, splice_total),
        "pct_NIC_of_splice": pct(nic, splice_total),
        "pct_NNC_of_splice": pct(nnc, splice_total),
        "pct_FSM_of_all":   pct(fsm, total_isoforms),
        "pct_ISM_of_all":   pct(ism, total_isoforms),
        "pct_NIC_of_all":   pct(nic, total_isoforms),
        "pct_NNC_of_all":   pct(nnc, total_isoforms),
    }

    if canon_true is not None:
        out["n_all_canonical"]   = canon_true
        out["pct_all_canonical"] = canon_pct

    return out


def infer_sample_name(path: Path) -> str:
    return path.stem


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Extract SQANTI3 structural-category metrics from one or more "
            "classification.txt files, broken down by sample and organism."
        )
    )
    ap.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="One or more SQANTI3 classification.txt files",
    )
    ap.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help=(
            "Optional sample names matching the order of input files. "
            "If omitted, inferred from filename stems."
        ),
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output TSV path (table is also printed to stdout).",
    )
    ap.add_argument(
        "--min-iso-exp",
        type=float,
        default=None,
        help=(
            "Optional: keep only rows with iso_exp >= this value before "
            "computing metrics (requires 'iso_exp' column)."
        ),
    )
    ap.add_argument(
        "--chrom-sep",
        default=CHROM_SEP,
        help=(
            f"Separator between organism prefix and chromosome name in the "
            f"'chrom' column (default: '{CHROM_SEP}')."
        ),
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    if args.samples is not None and len(args.samples) != len(args.files):
        raise SystemExit(
            f"ERROR: --samples provided {len(args.samples)} names but "
            f"{len(args.files)} input files were given."
        )

    rows = []

    for i, path in enumerate(args.files):
        if not path.exists():
            raise SystemExit(f"ERROR: File not found: {path}")

        sample = (
            args.samples[i]
            if args.samples is not None
            else infer_sample_name(path)
        )
        print(f"Processing: {path}  (sample='{sample}')", flush=True)

        df = read_classification(path)
        require_cols(df, path, ["structural_category", "exons", "chrom"])

        # Optional expression filter
        if args.min_iso_exp is not None:
            if "iso_exp" not in df.columns:
                raise SystemExit(
                    f"ERROR: --min-iso-exp used but 'iso_exp' column not "
                    f"found in {path}"
                )
            before = len(df)
            df = df[df["iso_exp"] >= args.min_iso_exp].copy()
            print(
                f"  iso_exp filter >= {args.min_iso_exp}: "
                f"{before} → {len(df)} isoforms retained"
            )

        # Extract organism from chrom column (handles mixed ___ and _ separators)
        df = df.copy()
        df["_organism"] = extract_organism(df["chrom"])

        organisms = sorted(df["_organism"].unique())
        print(f"  Organisms detected: {organisms}")
        print(f"  Structural categories found: {sorted(df['structural_category'].dropna().unique())}")

        for organism in organisms:
            subset = df[df["_organism"] == organism]
            metrics = compute_metrics_for_group(subset)
            rows.append({
                "sample":   sample,
                "organism": organism,
                "file":     str(path),
                **metrics,
            })

    out_df = pd.DataFrame(rows)

    # Column order: sample, organism first, then metrics
    priority = ["sample", "organism", "file"]
    cols = priority + [c for c in out_df.columns if c not in priority]
    out_df = out_df[cols]

    # Print
    print("\n=== SQANTI3 metrics — per sample × organism ===")
    if args.min_iso_exp is not None:
        print(f"Filter applied: iso_exp >= {args.min_iso_exp}")
    print()
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(out_df.to_string(index=False))

    if args.out is not None:
        out_df.to_csv(args.out, sep="\t", index=False)
        print(f"\nWrote TSV: {args.out}")


if __name__ == "__main__":
    main()
