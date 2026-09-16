#!/usr/bin/env python3
"""
Metagene coverage pipeline using deepTools.
Generates a 2x4 figure:
  - Row 1: metagene profiles for top 10% SHORTEST genes, one panel per category
  - Row 2: metagene profiles for top 10% LONGEST genes, one panel per category

Dependencies:
    deeptools, pandas, numpy, matplotlib

Install:
    pip install deeptools pandas numpy matplotlib
"""

__version__ = "1.0.0"

import os
import subprocess
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gzip
import pyBigWig

# ---------------------------------------------------------------------------
# 0. USER CONFIGURATION
# ---------------------------------------------------------------------------

GTF_FILE   = "data/genes.gtf.gz"
OUTPUT_DIR = "metagene_output"

# BAM files: dict of {category_name: [rep1.bam, rep2.bam, rep3.bam]}
# Replace with your own paths (relative or absolute).
BAMS = {
    "10X-3PRIME": [
        "data/bam/10X-3PRIME_REP1.genome.dedup.bam",
        "data/bam/10X-3PRIME_REP2.genome.dedup.bam",
        "data/bam/10X-3PRIME_REP3.genome.dedup.bam",
    ],
    "10X-5PRIME": [
        "data/bam/10X_5PRIME_REP1.genome.dedup.bam",
        "data/bam/10X_5PRIME_REP2.genome.dedup.bam",
        "data/bam/10X_5PRIME_REP3.genome.dedup.bam",
    ],
    "ARGENTAG": [
        "data/bam/ARGENTAG_REP1.genome.dedup.bam",
        "data/bam/ARGENTAG_REP2.genome.dedup.bam",
        "data/bam/ARGENTAG_REP3.genome.dedup.bam",
    ],
    "PARSE": [
        "data/bam/PARSE_REP1.genome.dedup.bam",
        "data/bam/PARSE_REP2.genome.dedup.bam",
        "data/bam/PARSE_REP3.genome.dedup.bam",
    ],
}

# ── COLOR CONFIG ────────────────────────────────────────────────────────────
CATEGORY_COLORS = {
    "10X-3PRIME": "#658b66",
    "10X-5PRIME": "#8aa783",
    "ARGENTAG":   "#2680a6",
    "PARSE":      "#875670",
}
# ────────────────────────────────────────────────────────────────────────────

# deepTools settings
THREADS   = 8
BINSIZE   = 10      # bamCoverage bin size in bp
NORMALIZE = "CPM"   # CPM, RPKM, BPM, RPGC, or None

# computeMatrix settings
BODY_LENGTH = 5000  # scaled gene body length in bp
UPSTREAM    = 50  # bp upstream of TSS
DOWNSTREAM  = 50  # bp downstream of TES

# Plot aesthetics
FIGURE_WIDTH  = 20
FIGURE_HEIGHT = 8
LINE_WIDTH    = 2.0
FONT_SIZE     = 10
TITLE_SIZE    = 11

# ---------------------------------------------------------------------------
# 1. PARSE GTF — exon union length for classification + gene span for BED
# ---------------------------------------------------------------------------

def compute_exon_union_lengths(gtf_file: str):
    """
    Parse GTF and compute exon-union length per gene.

    WHY TWO THINGS:
      - gene_df contains gene_start/gene_end (= full genomic span from first
        to last exon) used for the BED file so computeMatrix sees the whole
        transcription unit and captures 5'->3' gradients
      - union_length is the sum of merged exon lengths used ONLY for
        classifying genes as short/long — this is based on actual coding
        content, not intron-inflated span which would be misleading

    Returns:
        gene_df  — one row per gene: gene_start, gene_end, union_length
        exon_df  — raw exon records (used only for union_length computation)
    """
    print("[1/5] Parsing GTF and computing exon union lengths...")
    rows = []

    opener = gzip.open if gtf_file.endswith(".gz") else open
    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            gene_id = None
            transcript_id = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1]
                    
                elif attr.startswith("transcript_id"):
                    transcript_id = attr.split('"')[1]
                    
            if gene_id is None or transcript_id is None:
                continue
            rows.append({
                "gene_id":       gene_id,
                "transcript_id": transcript_id,  # ← add this
                "chrom":         chrom,
                "strand":        strand,
                "start":         int(start) - 1,
                "end":           int(end),
            })

    if not rows:
        raise ValueError("No exon records found in GTF. Check file format.")

    exon_df = pd.DataFrame(rows)

    results = []
    for (gene_id, chrom, strand), grp in exon_df.groupby(
            ["gene_id", "chrom", "strand"]):
        intervals = grp[["start", "end"]].sort_values("start").values.tolist()
        # Merge overlapping exons to avoid double-counting
        merged = [list(intervals[0])]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        union_length = sum(e - s for s, e in merged)
        results.append({
            "gene_id":      gene_id,
            "chrom":        chrom,
            "strand":       strand,
            "gene_start":   merged[0][0],   # full span start (for BED)
            "gene_end":     merged[-1][1],  # full span end   (for BED)
            "union_length": union_length,   # exon-only length (for classification)
        })

    gene_df = pd.DataFrame(results)
    print(f"    Found {len(gene_df)} genes with exon annotations.")
    return gene_df, exon_df


# ---------------------------------------------------------------------------
# 2. SELECT TOP/BOTTOM 10% BY EXON UNION LENGTH, WRITE GENE-SPAN BED FILES
# ---------------------------------------------------------------------------

def select_longest_transcripts(exon_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each gene, retain only the exons belonging to the transcript
    with the greatest total exonic length (longest transcript).

    This avoids the isoform-mixing problem where merged exon intervals
    from multiple transcripts create composite regions that don't
    correspond to any real transcript.
    """
    print("    Selecting longest transcript per gene...")

    transcript_lengths = []
    for (gene_id, transcript_id), grp in exon_df.groupby(["gene_id", "transcript_id"]):
        length = (grp["end"] - grp["start"]).sum()
        transcript_lengths.append({
            "gene_id":       gene_id,
            "transcript_id": transcript_id,
            "tx_length":     length,
        })

    tx_df = pd.DataFrame(transcript_lengths)
    # Pick the longest transcript per gene; ties broken by transcript_id (alphabetical)
    longest = (tx_df.sort_values("tx_length", ascending=False)
                    .groupby("gene_id", sort=False)
                    .first()
                    .reset_index()[["gene_id", "transcript_id"]])

    exon_df_filtered = exon_df.merge(longest, on=["gene_id", "transcript_id"])
    print(f"    Retained {exon_df_filtered['gene_id'].nunique()} genes "
          f"({exon_df_filtered['transcript_id'].nunique()} transcripts)")
    return exon_df_filtered

def filter_genes_by_coverage(gene_df: pd.DataFrame, exon_df: pd.DataFrame,
                              avg_bws: dict,
                              min_mean_coverage: float = 10.0) -> pd.DataFrame:
    """
    Retain only genes with mean BigWig coverage >= min_mean_coverage
    over EXONIC regions in ALL categories.

    Exon intervals are used instead of gene span to avoid penalising
    genes with long introns — intronic coverage is expected to be low
    or zero for spliced-RNA technologies (10X, PARSE, ARGENTAG).
    """

    print(f"[2b] Filtering genes by exonic coverage (min mean = {min_mean_coverage})...")
    before = len(gene_df)

    # Pre-compute merged exon intervals per gene for fast lookup
    exon_intervals = {}
    for (gene_id, chrom, strand), grp in exon_df.groupby(["gene_id", "chrom", "strand"]):
        intervals = grp[["start", "end"]].sort_values("start").values.tolist()
        merged = [list(intervals[0])]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        exon_intervals[gene_id] = (chrom, merged)

    keep = np.ones(len(gene_df), dtype=bool)

    for category, bw_path in avg_bws.items():
        bw = pyBigWig.open(bw_path)
        means = []
        for _, row in gene_df.iterrows():
            gene_id = row["gene_id"]
            if gene_id not in exon_intervals:
                means.append(0.0)
                continue
            chrom, intervals = exon_intervals[gene_id]
            # Average coverage across all merged exon intervals
            vals = []
            for s, e in intervals:
                try:
                    v = bw.stats(chrom, s, e, type="mean")[0]
                    if v is not None:
                        vals.append(v)
                except RuntimeError:
                    pass  # contig not in BigWig
            means.append(np.mean(vals) if vals else 0.0)
        bw.close()

        category_pass = np.array(means) >= min_mean_coverage
        keep &= category_pass
        print(f"    {category}: {category_pass.sum()} / {len(gene_df)} genes pass")

    gene_df_filtered = gene_df[keep].copy()
    print(f"    Kept {len(gene_df_filtered)} / {before} genes passing all categories")
    return gene_df_filtered


def write_bed_files(gene_df: pd.DataFrame, exon_df: pd.DataFrame,
                    output_dir: str, min_body_bp: int = 250,
                    region_mode: str = "exon"):
    """
    region_mode = "exon"      → one BED entry per merged exon interval
                                (longest transcript only, no introns)
    region_mode = "gene_body" → one BED entry per gene (TSS→TES, includes introns)
    """
    print(f"[2/5] Writing BED files (mode: {region_mode})...")

    before = len(gene_df)
    gene_df = gene_df[(gene_df["gene_end"] - gene_df["gene_start"]) >= min_body_bp].copy()
    print(f"    Excluded {before - len(gene_df)} genes with body < {min_body_bp}bp "
          f"({len(gene_df)} remaining)")

    p10 = np.percentile(gene_df["union_length"], 10)
    p90 = np.percentile(gene_df["union_length"], 90)
    short_ids = set(gene_df[gene_df["union_length"] <= p10]["gene_id"])
    long_ids  = set(gene_df[gene_df["union_length"] >= p90]["gene_id"])

    print(f"    Short genes (<= {p10:.0f} bp exon union): {len(short_ids)}")
    print(f"    Long  genes (>= {p90:.0f} bp exon union): {len(long_ids)}")

    os.makedirs(output_dir, exist_ok=True)

    def make_bed(gene_ids, path):
        rows = []
        if region_mode == "gene_body":
            # One entry per gene: full TSS→TES span
            subset = gene_df[gene_df["gene_id"].isin(gene_ids)]
            for _, row in subset.iterrows():
                rows.append({
                    "chrom":  row["chrom"],
                    "start":  row["gene_start"],
                    "end":    row["gene_end"],
                    "name":   row["gene_id"],
                    "score":  0,
                    "strand": row["strand"],
                })
        else:
            # One entry per merged exon interval (longest transcript)
            subset = exon_df[exon_df["gene_id"].isin(gene_ids)]
            for (gene_id, chrom, strand), grp in subset.groupby(
                    ["gene_id", "chrom", "strand"]):
                intervals = grp[["start", "end"]].sort_values("start").values.tolist()
                merged = [list(intervals[0])]
                for s, e in intervals[1:]:
                    if s <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], e)
                    else:
                        merged.append([s, e])
                for s, e in merged:
                    rows.append({
                        "chrom":  chrom,
                        "start":  s,
                        "end":    e,
                        "name":   gene_id,
                        "score":  0,
                        "strand": strand,
                    })

        bed = pd.DataFrame(rows)
        bed = bed[bed["end"] > bed["start"]]
        bed.to_csv(path, sep="\t", header=False, index=False)
        print(f"    Written {len(bed)} intervals -> {path}")

    suffix = "body" if region_mode == "gene_body" else "exons"
    short_bed = os.path.join(output_dir, f"short_genes_{suffix}.bed")
    long_bed  = os.path.join(output_dir, f"long_genes_{suffix}.bed")
    make_bed(short_ids, short_bed)
    make_bed(long_ids,  long_bed)

    return short_bed, long_bed

# ---------------------------------------------------------------------------
# 3. bamCoverage -> BigWig, then bigwigAverage per category
# ---------------------------------------------------------------------------

def run_bamCoverage(bams: dict, output_dir: str,
                    binsize: int, normalize: str, threads: int) -> dict:
    print("[3/5] Running bamCoverage and averaging replicates per category...")
    bw_dir = os.path.join(output_dir, "bigwigs")
    os.makedirs(bw_dir, exist_ok=True)

    avg_bws = {}
    for category, bam_list in bams.items():
        rep_bws = []
        for i, bam in enumerate(bam_list, 1):
            bw_out = os.path.join(bw_dir, f"{category}_rep{i}.bw")
            cmd = [
                "bamCoverage",
                "-b", bam,
                "-o", bw_out,
                "--binSize",            str(binsize),
                "--normalizeUsing",     normalize,
                "--numberOfProcessors", str(threads),
                "--skipNonCoveredRegions",
            ]
            print(f"    {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            rep_bws.append(bw_out)

        avg_bw = os.path.join(bw_dir, f"{category}_avg.bw")
        if len(rep_bws) == 1:
            import shutil
            shutil.copy(rep_bws[0], avg_bw)
            print(f"    Single replicate — copied to {avg_bw}")
        else:
            cmd = [
                "bigwigAverage",
                "-b", *rep_bws,
                "-o", avg_bw,
                "--numberOfProcessors", str(threads),
            ]
            print(f"    Averaging {len(rep_bws)} replicates -> {avg_bw}")
            subprocess.run(cmd, check=True)
        avg_bws[category] = avg_bw

    return avg_bws


def get_existing_bigwigs(output_dir: str) -> dict:
    """Return paths to already-computed average BigWigs (for --skip-bam)."""
    bw_dir = os.path.join(output_dir, "bigwigs")
    avg_bws = {}
    for category in BAMS.keys():
        bw = os.path.join(bw_dir, f"{category}_avg.bw")
        if not os.path.exists(bw):
            raise FileNotFoundError(
                f"Expected BigWig not found: {bw}\n"
                f"Run without --skip-bam first, or check your output directory."
            )
        avg_bws[category] = bw
        print(f"    Found existing BigWig: {bw}")
    return avg_bws


# ---------------------------------------------------------------------------
# 4. computeMatrix + plotProfile --outFileNameData per category x gene set
# ---------------------------------------------------------------------------

def run_compute_matrix(avg_bws: dict, bed_file: str, label: str,
                       output_dir: str, body_length: int,
                       upstream: int, downstream: int, threads: int) -> dict:
    """
    Run computeMatrix scale-regions on gene body BED intervals.

    NOTE: --skipZeros is intentionally OMITTED.
    For snRNA-seq / snONT data, introns genuinely have low/zero coverage.
    Skipping zeros would remove these regions and distort the 5'->3' gradient
    that we are specifically trying to visualise.

    Returns {category: profile_tsv_path}
    """
    matrix_dir = os.path.join(output_dir, f"matrices_{label}")
    os.makedirs(matrix_dir, exist_ok=True)

    profile_tables = {}
    for category, bw in avg_bws.items():
        matrix_file  = os.path.join(matrix_dir, f"{category}_{label}.gz")
        profile_data = os.path.join(matrix_dir, f"{category}_{label}_profile.tab")
        dummy_plot   = os.path.join(matrix_dir, f"{category}_{label}_profile.png")

        cmd = [
            "computeMatrix", "scale-regions",
            "-S", bw,
            "-R", bed_file,
            "--regionBodyLength",   str(body_length),
            "--upstream",           str(upstream),
            "--downstream",         str(downstream),
            # No --skipZeros: intron zeros are real signal for snONT
            "--numberOfProcessors", str(threads),
            "-o", matrix_file,
            "--samplesLabel", category,
        ]
        print(f"    computeMatrix [{label}] {category}")
        subprocess.run(cmd, check=True)

        cmd = [
            "plotProfile",
            "-m", matrix_file,
            "-o", dummy_plot,
            "--outFileNameData", profile_data,
            "--averageType", "mean",
        ]
        subprocess.run(cmd, check=True)
        profile_tables[category] = profile_data

    return profile_tables


def run_plotprofile_only(output_dir: str, categories: list) -> tuple:
    """Re-run only plotProfile on existing .gz matrices (for --from-matrices)."""
    short_tables = {}
    long_tables  = {}

    for label, tables in [("short", short_tables), ("long", long_tables)]:
        matrix_dir = os.path.join(output_dir, f"matrices_{label}")
        for category in categories:
            matrix_file  = os.path.join(matrix_dir, f"{category}_{label}.gz")
            profile_data = os.path.join(matrix_dir, f"{category}_{label}_profile.tab")
            dummy_plot   = os.path.join(matrix_dir, f"{category}_{label}_profile.png")

            if not os.path.exists(matrix_file):
                raise FileNotFoundError(
                    f"Matrix not found: {matrix_file}\n"
                    f"Run without --from-matrices first."
                )
            cmd = [
                "plotProfile",
                "-m", matrix_file,
                "-o", dummy_plot,
                "--outFileNameData", profile_data,
                "--averageType", "mean",
            ]
            print(f"    plotProfile [{label}] {category}")
            subprocess.run(cmd, check=True)
            tables[category] = profile_data

    return short_tables, long_tables


def get_existing_profile_tables(output_dir: str, categories: list) -> tuple:
    """Return paths to existing .tab files without re-running anything (for --from-plotprofile)."""
    short_tables = {}
    long_tables  = {}

    for label, tables in [("short", short_tables), ("long", long_tables)]:
        matrix_dir = os.path.join(output_dir, f"matrices_{label}")
        for category in categories:
            tab = os.path.join(matrix_dir, f"{category}_{label}_profile.tab")
            if not os.path.exists(tab):
                raise FileNotFoundError(
                    f"Profile table not found: {tab}\n"
                    f"Run with --from-matrices (not --from-plotprofile) first."
                )
            tables[category] = tab
            print(f"    Found existing profile: {tab}")

    return short_tables, long_tables


# ---------------------------------------------------------------------------
# PARSE plotProfile OUTPUT TABLE
# ---------------------------------------------------------------------------

def parse_profile_table(tsv_path: str):
    """
    Parse plotProfile --outFileNameData output robustly.

    deeptools format (2-line header):
      Line 0: "bin labels  <empty>  -upstream  ...  TSS  ...  TES  ...  +downstream"
      Line 1: "bins        1.0      2.0        ...  (numeric bin indices)"
      Line 2+: "<sample_name>  val  val  val  ..."

    Returns:
        x -- bin indices as float array
        y -- mean coverage across all data rows, non-numeric fields skipped
    """
    with open(tsv_path) as fh:
        lines = fh.readlines()

    # Row 1: bin indices — skip first field ("bins"), drop empty strings
    x_raw = lines[1].strip().split("\t")[1:]
    x = np.array([v for v in x_raw if v.strip() != ""], dtype=float)

    # Row 2+: data — skip any non-numeric fields (sample name in col 0)
    data_rows = []
    for line in lines[2:]:
        if not line.strip():
            continue
        numeric = []
        for v in line.strip().split("\t"):
            try:
                numeric.append(float(v))
            except ValueError:
                continue
        if numeric:
            data_rows.append(numeric)

    if not data_rows:
        raise ValueError(f"No numeric data rows found in {tsv_path}")

    y = np.mean(data_rows, axis=0)
    y = y[:len(x)]  # trim trailing empties if any

    return x, y


def minmax_normalize(y: np.ndarray) -> np.ndarray:
    """
    Scale y to [0, 1].

    WHY: Different sequencing technologies have very different absolute
    coverage depths. Min-max normalization makes profile SHAPES directly
    comparable across categories without being dominated by depth differences.
    The y-axis then represents relative coverage (0 = min, 1 = max of that
    profile), which is what matters for comparing 5'->3' gradients.
    """
    ymin, ymax = y.min(), y.max()
    if ymax - ymin < 1e-10:
        return np.zeros_like(y)
    return (y - ymin) / (ymax - ymin)


# ---------------------------------------------------------------------------
# 5. ASSEMBLE 2×N MATPLOTLIB FIGURE
# ---------------------------------------------------------------------------

def build_figure(short_tables: dict, long_tables: dict,
                 category_colors: dict, output_dir: str,
                 upstream: int, downstream: int, body_length: int):
    print("[5/5] Assembling figure...")

    categories  = list(short_tables.keys())
    n_cats      = len(categories)
    total_span  = upstream + body_length + downstream

    fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig.suptitle("Metagene gene-body coverage profiles (min-max normalised)",
                 fontsize=TITLE_SIZE + 2, y=1.02)

    gs = gridspec.GridSpec(2, n_cats, figure=fig, hspace=0.5, wspace=0.35)

    row_configs = [
        (short_tables, "Short genes (bottom 10% gene body)"),
        (long_tables,  "Long genes (top 10% gene body)"),
    ]

    for row_idx, (tables, row_label) in enumerate(row_configs):
        for col_idx, category in enumerate(categories):
            ax     = fig.add_subplot(gs[row_idx, col_idx])
            x, y   = parse_profile_table(tables[category])
            y_norm = minmax_normalize(y)
            color  = category_colors.get(category, "steelblue")

            ax.plot(x, y_norm, color=color, linewidth=LINE_WIDTH)
            ax.fill_between(x, y_norm, alpha=0.12, color=color)

            # TSS and TES bin positions
            n_bins  = len(x)
            tss_bin = round(n_bins * upstream / total_span)
            tes_bin = round(n_bins * (upstream + body_length) / total_span) - 1

            ax.axvline(x=x[tss_bin], color="gray", linestyle="--", linewidth=0.8)
            ax.axvline(x=x[tes_bin], color="gray", linestyle="--", linewidth=0.8)

            #tick_positions = [x[0], x[tss_bin], x[tes_bin], x[-1]]
            #tick_labels    = [f"-{upstream}bp", "TSS", "TES", f"+{downstream}bp"]
            tick_positions = [x[tss_bin], x[tes_bin]]
            tick_labels    = ["TSS", "TES"]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, fontsize=FONT_SIZE - 2,
                               rotation=15, ha="right")
            ax.set_ylim(-0.05, 1.15)
            ax.tick_params(labelsize=FONT_SIZE - 1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row_idx == 0:
                ax.set_title(category, fontsize=TITLE_SIZE,
                             fontweight="bold", color=color)

            if col_idx == 0:
                ax.set_ylabel(
                    f"{row_label}\nNormalised coverage (0-1)",
                    fontsize=FONT_SIZE - 1,
                )
            else:
                ax.set_ylabel("Normalised coverage (0-1)", fontsize=FONT_SIZE - 1)

            ax.set_xlabel("Genomic position", fontsize=FONT_SIZE - 1)

    os.makedirs(output_dir, exist_ok=True)
    for ext in ["png", "pdf", "svg"]:
        out_path = os.path.join(output_dir, f"metagene_genebody_2xN.{ext}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"    Saved: {out_path}")

    plt.close(fig)

def build_overlaid_figure(short_tables: dict, long_tables: dict,
                           category_colors: dict, output_dir: str,
                           upstream: int, downstream: int, body_length: int):
    """
    Alternative output: one panel per gene-length class, all category curves
    overlaid in each panel. Produces two subplots side by side:
      Left  — short genes (bottom 10%)
      Right — long genes  (top 10%)
    """
    print("[5b/5] Assembling overlaid figure...")

    categories = list(short_tables.keys())
    total_span = upstream + body_length + downstream

    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT / 1.4),
                             sharey=True)
    fig.suptitle("Metagene profiles — all platforms overlaid (min-max normalised)",
                 fontsize=TITLE_SIZE + 2, y=1.02)

    panel_configs = [
        (axes[0], short_tables, "Short genes (bottom 10% exon union)"),
        (axes[1], long_tables,  "Long genes (top 10% exon union)"),
    ]

    for ax, tables, panel_title in panel_configs:
        for category in categories:
            x, y   = parse_profile_table(tables[category])
            y_norm = minmax_normalize(y)
            color  = category_colors.get(category, "steelblue")

            ax.plot(x, y_norm, color=color, linewidth=LINE_WIDTH, label=category)
            ax.fill_between(x, y_norm, alpha=0.08, color=color)

        # TSS / TES markers (derived from the last x array — all categories share bins)
        n_bins  = len(x)
        tss_bin = round(n_bins * upstream / total_span)
        tes_bin = round(n_bins * (upstream + body_length) / total_span) - 1

        ax.axvline(x=x[tss_bin], color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(x=x[tes_bin], color="gray", linestyle="--", linewidth=0.8)

        ax.set_xticks([x[tss_bin], x[tes_bin]])
        ax.set_xticklabels(["TSS", "TES"], fontsize=FONT_SIZE - 2,
                           rotation=15, ha="right")
        ax.set_ylim(-0.05, 1.15)
        ax.tick_params(labelsize=FONT_SIZE - 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title(panel_title, fontsize=TITLE_SIZE, fontweight="bold")
        ax.set_xlabel("Genomic position", fontsize=FONT_SIZE - 1)

    axes[0].set_ylabel("Normalised coverage (0–1)", fontsize=FONT_SIZE - 1)

    # Single shared legend, placed outside the right panel
    handles, labels = axes[0].get_legend_handles_labels()
    axes[1].legend(handles, labels,
                   fontsize=FONT_SIZE - 1,
                   frameon=False,
                   loc="upper left",
                   bbox_to_anchor=(1.02, 1))

    fig.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    for ext in ["png", "pdf", "svg"]:
        out_path = os.path.join(output_dir, f"metagene_overlaid_1x2.{ext}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"    Saved: {out_path}")

    plt.close(fig)

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Metagene gene-body coverage profiles: short/long x N categories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Entry points (mutually exclusive, from slowest to fastest):
  [default]           Full run from BAM files
  --skip-bam          Skip bamCoverage; reuse existing BigWigs
  --from-matrices     Skip bamCoverage + computeMatrix; rerun plotProfile + figure
  --from-plotprofile  Skip everything; reuse existing .tab files and redraw figure
        """
    )
    parser.add_argument("--gtf",        default=GTF_FILE)
    parser.add_argument("--outdir",     default=OUTPUT_DIR)
    parser.add_argument("--threads",    default=THREADS,      type=int)
    parser.add_argument("--binsize",    default=BINSIZE,      type=int)
    parser.add_argument("--normalize",  default=NORMALIZE,
                        choices=["CPM", "RPKM", "BPM", "RPGC", "None"])
    parser.add_argument("--body",       default=BODY_LENGTH,  type=int,
                        help="Scaled gene body length for computeMatrix (bp)")
    parser.add_argument("--upstream",   default=UPSTREAM,     type=int)
    parser.add_argument("--downstream", default=DOWNSTREAM,   type=int)
    parser.add_argument("--min-body",   default=200,           type=int,
                        help="Exclude genes with gene span shorter than this (bp, default: 200)")
    parser.add_argument("--min-coverage", default=5.0, type=float,
                    help="Min mean CPM coverage across gene body in all categories (default: 5)")
    parser.add_argument("--region-mode", default="exon",
                    choices=["exon", "gene_body"],
                    help="'exon': plot longest-transcript exons only (default). "
                         "'gene_body': plot full TSS->TES span (includes introns).")
    
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--skip-bam",         action="store_true",
                      help="Reuse existing BigWigs; redo BED + computeMatrix + figure")
    mode.add_argument("--from-matrices",    action="store_true",
                      help="Reuse existing .gz matrices; redo plotProfile + figure")
    mode.add_argument("--from-plotprofile", action="store_true",
                      help="Reuse existing .tab profile files; redraw figure only")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    categories = list(BAMS.keys())

    if args.from_plotprofile:
        print("[Resuming from existing .tab profile files — figure only]")
        short_tables, long_tables = get_existing_profile_tables(args.outdir, categories)

    elif args.from_matrices:
        print("[Resuming from existing .gz matrices — plotProfile + figure]")
        short_tables, long_tables = run_plotprofile_only(args.outdir, categories)

    elif args.skip_bam:
        print("[Resuming from existing BigWigs — BED + computeMatrix + figure]")
        avg_bws = get_existing_bigwigs(args.outdir)
        
        gene_df, exon_df = compute_exon_union_lengths(args.gtf)

        if args.region_mode == "exon":
            exon_df = select_longest_transcripts(exon_df)

        gene_df = filter_genes_by_coverage(gene_df, exon_df, avg_bws,
                                        min_mean_coverage=args.min_coverage)
        short_bed, long_bed = write_bed_files(gene_df, exon_df, args.outdir,
                                            args.min_body, args.region_mode)

        print("[4/5] Running computeMatrix + plotProfile per category...")
        short_tables = run_compute_matrix(
            avg_bws, short_bed, "short",
            args.outdir, args.body, args.upstream, args.downstream, args.threads
        )
        long_tables = run_compute_matrix(
            avg_bws, long_bed, "long",
            args.outdir, args.body, args.upstream, args.downstream, args.threads
        )

    else:
        gene_df, exon_df = compute_exon_union_lengths(args.gtf)
        avg_bws = run_bamCoverage(
            BAMS, args.outdir, args.binsize, args.normalize, args.threads
        )
        gene_df, exon_df = compute_exon_union_lengths(args.gtf)

        if args.region_mode == "exon":
            exon_df = select_longest_transcripts(exon_df)

        gene_df = filter_genes_by_coverage(gene_df, exon_df, avg_bws,
                                        min_mean_coverage=args.min_coverage)
        short_bed, long_bed = write_bed_files(gene_df, exon_df, args.outdir,
                                            args.min_body, args.region_mode)

        print("[4/5] Running computeMatrix + plotProfile per category...")
        short_tables = run_compute_matrix(
            avg_bws, short_bed, "short",
            args.outdir, args.body, args.upstream, args.downstream, args.threads
        )
        long_tables = run_compute_matrix(
            avg_bws, long_bed, "long",
            args.outdir, args.body, args.upstream, args.downstream, args.threads
        )

    build_figure(
        short_tables, long_tables,
        CATEGORY_COLORS, args.outdir,
        args.upstream, args.downstream, args.body
    )

    build_overlaid_figure(
        short_tables, long_tables,
        CATEGORY_COLORS, args.outdir,
        args.upstream, args.downstream, args.body
    )

    print(f"\n[Done] Final figure saved to: {args.outdir}/metagene_genebody_2xN.png/.pdf/.svg")


if __name__ == "__main__":
    main()
