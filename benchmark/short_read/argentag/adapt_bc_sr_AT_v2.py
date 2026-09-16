#!/usr/bin/env python3
import os
import sys
import json
import gzip
import glob
import shutil
from collections import Counter

def ensure_dir(directory):
    """Creates directory if it does not exist."""
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

def validate_inputs(input_dir):
    """Checks if required R1 and R2 files exist in the input directory."""
    r1_exists = any(os.path.isfile(os.path.join(input_dir, f)) for f in ["R1.fastq", "R1.fastq.gz"])
    r2_exists = any(os.path.isfile(os.path.join(input_dir, f)) for f in ["R2.fastq", "R2.fastq.gz"])
    
    if not r1_exists or not r2_exists:
        print(f"Error: Missing required R1 or R2 files in {input_dir}")
        sys.exit(1)
    print("[*] Input validation successful.")

def run_step1_count(input_file):
    """Step 1: Counts barcodes from original R1 (handling .gz if needed)."""
    print(f"[*] Step 1: Counting barcodes...")
    counts = Counter()
    open_func = gzip.open if input_file.endswith('.gz') else open
    with open_func(input_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1: # Sequence line
                bc = line.strip()[:16]
                counts[bc] += 1
    return counts

def run_step2_mapping(counts, output_dir, whitelist_path):
    """Step 2: Maps observed BCs to whitelist and saves files."""
    print(f"[*] Step 2: Mapping barcodes to whitelist...")
    ensure_dir(output_dir)
    
    # Sort barcodes by frequency
    count_barcodes = [bc for bc, _ in counts.most_common()]
    
    bc_counts_path = os.path.join(output_dir, "bc_counts.txt")
    with open(bc_counts_path, 'w') as f_bc:
        for bc, freq in counts.most_common():
            f_bc.write(f"{bc}\t{freq}\n")

    with gzip.open(whitelist_path, 'rt') as f:
        ref_barcodes = [line.strip() for line in f]
            
    barcode_map = {obs: ref_barcodes[i] for i, obs in enumerate(count_barcodes) if i < len(ref_barcodes)}
    
    json_out = os.path.join(output_dir, 'barcode_mapping.json')
    with open(json_out, 'w') as f:
        json.dump(barcode_map, f, indent=2)
    return barcode_map

def process_r1_final(input_file, output_file, barcode_map):
    """
    Combines Step 0 and Step 3:
    Fixes header of R1 from taggy, AND transforms sequence/quality in one pass to avoid duplication.
    """
    print(f"[*] Step 3: Transforming R1 (Fixing headers + Mapping barcodes)...")
    ensure_dir(os.path.dirname(output_file))
    
    open_func = gzip.open if input_file.endswith('.gz') else open
    with open_func(input_file, 'rt') as fin, open(output_file, 'w') as fout:
        while True:
            header = fin.readline().strip()
            if not header: break
            seq = fin.readline().strip()
            plus = fin.readline().strip()
            qual = fin.readline().strip()
            
            # Step 0 Logic: Fix header
            header = header.replace(" 2:N:0:", " 1:N:0:") # This error comes from taggy 3.2.5d
            
            # Step 3 Logic: Transform
            if len(seq) >= 28:
                virtual_bc = seq[:16]
                if virtual_bc in barcode_map:
                    new_seq = barcode_map[virtual_bc] + seq[-12:]
                    new_qual = ("I" * 16) + qual[-12:]
                    fout.write(f"{header}\n{new_seq}\n{plus}\n{new_qual}\n")
                    continue
            
            fout.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

def generate_indices(r1_file, i1_file, i2_file):
    """Step 4: Generates I1 and I2 index files from the processed R1."""
    print(f"[*] Step 4: Generating Index files (I1 and I2)...")
    ensure_dir(os.path.dirname(i1_file))
    
    # I1: CCCACCACAA with 1:N:0 (R1 header)
    # I2: AAGCGGAGGG with 2:N:0 (R2 header)
    
    with open(r1_file, 'r') as fin, open(i1_file, 'w') as f_i1, open(i2_file, 'w') as f_i2:
        for i, line in enumerate(fin):
            mod = i % 4
            if mod == 0:  # Header
                header_r1 = line.strip()
                header_i2 = header_r1.replace(" 1:N:0:", " 2:N:0:")
                f_i1.write(header_r1 + "\n")
                f_i2.write(header_i2 + "\n")
            elif mod == 1:  # Sequence
                f_i1.write("CCCACCACAA\n")
                f_i2.write("AAGCGGAGGG\n")
            elif mod == 2:  # Plus
                f_i1.write("+\n")
                f_i2.write("+\n")
            elif mod == 3:  # Quality
                f_i1.write("IIIIIIIIII\n")
                f_i2.write("IIIIIIIIII\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 process_barcode_pipeline.py <input_dir> <output_dir> <sample_name> <whitelist_path>")
        sys.exit(1)

    in_base = sys.argv[1]
    out_base = sys.argv[2]
    sample_name = sys.argv[3]
    whitelist_path = sys.argv[4]

    validate_inputs(in_base)

    # Path Setup
    counts_dir = os.path.join(out_base, "01_counts")
    mapping_dir = os.path.join(out_base, "02_mapping")
    final_dir = os.path.join(out_base, "03_cellranger_ready")

    # Determine input R1 path
    raw_r1 = os.path.join(in_base, "R1.fastq")
    if not os.path.exists(raw_r1):
        raw_r1 = os.path.join(in_base, "R1.fastq.gz")

    # Step 1: Count directly from source
    counts = run_step1_count(raw_r1)

    # Step 2: Create map
    barcode_map = run_step2_mapping(counts, mapping_dir, whitelist_path)

    # Step 3: Finalize R1 (Header fix + BC mapping in one go)
    final_r1_path = os.path.join(final_dir, f"{sample_name}_S1_L001_R1_001.fastq")
    process_r1_final(raw_r1, final_r1_path, barcode_map)

    # Step 4: Finalize R2 and Indices
    print(f"[*] Finalizing Cell Ranger folder...")
    
    # Process R2: Decompress if necessary while copying to final destination
    raw_r2 = os.path.join(in_base, "R2.fastq")
    if not os.path.exists(raw_r2):
        raw_r2 = os.path.join(in_base, "R2.fastq.gz")
    
    dest_r2 = os.path.join(final_dir, f"{sample_name}_S1_L001_R2_001.fastq")
    if raw_r2.endswith('.gz'):
        with gzip.open(raw_r2, 'rt') as f_in, open(dest_r2, 'w') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(raw_r2, dest_r2)

    # Generate I1 and I2
    i1_path = os.path.join(final_dir, f"{sample_name}_S1_L001_I1_001.fastq")
    i2_path = os.path.join(final_dir, f"{sample_name}_S1_L001_I2_001.fastq")
    generate_indices(final_r1_path, i1_path, i2_path)

    print(f"\n[SUCCESS] Pipeline completed. Output in: {final_dir}")

if __name__ == "__main__":
    main()
