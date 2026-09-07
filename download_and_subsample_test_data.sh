#!/usr/bin/env bash
# ==============================================================================
# Script: download_and_subsample_test_data.sh
# Description: Downloads Rep1 FASTQ from ENA (PRJEB125152) and subsamples reads
#              using seqkit for testing scLR-multi.
# ==============================================================================

set -euo pipefail

DEST_DIR="${1:-test_data}"
mkdir -p "${DEST_DIR}"

declare -A SAMPLES=(
    ["10X-3prime"]="ERR17793697 1000000 10X-3prime.fastq.gz"
    ["10X-5prime"]="ERR17793700 1000000 10X-5prime.fastq.gz"
    ["Argentag"]="ERR17793703 1000000 Argentag.fastq.gz"
    ["Parse"]="ERR17793706 2000000 Parse.fastq.gz"
)

usage() {
    echo "Usage: $0 [TARGET] [OUTPUT_DIR]"
    echo ""
    echo "Targets:"
    echo "  10X-3prime   - Download & subsample 10X 3' Rep1 (1M reads, ERR17793697)"
    echo "  10X-5prime   - Download & subsample 10X 5' Rep1 (1M reads, ERR17793700)"
    echo "  Argentag     - Download & subsample Argentag Rep1 (1M reads, ERR17793703)"
    echo "  Parse        - Download & subsample Parse Rep1 (2M reads, ERR17793706)"
    echo "  all          - Download & subsample all test datasets"
    echo ""
    echo "Example:"
    echo "  $0 10X-3prime test_data"
    exit 1
}

TARGET="${1:-all}"
if [[ "$#" -ge 2 ]]; then
    DEST_DIR="$2"
else
    DEST_DIR="test_data"
fi

mkdir -p "${DEST_DIR}"

process_sample() {
    local PLATFORM="$1"
    local INFO="${SAMPLES[$PLATFORM]}"
    local ACCESSION=$(echo "$INFO" | awk '{print $1}')
    local READS=$(echo "$INFO" | awk '{print $2}')
    local OUT_NAME=$(echo "$INFO" | awk '{print $3}')

    # Construct ENA FTP URL (e.g. ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR177/097/ERR17793697/ERR17793697.fastq.gz)
    local VOL="${ACCESSION:0:6}"
    local SUFFIX=$(printf "%03d" $(( 10#${ACCESSION: -3} % 1000 )))
    local ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${VOL}/0${SUFFIX: -2}/${ACCESSION}/${ACCESSION}.fastq.gz"
    # Fallback to standard HTTP download URL
    local HTTP_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/${VOL}/0${SUFFIX: -2}/${ACCESSION}/${ACCESSION}.fastq.gz"

    local RAW_FASTQ="${DEST_DIR}/${ACCESSION}.fastq.gz"
    local FINAL_FASTQ="${DEST_DIR}/${OUT_NAME}"

    echo "=========================================================="
    echo "Processing Platform : ${PLATFORM}"
    echo "ENA Accession       : ${ACCESSION}"
    echo "Target reads        : ${READS}"
    echo "Output FASTQ        : ${FINAL_FASTQ}"
    echo "=========================================================="

    if [[ -f "${FINAL_FASTQ}" ]]; then
        echo "--> Output file ${FINAL_FASTQ} already exists. Skipping."
        return 0
    fi

    if [[ ! -f "${RAW_FASTQ}" ]]; then
        echo "--> Downloading ${ACCESSION} from ENA..."
        if command -v wget &> /dev/null; then
            wget -c "${HTTP_URL}" -O "${RAW_FASTQ}" || wget -c "${ENA_URL}" -O "${RAW_FASTQ}"
        elif command -v curl &> /dev/null; then
            curl -L -C - "${HTTP_URL}" -o "${RAW_FASTQ}" || curl -L -C - "${ENA_URL}" -o "${RAW_FASTQ}"
        else
            echo "Error: Neither wget nor curl found on PATH." >&2
            exit 1
        fi
    else
        echo "--> Raw file ${RAW_FASTQ} already downloaded."
    fi

    echo "--> Subsampling ${READS} reads with seqkit (seed: 42)..."
    if command -v seqkit &> /dev/null; then
        seqkit sample -n "${READS}" -s 42 "${RAW_FASTQ}" -o "${FINAL_FASTQ}"
    elif command -v singularity &> /dev/null; then
        echo "Using seqkit via Singularity container..."
        singularity exec docker://quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0 \
            seqkit sample -n "${READS}" -s 42 "${RAW_FASTQ}" -o "${FINAL_FASTQ}"
    else
        echo "Error: seqkit is not installed and Singularity is unavailable." >&2
        exit 1
    fi

    echo "--> Finished creating ${FINAL_FASTQ}!"
}

if [[ "${TARGET}" == "all" ]]; then
    for PLATFORM in "10X-3prime" "10X-5prime" "Argentag" "Parse"; do
        process_sample "${PLATFORM}"
    done
elif [[ -n "${SAMPLES[${TARGET}]:-}" ]]; then
    process_sample "${TARGET}"
else
    echo "Error: Unknown target '${TARGET}'"
    usage
fi

echo ""
echo "Done! Subsampled test datasets are ready in '${DEST_DIR}'."
