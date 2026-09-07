#!/usr/bin/env bash
# ==============================================================================
# Script: download_and_subsample_test_data.sh
# Description: Downloads Rep1 FASTQ from ENA (PRJEB125152) and subsamples reads
#              using seqkit for testing scLR-multi.
# ==============================================================================

set -euo pipefail

declare -A SAMPLES=(
    ["10X-3prime"]="ERR17793697 1000000 10X-3prime.fastq.gz"
    ["10X-5prime"]="ERR17793700 1000000 10X-5prime.fastq.gz"
    ["Argentag"]="ERR17793703 1000000 Argentag.fastq.gz"
    ["Parse"]="ERR17793706 2000000 Parse.fastq.gz"
)

usage() {
    echo "Usage: $0 [TARGET] [OUTPUT_DIR] [--dryrun]"
    echo ""
    echo "Targets:"
    echo "  10X-3prime   - Download & subsample 10X 3' Rep1 (1M reads, ERR17793697)"
    echo "  10X-5prime   - Download & subsample 10X 5' Rep1 (1M reads, ERR17793700)"
    echo "  Argentag     - Download & subsample Argentag Rep1 (1M reads, ERR17793703)"
    echo "  Parse        - Download & subsample Parse Rep1 (2M reads, ERR17793706)"
    echo "  all          - Download & subsample all test datasets (default)"
    echo ""
    echo "Options:"
    echo "  --dryrun, -n, dryrun   Skip downloading/subsampling and only generate samplesheets & execute.sh"
    echo ""
    echo "Examples:"
    echo "  $0 Argentag test_data"
    echo "  $0 all test_data --dryrun"
    echo "  $0 --dryrun"
    exit 1
}

TARGET="all"
DEST_DIR="test_data"
DRY_RUN=false

# Parse arguments flexibly
for arg in "$@"; do
    case "$arg" in
        --dryrun|-n|dryrun)
            DRY_RUN=true
            ;;
        -h|--help|help)
            usage
            ;;
        10X-3prime|10X-5prime|Argentag|Parse|all)
            TARGET="$arg"
            ;;
        *)
            if [[ -z "${CUSTOM_DEST:-}" && "$arg" != -* ]]; then
                DEST_DIR="$arg"
                CUSTOM_DEST=1
            else
                echo "Error: Unknown argument '$arg'"
                usage
            fi
            ;;
    esac
done

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

    if [[ "${DRY_RUN}" == true ]]; then
        echo "--> [DRYRUN] Skipping download and subsampling for ${PLATFORM}."
        return 0
    fi

    if [[ -f "${FINAL_FASTQ}" ]]; then
        echo "--> Output file ${FINAL_FASTQ} already exists. Skipping."
        return 0
    fi

    # Number of FASTQ lines = reads * 4
    local LINES=$(( READS * 4 ))
    local TEMP_FASTQ="${FINAL_FASTQ}.tmp.gz"

    echo "--> Streaming and extracting first ${READS} reads directly from ENA..."

    # Use curl + gzip pipeline to extract only the needed lines and terminate download
    set +e
    curl -sSL --fail "${HTTP_URL}" | gzip -dc 2>/dev/null | head -n "${LINES}" | gzip > "${TEMP_FASTQ}"
    local CURL_STATUS=$?
    set -e

    # If stream produced the file, finalize it
    if [[ -s "${TEMP_FASTQ}" ]]; then
        mv "${TEMP_FASTQ}" "${FINAL_FASTQ}"
        echo "--> Successfully created ${FINAL_FASTQ} ($(du -h "${FINAL_FASTQ}" | awk '{print $1}'))!"
    else
        echo "--> HTTP stream failed or was empty, trying FTP URL: ${ENA_URL}..."
        set +e
        curl -sSL --fail "${ENA_URL}" | gzip -dc 2>/dev/null | head -n "${LINES}" | gzip > "${TEMP_FASTQ}"
        set -e
        if [[ -s "${TEMP_FASTQ}" ]]; then
            mv "${TEMP_FASTQ}" "${FINAL_FASTQ}"
            echo "--> Successfully created ${FINAL_FASTQ} ($(du -h "${FINAL_FASTQ}" | awk '{print $1}'))!"
        else
            echo "Error: Failed to stream from ${HTTP_URL} and ${ENA_URL}" >&2
            rm -f "${TEMP_FASTQ}"
            exit 1
        fi
    fi
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

print_run_command() {
    local PLATFORM="$1"
    local SAMPLESHEET_NAME="samplesheet_${PLATFORM}.csv"
    local PARAMS_NAME="params_${PLATFORM}.yaml"

    echo "----------------------------------------------------------"
    echo "To run ${PLATFORM}:"
    echo "   nextflow run main.nf \\"
    echo "      -profile singularity,newcrg \\"
    echo "      -params-file test_input/${PARAMS_NAME} \\"
    echo "      --input ${DEST_DIR}/${SAMPLESHEET_NAME} \\"
    echo "      --outdir results_${PLATFORM}"
    echo "----------------------------------------------------------"
}

EXECUTE_SCRIPT="execute.sh"
echo "#!/usr/bin/env bash" > "${EXECUTE_SCRIPT}"
echo "# Auto-generated execution script for scLR-multi test runs" >> "${EXECUTE_SCRIPT}"
echo "set -euo pipefail" >> "${EXECUTE_SCRIPT}"
echo "" >> "${EXECUTE_SCRIPT}"

echo ""
echo "=========================================================="
echo "Done! Subsampled test datasets are ready in '${DEST_DIR}'."
echo "=========================================================="
echo ""
echo "Pipeline execution commands (written to ${EXECUTE_SCRIPT}):"
echo ""

if [[ "${TARGET}" == "all" ]]; then
    for PLATFORM in "10X-3prime" "10X-5prime" "Argentag" "Parse"; do
        # Also generate the samplesheet in DEST_DIR automatically
        local_sample_file="${DEST_DIR}/samplesheet_${PLATFORM}.csv"
        fastq_file="${DEST_DIR}/${PLATFORM}.fastq.gz"
        echo "sample,fastq,cell_count" > "${local_sample_file}"
        echo "${PLATFORM},${fastq_file},10000" >> "${local_sample_file}"
        print_run_command "${PLATFORM}"

        PARAMS_NAME="params_${PLATFORM}.yaml"
        cat <<EOF >> "${EXECUTE_SCRIPT}"
nextflow run main.nf \\
      -profile singularity,newcrg \\
      -params-file test_input/${PARAMS_NAME} \\
      --input ${DEST_DIR}/samplesheet_${PLATFORM}.csv \\
      --outdir results_${PLATFORM}

EOF
    done
else
    local_sample_file="${DEST_DIR}/samplesheet_${TARGET}.csv"
    fastq_file="${DEST_DIR}/${TARGET}.fastq.gz"
    echo "sample,fastq,cell_count" > "${local_sample_file}"
    echo "${TARGET},${fastq_file},10000" >> "${local_sample_file}"
    print_run_command "${TARGET}"

    PARAMS_NAME="params_${TARGET}.yaml"
    cat <<EOF >> "${EXECUTE_SCRIPT}"
nextflow run main.nf \\
      -profile singularity,newcrg \\
      -params-file test_input/${PARAMS_NAME} \\
      --input ${DEST_DIR}/samplesheet_${TARGET}.csv \\
      --outdir results_${TARGET}
EOF
fi

chmod +x "${EXECUTE_SCRIPT}"
echo "Now you can execute by sourcing ${EXECUTE_SCRIPT}"


