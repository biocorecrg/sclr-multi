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

# Ensembl v111 reference URLs
ENSEMBL_RELEASE="111"
declare -A REF_FASTA_URLS=(
    ["human"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    ["mouse"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    ["dog"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz"
)

declare -A REF_GTF_URLS=(
    ["human"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"
    ["mouse"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz"
    ["dog"]="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.${ENSEMBL_RELEASE}.gtf.gz"
)

usage() {
    echo "Usage: $0 [TARGET] [OUTPUT_DIR] [--dryrun] [--skip-ref] [--skip-data]"
    echo ""
    echo "Targets:"
    echo "  10X-3prime   - Download & subsample 10X 3' Rep1 (1M reads, ERR17793697)"
    echo "  10X-5prime   - Download & subsample 10X 5' Rep1 (1M reads, ERR17793700)"
    echo "  Argentag     - Download & subsample Argentag Rep1 (1M reads, ERR17793703)"
    echo "  Parse        - Download & subsample Parse Rep1 (2M reads, ERR17793706)"
    echo "  references   - Download & concatenate reference genomes (Human, Mouse, Dog) & GTFs only"
    echo "  all          - Download & subsample all test datasets and references (default)"
    echo ""
    echo "Options:"
    echo "  --dryrun, -n, dryrun   Skip downloading/subsampling and only generate samplesheets & execute.sh"
    echo "  --skip-data            Skip FASTQ dataset download & subsampling (only download reference databases)"
    echo "  --skip-ref             Skip reference genome and annotation download"
    echo ""
    echo "Examples:"
    echo "  $0 --skip-data test_data"
    echo "  $0 references test_data"
    echo "  $0 Argentag test_data"
    echo "  $0 all test_data --dryrun"
    echo "  $0 --dryrun"
    exit 1
}

TARGET="all"
DEST_DIR="test_data"
DRY_RUN=false
DOWNLOAD_REF=true
DOWNLOAD_DATA=true

# Parse arguments flexibly
for arg in "$@"; do
    case "$arg" in
        --dryrun|-n|dryrun)
            DRY_RUN=true
            ;;
        --skip-ref|--skip-references|skip-ref)
            DOWNLOAD_REF=false
            ;;
        --skip-data|--skip-sample|--skip-samples|skip-data)
            DOWNLOAD_DATA=false
            ;;
        -h|--help|help)
            usage
            ;;
        10X-3prime|10X-5prime|Argentag|Parse|references|all)
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

# ------------------------------------------------------------------------------
# Download and concatenate reference genomes and GTF annotations
# ------------------------------------------------------------------------------
process_references() {
    local COMBINED_FASTA="${DEST_DIR}/genome.fa"
    local COMBINED_GTF="${DEST_DIR}/genes.gtf.gz"

    echo "=========================================================="
    echo "Processing Reference Genomes & Annotations (Ensembl v${ENSEMBL_RELEASE})"
    echo "Species             : Human (GRCh38), Mouse (GRCm39), Dog (ROS_Cfam_1.0)"
    echo "Target FASTA        : ${COMBINED_FASTA}"
    echo "Target GTF          : ${COMBINED_GTF}"
    echo "=========================================================="

    if [[ "${DRY_RUN}" == true ]]; then
        echo "--> [DRYRUN] Skipping reference download and concatenation."
        return 0
    fi

    # Define species prefix mapping
    declare -A PREFIXES=(
        ["human"]="hum_"
        ["mouse"]="mou_"
        ["dog"]="dog_"
    )

    # 1. FASTA concatenation with prefixing
    if [[ -f "${COMBINED_FASTA}" ]]; then
        echo "--> Combined genome FASTA ${COMBINED_FASTA} already exists. Skipping."
    else
        echo "--> Downloading, prefixing chromosome names, and concatenating reference FASTA files..."
        local TMP_FASTA="${COMBINED_FASTA}.tmp"
        rm -f "${TMP_FASTA}"
        for SPECIES in "human" "mouse" "dog"; do
            local URL="${REF_FASTA_URLS[$SPECIES]}"
            local PREFIX="${PREFIXES[$SPECIES]}"
            echo "   -> Downloading and prefixing ${SPECIES} FASTA (prefix: '${PREFIX}'): ${URL}"
            curl -sSL --fail "${URL}" | gzip -dc | sed "s/^>/>${PREFIX}/" >> "${TMP_FASTA}"
        done
        mv "${TMP_FASTA}" "${COMBINED_FASTA}"
        echo "--> Successfully created combined FASTA ${COMBINED_FASTA} ($(du -h "${COMBINED_FASTA}" | awk '{print $1}'))!"
    fi

    # 2. GTF concatenation with prefixing
    if [[ -f "${COMBINED_GTF}" ]]; then
        echo "--> Combined GTF ${COMBINED_GTF} already exists. Skipping."
    else
        echo "--> Downloading, prefixing chromosome names, and concatenating GTF annotation files..."
        local TMP_GTF="${COMBINED_GTF}.tmp"
        rm -f "${TMP_GTF}"
        for SPECIES in "human" "mouse" "dog"; do
            local URL="${REF_GTF_URLS[$SPECIES]}"
            local PREFIX="${PREFIXES[$SPECIES]}"
            echo "   -> Downloading and prefixing ${SPECIES} GTF (prefix: '${PREFIX}'): ${URL}"
            curl -sSL --fail "${URL}" | gzip -dc | sed -e "/^#/!s/^/${PREFIX}/" >> "${TMP_GTF}"
        done
        gzip -c "${TMP_GTF}" > "${COMBINED_GTF}"
        rm -f "${TMP_GTF}"
        echo "--> Successfully created combined GTF ${COMBINED_GTF} ($(du -h "${COMBINED_GTF}" | awk '{print $1}'))!"
    fi
}

if [[ "${DOWNLOAD_REF}" == true && ( "${TARGET}" == "all" || "${TARGET}" == "references" ) ]]; then
    process_references
fi

if [[ "${TARGET}" == "references" ]]; then
    echo "Reference preparation complete."
    exit 0
fi

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

    if [[ "${DRY_RUN}" == true || "${DOWNLOAD_DATA}" == false ]]; then
        echo "--> [SKIP-DATA] Skipping FASTQ download and subsampling for ${PLATFORM}."
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

declare -A PARAMS_FILES=(
    ["10X-3prime"]="params_10X_3prime.yaml"
    ["10X-5prime"]="params_10X_5prime.yaml"
    ["Argentag"]="params_argentag.yaml"
    ["Parse"]="params_parse.yaml"
)

print_run_command() {
    local PLATFORM="$1"
    local SAMPLESHEET_NAME="samplesheet_${PLATFORM}.csv"
    local PARAMS_NAME="${PARAMS_FILES[$PLATFORM]}"

    echo "----------------------------------------------------------"
    echo "To run ${PLATFORM}:"
    echo "   nextflow run main.nf \\"
    echo "      -profile singularity,newcrg \\"
    echo "      -params-file ${DEST_DIR}/${PARAMS_NAME} \\"
    echo "      --genome_fasta ${DEST_DIR}/genome.fa \\"
    echo "      --gtf ${DEST_DIR}/genes.gtf.gz \\"
    echo "      --input ${DEST_DIR}/${SAMPLESHEET_NAME} \\"
    echo "      --outdir results_${PLATFORM}"
    echo "----------------------------------------------------------"
}

generate_params_file() {
    local PLATFORM="$1"
    local PARAMS_NAME="${PARAMS_FILES[$PLATFORM]}"
    local PARAMS_PATH="${DEST_DIR}/${PARAMS_NAME}"
    local PLATFORM_VAL="10X"
    local BC_FORMAT=""

    if [[ "${PLATFORM}" == "10X-3prime" ]]; then
        PLATFORM_VAL="10X"
        BC_FORMAT="10X_3v4"
    elif [[ "${PLATFORM}" == "10X-5prime" ]]; then
        PLATFORM_VAL="10X"
        BC_FORMAT="10X_5v3"
    elif [[ "${PLATFORM}" == "Argentag" ]]; then
        PLATFORM_VAL="Argentag"
        BC_FORMAT=""
    elif [[ "${PLATFORM}" == "Parse" ]]; then
        PLATFORM_VAL="Parse"
        BC_FORMAT=""
    fi

    cat <<EOF > "${PARAMS_PATH}"
## Input options
input                               : "samplesheet.csv"
platform                            : "${PLATFORM_VAL}"
outdir                              : "./results"

## References
genome                              : null
genome_fasta                        : "genome.fa"
transcript_fasta                    : null
gtf                                 : "genes.gtf.gz"
igenomes_base                       : 's3://ngi-igenomes/igenomes/'
igenomes_ignore                     : true
fasta_delimiter                     : ""

## Fastq Options
split_amount                        : 12000000

## Read Trimming Options
min_length                          : 0
min_q_score                         : 0
skip_trimming                       : true

## 10X cell barcode options
whitelist                           : null
barcode_format                      : "${BC_FORMAT}"

## Deduplication
dedup_tool                          : "umitools"

## Parse related options
parse_generate_pe                   : "--chemistry v3 --l1dist 3 --l2dist 2"
parse_spipe                         : "--chemistry v3 --kit WT_mini"

## Argentag related options
argentag_taggy_demux                : "--out-fmt='flames' --trim-TSO --preserve --trim-poly 'normal' --orient 'sense' --keep-failed"

## Library strandness option
stranded                            : null

## Mapping
skip_save_minimap2_index            : false
kmer_size                           : 14

## Analysis options
retain_introns                      : true
quantifier                          : "isoquant"

## Skipping options:
skip_qc                             : false
skip_seurat                         : false
skip_nanoplot                       : true
skip_toulligqc                      : true
skip_fastqc                         : false
skip_bam_nanocomp                   : true
EOF
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
        # Also generate the samplesheet and params yaml in DEST_DIR automatically
        local_sample_file="${DEST_DIR}/samplesheet_${PLATFORM}.csv"
        fastq_file="${DEST_DIR}/${PLATFORM}.fastq.gz"
        echo "sample,fastq,cell_count" > "${local_sample_file}"
        echo "${PLATFORM},${fastq_file},10000" >> "${local_sample_file}"
        generate_params_file "${PLATFORM}"
        print_run_command "${PLATFORM}"

        PARAMS_NAME="${PARAMS_FILES[$PLATFORM]}"
        cat <<EOF >> "${EXECUTE_SCRIPT}"
nextflow run main.nf \\
      -profile singularity,newcrg \\
      -params-file ${DEST_DIR}/${PARAMS_NAME} \\
      --genome_fasta ${DEST_DIR}/genome.fa \\
      --gtf ${DEST_DIR}/genes.gtf.gz \\
      --input ${DEST_DIR}/samplesheet_${PLATFORM}.csv \\
      --outdir results_${PLATFORM}

EOF
    done
else
    local_sample_file="${DEST_DIR}/samplesheet_${TARGET}.csv"
    fastq_file="${DEST_DIR}/${TARGET}.fastq.gz"
    echo "sample,fastq,cell_count" > "${local_sample_file}"
    echo "${TARGET},${fastq_file},10000" >> "${local_sample_file}"
    generate_params_file "${TARGET}"
    print_run_command "${TARGET}"

    PARAMS_NAME="${PARAMS_FILES[$TARGET]}"
    cat <<EOF >> "${EXECUTE_SCRIPT}"
nextflow run main.nf \\
      -profile singularity,newcrg \\
      -params-file ${DEST_DIR}/${PARAMS_NAME} \\
      --genome_fasta ${DEST_DIR}/genome.fa \\
      --gtf ${DEST_DIR}/genes.gtf.gz \\
      --input ${DEST_DIR}/samplesheet_${TARGET}.csv \\
      --outdir results_${TARGET}
EOF
fi

chmod +x "${EXECUTE_SCRIPT}"
echo "Now you can execute by sourcing ${EXECUTE_SCRIPT}"


