/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { NANOFILT                                                    } from "../modules/local/nanofilt"
include { SPLIT_FILE                                                  } from "../modules/local/split_file"
include { SPLIT_FILE as SPLIT_FILE_BC_FASTQ                           } from "../modules/local/split_file"
include { SPLIT_FILE as SPLIT_FILE_BC_CSV                             } from "../modules/local/split_file"
include { SPLIT_FILE as SPLIT_FILE_PARSE                              } from "../modules/local/split_file"
include { SPLIT_FILE as SPLIT_FILE_ARGENTAG                           } from "../modules/local/split_file"
include { BLAZE                                                       } from "../modules/local/blaze"
include { PREEXTRACT_FASTQ                                            } from "../modules/local/preextract_fastq.nf"
include { READ_COUNTS as READ_COUNTS_GENOME                           } from "../modules/local/read_counts.nf"
include { READ_COUNTS as READ_COUNTS_TRANSCRIPT                       } from "../modules/local/read_counts.nf"
include { CORRECT_BARCODES                                            } from "../modules/local/correct_barcodes"
include { UCSC_GTFTOGENEPRED                                          } from "../modules/local/ucsc_gtftogenepred"
include { UCSC_GENEPREDTOBED                                          } from "../modules/local/ucsc_genepredtobed"

include { GENERATE_PE                                                 } from "../modules/local/generate_PE/main"
include { SPLITPIPE_PRE                                               } from "../modules/local/split_pipe_pre/main"

include { ARGENTAG_SPLIT                                              } from "../modules/local/argentag/split/main"
include { ARGENTAG_TAGGY_DEMUX                                        } from "../modules/local/argentag/taggy_demux/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_REFERENCE_FILES                                     } from "../subworkflows/local/prepare_reference_files"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_GENOME     } from "../subworkflows/local/process_longread_scrna"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_TRANSCRIPT } from "../subworkflows/local/process_longread_scrna"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//
include { PIGZ_UNCOMPRESS as GUNZIP_FASTQ                             } from "../modules/nf-core/pigz/uncompress/main"
include { PIGZ_UNCOMPRESS as GUNZIP_WHITELIST                         } from "../modules/nf-core/pigz/uncompress/main"
include { PIGZ_COMPRESS as PIGZ_COMPRESS_10X                          } from "../modules/nf-core/pigz/compress/main"
include { PIGZ_COMPRESS as PIGZ_COMPRESS_ARGENTAG                     } from "../modules/nf-core/pigz/compress/main"
include { NANOCOMP as NANOCOMP_FASTQ                                  } from "../modules/nf-core/nanocomp/main"
include { MULTIQC as MULTIQC_RAWQC                                    } from "../modules/nf-core/multiqc/main"
include { MULTIQC as MULTIQC_FINALQC                                  } from "../modules/nf-core/multiqc/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                                 } from "../modules/nf-core/custom/dumpsoftwareversions/main"
include { CAT_CAT                                                     } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_PREEXTRACT                               } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_BARCODE                                  } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_ARGENTAG                                 } from "../modules/nf-core/cat/cat/main"
include { CAT_FASTQ                                                   } from "../modules/nf-core/cat/fastq/main"
include { CAT_FASTQ as CAT_FASTQ_PARSE                                } from "../modules/nf-core/cat/fastq/main"
include { paramsSummaryMap                                            } from "plugin/nf-schema"

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/subworkflows
 */
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM         } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM        } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT     } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { paramsSummaryMultiqc                                        } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML                                      } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText                                      } from "../subworkflows/local/utils_nfcore_scnanoseq_pipeline"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCNANOSEQ {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    /// Initialize variables
    def GENOME_QUANT_OPTS = ['isoquant']
    def TRANSCRIPT_QUANT_OPTS = ['oarfish']

    def blaze_whitelist = getWhitelistFile()

    def quantifiers = params.quantifier.split(',').collect { it.trim() }

    // Now GENOME_QUANT_OPTS and TRANSCRIPT_QUANT_OPTS should be visible
    def genome_quants = quantifiers.findAll { quantifier ->
        quantifier in GENOME_QUANT_OPTS
    }

    def transcript_quants = quantifiers.findAll { quantifier ->
        quantifier in TRANSCRIPT_QUANT_OPTS
    }

    log.info("DEBUG: params.quantifier = ${params.quantifier}")
    log.info("DEBUG: quantifiers = ${quantifiers}")
    log.info("DEBUG: genome_quants = ${genome_quants}")
    log.info("DEBUG: transcript_quants = ${transcript_quants}")

    //// Config file declarations:
    ch_multiqc_config = channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)


    ch_versions = channel.empty()
    ch_multiqc_report = channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_samplesheet
        .branch { meta, fastq ->
            single: fastq.size() == 1
            return [meta, fastq.flatten()]
            multiple: fastq.size() > 1
            return [meta, fastq.flatten()]
        }
        .set { ch_fastqs }

    //
    // MODULE: Combine fastqs from the same sample
    //
    CAT_FASTQ(ch_fastqs.multiple).reads.mix(ch_fastqs.single).set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot, ToulligQC and FastQC - pre-trim QC
    //

    ch_fastqc_multiqc_pretrim = channel.empty()
    if (!params.skip_qc) {

        FASTQC_NANOPLOT_PRE_TRIM(ch_cat_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc, params.skip_nanoq)

        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_PRE_TRIM.out.nanoq_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
        ch_nanostat_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_txt.ifEmpty([])
        ch_nanoq_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.nanoq_stats.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = channel.empty()
    ch_nanocomp_fastq_txt = channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ(
            ch_cat_fastq.map { meta, fastq ->
                //[ [ 'id': 'nanocomp_fastq.' ] , it ]
                [[id: meta.id ?: 'nanocomp_fastq'], fastq]
            }
        )

        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.report_html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.stats_txt

        ch_versions = ch_versions.mix(NANOCOMP_FASTQ.out.versions.first())
    }

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES(
        params.genome_fasta ? file(params.genome_fasta) : [],
        params.transcript_fasta ? file(params.transcript_fasta) : [],
        params.gtf ? file(params.gtf) : [],
    )

    genome_fasta = PREPARE_REFERENCE_FILES.out.prepped_genome_fasta
    genome_fai = PREPARE_REFERENCE_FILES.out.genome_fai
    transcript_fasta = PREPARE_REFERENCE_FILES.out.prepped_transcript_fasta
    transcript_fai = PREPARE_REFERENCE_FILES.out.transcript_fai
    gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf

    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)

    //
    // MODULE: Generate bed file from input gtf for rseqc
    //

    // come back to this once intron work is finished (likely input will be fine)
    ch_pred = channel.empty()
    ch_rseqc_bed = channel.empty()
    if (!params.skip_qc) {
        UCSC_GTFTOGENEPRED(gtf)
        ch_pred = UCSC_GTFTOGENEPRED.out.genepred
        ch_refflat = UCSC_GTFTOGENEPRED.out.refflat
        ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

        UCSC_GENEPREDTOBED(ch_pred)
        ch_rseqc_bed = UCSC_GENEPREDTOBED.out.bed
        ch_versions = ch_versions.mix(UCSC_GENEPREDTOBED.out.versions)
    }

    //
    // MODULE: Unzip fastq
    //
    GUNZIP_FASTQ(ch_cat_fastq)
    ch_unzipped_fastqs = GUNZIP_FASTQ.out.file
    ch_versions = ch_versions.mix(GUNZIP_FASTQ.out.versions)

    //
    // MODULE: Trim and filter reads
    //
    ch_fastqc_multiqc_postrim = channel.empty()
    ch_trimmed_reads_combined = channel.empty()

    // NOTE: this block of code, doesn't have much sense as it splits and later concatenates without any filtering?
    if (!params.skip_trimming) {
        //
        // MODULE: Split fastq
        //
        ch_fastqs = ch_unzipped_fastqs

        if (params.split_amount > 0) {
            SPLIT_FILE(ch_unzipped_fastqs, '.fastq', params.split_amount)

            // Temporarily change the meta object so that the id is present on the
            // fastq to prevent duplicated names
            SPLIT_FILE.out.split_files
                .transpose()
                .set { ch_fastqs }

            ch_versions = ch_versions.mix(SPLIT_FILE.out.versions)
        }

        ch_trimmed_reads = ch_fastqs
        if (!params.skip_trimming) {
            // This if here doesnt have sense, as if skip_trimming was true, it won't enter this if entry

            NANOFILT(ch_fastqs)
            ch_trimmed_reads = NANOFILT.out.reads
            ch_versions = ch_versions.mix(NANOFILT.out.versions)
        }

        // If the fastqs were split, combine them together
        ch_trimmed_reads_combined = ch_trimmed_reads
        if (params.split_amount > 0) {
            CAT_CAT(ch_trimmed_reads.groupTuple())
            ch_trimmed_reads_combined = CAT_CAT.out.file_out
        }

        //
        // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
        //
        if (!params.skip_qc) {

            //
            // MODULE: Run qc on the post trimmed reads
            //
            FASTQC_NANOPLOT_POST_TRIM(ch_trimmed_reads_combined, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc, params.skip_nanoq)

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_nanostat_posttrim = FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_txt.ifEmpty([])

            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.toulligqc_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.fastqc_version.first().ifEmpty(null))
            ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_TRIM.out.nanoq_version.first().ifEmpty(null))
        }
    }
    else {
        ch_trimmed_reads_combined = ch_unzipped_fastqs
    }

    // DEMULTIPLEXING MODULES:

    if (params.platform == "10X") {
        //
        // MODULE: Unzip whitelist
        //

        // NOTE: Blaze does not support '.gzip'
        ch_blaze_whitelist = blaze_whitelist

        if (blaze_whitelist.endsWith('.gz')) {

            GUNZIP_WHITELIST([[:], blaze_whitelist])

            ch_blaze_whitelist = GUNZIP_WHITELIST.out.file.map { meta, whitelist ->
                [whitelist]
            }

            ch_versions = ch_versions.mix(GUNZIP_WHITELIST.out.versions)
        }

        //
        // MODULE: Generate whitelist
        //

        BLAZE(ch_trimmed_reads_combined, ch_blaze_whitelist)

        ch_putative_bc = BLAZE.out.putative_bc
        ch_gt_whitelist = BLAZE.out.whitelist
        ch_whitelist_bc_count = BLAZE.out.bc_count
        ch_versions = ch_versions.mix(BLAZE.out.versions)

        ch_split_bc_fastqs = ch_trimmed_reads_combined
        ch_split_bc = ch_putative_bc
        if (params.split_amount > 0) {
            SPLIT_FILE_BC_FASTQ(ch_trimmed_reads_combined, '.fastq', params.split_amount)

            SPLIT_FILE_BC_FASTQ.out.split_files
                .transpose()
                .set { ch_split_bc_fastqs }

            ch_versions = ch_versions.mix(SPLIT_FILE_BC_FASTQ.out.versions)

            SPLIT_FILE_BC_CSV(ch_putative_bc, '.csv', (params.split_amount / 4))
            SPLIT_FILE_BC_CSV.out.split_files
                .transpose()
                .set { ch_split_bc }
        }


        //
        // MODULE: Extract barcodes
        //

        PREEXTRACT_FASTQ(ch_split_bc_fastqs.join(ch_split_bc), params.barcode_format)
        ch_barcode_info = PREEXTRACT_FASTQ.out.barcode_info
        ch_preextract_fastq = PREEXTRACT_FASTQ.out.extracted_fastq

        //
        // MODULE: Correct Barcodes
        //

        CORRECT_BARCODES(
            ch_barcode_info.combine(ch_gt_whitelist, by: 0).combine(ch_whitelist_bc_count, by: 0)
        )
        ch_corrected_bc_file = CORRECT_BARCODES.out.corrected_bc_info
        ch_versions = ch_versions.mix(CORRECT_BARCODES.out.versions)

        ch_extracted_fastq = ch_preextract_fastq
        ch_corrected_bc_info = ch_corrected_bc_file

        if (params.split_amount > 0) {
            //
            // MODULE: Cat Preextract
            //
            CAT_CAT_PREEXTRACT(ch_preextract_fastq.groupTuple())
            ch_cat_preextract_fastq = CAT_CAT_PREEXTRACT.out.file_out

            //
            // MODULE: Cat barcode file
            //
            CAT_CAT_BARCODE(ch_corrected_bc_file.groupTuple())
            ch_corrected_bc_info = CAT_CAT_BARCODE.out.file_out

            //
            // MODULE: Zip the reads
            //
            PIGZ_COMPRESS_10X(ch_cat_preextract_fastq)
            ch_extracted_fastq = PIGZ_COMPRESS_10X.out.archive
            ch_versions = ch_versions.mix(PIGZ_COMPRESS_10X.out.versions)
        }
    }
    else if (params.platform == "Parse") {

        // Split fastq files:
        if (params.split_amount > 0) {
            SPLIT_FILE_PARSE(ch_trimmed_reads_combined, '.fq', params.split_amount)

            // Temporarily change the meta object so that the id is present on the
            // fastq to prevent duplicated names
            SPLIT_FILE_PARSE.out.split_files
                .transpose()
                .set { ch_split_parse }

            ch_versions = ch_versions.mix(SPLIT_FILE_PARSE.out.versions)
        }

        // Generate artificial PE reads:
        ch_pe_parse = GENERATE_PE(ch_split_parse)
        ch_versions = ch_versions.mix(ch_pe_parse.versions)

        // Concatenate PE reads:
        ch_to_concat = ch_pe_parse.out
            .groupTuple(by: 0)
            .map { meta, reads ->
                def read_list = reads.flatten()
                [meta, read_list]
            }

        ch_concatenated = CAT_FASTQ_PARSE(ch_to_concat)
        ch_versions = ch_versions.mix(ch_concatenated.versions)

        // Run spipe to generate the corrected fastq:
        spipe_params_file = file("${baseDir}/assets/parfile_parse.txt")
        spipe_mock_genome_dir = "${baseDir}/assets/mock_genome_parse"

        ch_cat_preextract_fastq = SPLITPIPE_PRE(ch_concatenated.reads, spipe_mock_genome_dir, spipe_params_file)
        ch_versions = ch_versions.mix(ch_cat_preextract_fastq.versions)

        ch_extracted_fastq = ch_cat_preextract_fastq.out
        ch_corrected_bc_info = channel.empty()
    }
    else if (params.platform == "Argentag") {

        // Split fastq files:
        if (params.split_amount > 0) {
            SPLIT_FILE_ARGENTAG(ch_trimmed_reads_combined, '.fq', params.split_amount)

            // Temporarily change the meta object so that the id is present on the
            // fastq to prevent duplicated names
            SPLIT_FILE_ARGENTAG.out.split_files
                .transpose()
                .set { ch_split_argentag }

            ch_versions = ch_versions.mix(SPLIT_FILE_ARGENTAG.out.versions)
        }

        // Run chimera splitting:
        ch_chimera_splitted = ARGENTAG_SPLIT(ch_split_argentag)
        ch_versions = ch_versions.mix(ch_chimera_splitted.versions)

        // Concatenate tagged reads per sample:
        ch_to_concat = ch_chimera_splitted.out
            .groupTuple(by: 0)
            .map { meta, reads ->
                def read_list = reads.flatten()
                [meta, read_list]
            }

        ch_concatenated = CAT_CAT_ARGENTAG(ch_to_concat)
        ch_versions = ch_versions.mix(ch_concatenated.versions)

        // Run taggy demux to generate the corrected fastq:
        ch_cat_preextract_fastq = ARGENTAG_TAGGY_DEMUX(ch_concatenated.file_out)
        ch_versions = ch_versions.mix(ch_cat_preextract_fastq.versions)

        // Compress the tagged fastq:
        PIGZ_COMPRESS_ARGENTAG(ch_cat_preextract_fastq.out)
        ch_extracted_fastq = PIGZ_COMPRESS_ARGENTAG.out.archive
        ch_versions = ch_versions.mix(PIGZ_COMPRESS_ARGENTAG.out.versions)

        ch_corrected_bc_info = channel.empty()
    }
    else {
        exit(1, "Single cell platform not recognized. You can choose either 10X, Parse or Argentag.\n")
    }

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = channel.empty()
    ch_read_counts = channel.empty()
    if (!params.skip_qc) {
        FASTQC_NANOPLOT_POST_EXTRACT(ch_extracted_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc, params.skip_nanoq)

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_nanostat_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_txt.ifEmpty([])
        ch_nanoq_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.nanoq_stats.ifEmpty([])
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_NANOPLOT_POST_EXTRACT.out.nanoq_version.first().ifEmpty(null))

        //
        // MODULE: Generate read counts
        //

        ch_pretrim_counts = channel.empty()
        ch_posttrim_counts = channel.empty()
        ch_postextract_counts = channel.empty()
        if (!params.skip_fastqc) {
            ch_pretrim_counts = ch_fastqc_multiqc_pretrim.collect { it[0] }
            ch_posttrim_counts = ch_fastqc_multiqc_postrim.collect { it[0] }
            ch_postextract_counts = ch_fastqc_multiqc_postextract.collect { it[0] }
        }
        else if (!params.skip_nanoplot) {
            ch_pretrim_counts = ch_nanostat_pretrim.collect { it[1] }
            ch_posttrim_counts = ch_nanostat_posttrim.collect { it[1] }
            ch_postextract_counts = ch_nanostat_postextract.collect { it[1] }
        }
    }

    //
    // SUBWORKFLOW: Align Long Read Data
    //

    ch_multiqc_finalqc_files = channel.empty()

    if (genome_quants) {
        PROCESS_LONGREAD_SCRNA_GENOME(
            genome_fasta,
            genome_fai,
            gtf,
            ch_extracted_fastq,
            ch_refflat,
            ch_rseqc_bed,
            ch_corrected_bc_info,
            genome_quants,
            params.dedup_tool,
            true,
            params.fasta_delimiter,
            params.platform,
            params.skip_save_minimap2_index,
            params.skip_qc,
            params.skip_rseqc,
            params.skip_bam_nanocomp,
            params.skip_seurat,
            params.skip_dedup,
        )
    }


    ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_GENOME.out.versions)

    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_flagstat.collect { it[1] }.ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_idxstats.collect { it[1] }.ifEmpty([])
    )
    if (!params.skip_rseqc) {
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_rseqc_read_dist.collect { it[1] }.ifEmpty([])
        )
    }
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_nanocomp_bam_txt.collect { it[1] }.ifEmpty([])
    )

    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.bc_tagged_flagstat.collect { it[1] }.ifEmpty([])
    )

    if (!params.skip_dedup) {
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_flagstat.collect { it[1] }.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_idxstats.collect { it[1] }.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.picard_metrics.collect { it[1] }.ifEmpty([])
        )
    }

    READ_COUNTS_GENOME(
        ch_pretrim_counts.ifEmpty([]),
        ch_postextract_counts.ifEmpty([]),
        ch_corrected_bc_info.collect { it[1] }.ifEmpty([]),
        PROCESS_LONGREAD_SCRNA_GENOME.out.bc_tagged_flagstat.collect { it[1] }.ifEmpty([]),
        PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_flagstat.collect { it[1] }.ifEmpty([]),
    )

    ch_read_counts_genome = READ_COUNTS_GENOME.out.read_counts
    ch_versions = ch_versions.mix(READ_COUNTS_GENOME.out.versions)

    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ch_read_counts_genome.collect().ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.gene_qc_stats.collect().ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        PROCESS_LONGREAD_SCRNA_GENOME.out.transcript_qc_stats.collect().ifEmpty([])
    )


    // oarfish expects deduplicated reads
    if (transcript_quants) {
        PROCESS_LONGREAD_SCRNA_TRANSCRIPT(
            transcript_fasta,
            transcript_fai,
            gtf,
            ch_extracted_fastq,
            ch_refflat,
            ch_rseqc_bed,
            ch_corrected_bc_info,
            transcript_quants,
            params.dedup_tool,
            false,
            params.fasta_delimiter,
            params.platform,
            params.skip_save_minimap2_index,
            params.skip_qc,
            true,
            true,
            params.skip_seurat,
            false,
        )

        ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_flagstat.collect { it[1] }.ifEmpty([])
        )
        if (!params.skip_rseqc) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_rseqc_read_dist.collect { it[1] }.ifEmpty([])
            )
        }
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_nanocomp_bam_txt.collect { it[1] }.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.bc_tagged_flagstat.collect { it[1] }.ifEmpty([])
        )

        if (!params.skip_dedup) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.dedup_flagstat.collect { it[1] }.ifEmpty([])
            )
        }

        READ_COUNTS_TRANSCRIPT(
            ch_pretrim_counts.ifEmpty([]),
            ch_postextract_counts.ifEmpty([]),
            ch_corrected_bc_info.collect { it[1] }.ifEmpty([]),
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_flagstat.collect { it[1] }.ifEmpty([]),
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.dedup_flagstat.collect { it[1] }.ifEmpty([]),
        )

        ch_read_counts_transcript = READ_COUNTS_TRANSCRIPT.out.read_counts
        ch_versions = ch_versions.mix(READ_COUNTS_TRANSCRIPT.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            ch_read_counts_transcript.collect().ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.transcript_qc_stats.collect().ifEmpty([])
        )
    }

    //
    // SOFTWARE_VERSIONS
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (!params.skip_qc && !params.skip_multiqc) {

        //
        // MODULE: MultiQC for raw data
        //

        ch_multiqc_rawqc_files = channel.empty()
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_config)
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.map { meta, txt -> txt }.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanoq_pretrim.collect { it[1] }.ifEmpty([]))

        //MULTIQC_RAWQC (
        //    ch_multiqc_rawqc_files.collect(),
        //    ch_multiqc_config,
        //    ch_multiqc_custom_config.collect().ifEmpty([]),
        //    ch_multiqc_logo.collect().ifEmpty([]),
        //    [],
        //    []
        //)

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_config)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_postextract_counts.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_rawqc_files.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_nanoq_postextract.collect { it[1] }.ifEmpty([]))

        MULTIQC_FINALQC(
            ch_multiqc_finalqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([]),
            [],
            [],
        )
        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_versions = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }

    emit:
    multiqc_report = ch_multiqc_report.toList()
    versions       = ch_versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER PRESETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Helper function - this is OK at top level
def getWhitelistFile() {
    if (params.whitelist) {
        return params.whitelist
    }

    def whitelist_map = ['10X_3v3': "${baseDir}/assets/whitelist/3M-february-2018.zip", '10X_5v2': "${baseDir}/assets/whitelist/737K-august-2016.txt.zip", '10X_3v4': "${baseDir}/assets/whitelist/3M-3pgex-may-2023_TRU.txt.zip", '10X_5v3': "${baseDir}/assets/whitelist/3M-5pgex-jan-2023.txt.zip"]

    def whitelist_path = whitelist_map[params.barcode_format]
    return whitelist_path ? file(whitelist_path) : null
}
