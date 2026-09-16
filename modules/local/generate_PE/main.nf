process GENERATE_PE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker://biocorecrg/parse_pe_ont:4.0'

    input:
    tuple val(meta),   path(reads)

    output:
    tuple val(meta), path("*_R{1,2}.fastq.gz")                       , emit: out
    path  "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def sample_name   = reads.name.replaceAll(/\.fq$/, '')

    """

    LR_generate_pairs.py \
        ${args} \
        --out_dir ./ \
        --fastq ${reads} \
        --new_fname ${sample_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_pe: 2.0

    END_VERSIONS

    """


}
