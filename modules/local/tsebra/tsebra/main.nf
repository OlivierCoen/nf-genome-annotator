nextflow.enable.types = true

process TSEBRA_TSEBRA {
    tag "$id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
     container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/660fcbc43ec4c60cdd3f2ec61f03865b63d8429188b5584b047da06eff4313b2/data':
        'community.wave.seqera.io/library/tsebra:1.1.2.5--8417f53cddae9ef5' }"

    input:
        record(
            id: String,
            gtfs: Iterable<Path>,
            hintsfiles: Iterable<Path>
        )

    output:
        record(
            id: id,
            merged_gtf: file("*.gtf")
        )

    topic:
        tuple('tsebra', id, file("*.tsv")) >> 'additional_results'
        // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
        tuple("${task.process}", 'tsebra', "1.1.2.5")          >> 'versions'

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "$id"
    def gtf_arg     = '-g ' + gtfs.collect { "$it" }.join(',')
    def hints_arg   = '-e ' + hintsfiles.collect { "$it" }.join(',')

    """
    tsebra.py \\
        $gtf_arg \\
        $hints_arg \\
        $args \\
        -o ${prefix}.gtf \\
        -s ${prefix}.tsv
    """
}
