process TSEBRA_TSEBRA {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/660fcbc43ec4c60cdd3f2ec61f03865b63d8429188b5584b047da06eff4313b2/data':
        'community.wave.seqera.io/library/tsebra:1.1.2.5--8417f53cddae9ef5' }"

    input:
    tuple val(meta), path(gtfs), path(hints_files)
    path keep_gtfs
    path config

    output:
    tuple val(meta), path("*.gtf"), emit: merged_gtf
    tuple val(meta), path("*.tsv"), emit: tsebra_scores

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('tsebra'), eval("1.1.2.5"),    topic: versions

    script:
    def args        = task.ext.args                                     ?: ''
    def prefix      = task.ext.prefix                                   ?: "${meta.id}"
    def gtf_arg     = '-g ' + gtfs.collect { "$it" }.join(',')
    def hints_arg   = '-e ' + hints_files.collect { "$it" }.join(',')
    def keep_arg    = keep_gtfs                                         ? ( '-k ' + keep_gtfs.collect { "$it" }.join(',') ) : ''
    def config_arg  = config                                            ? "-c $config"                                      : ''

    """
    tsebra.py \\
        $gtf_arg \\
        $hints_arg \\
        $keep_arg \\
        $config_arg \\
        $args \\
        -o ${prefix}.gtf \\
        -s ${prefix}.tsv
    """

    stub:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def VERSION     = '1.1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.gtf
    touch ${prefix}.tsv
    """
}
