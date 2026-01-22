process AGAT_SPMERGEANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d0194019d3fcefea0875ce4703df12dd2244cf9b8932b935197e9063cbc118ae/data':
        'community.wave.seqera.io/library/agat:1.6.1--d39db4f54af12afb' }"

    input:
    tuple val(meta), path(gff)
    path config

    output:
    tuple val(meta), path("*_prepared.gff"), emit: gff
    tuple val(meta), path("agat.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_sp_merge_annotations.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    """
    agat_sp_merge_annotations.pl \\
        --gff ${gff} \\
        ${config_param} \\
        --output ${prefix}_prepared.gff \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_prepared.gff
    """
}
