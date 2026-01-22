process AGAT_SQREMOVEREDUNDANTENTRIES {
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
    tuple val(meta), path("*_redundant_entries_removes.gff"), emit: gff
    tuple val(meta), path("agat.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_sq_remove_redundant_entries.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    """
    agat_sq_remove_redundant_entries.pl \\
        --gff $gff \\
        ${config_param} \\
        --output ${prefix}_redundant_entries_removes.gff \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_redundant_entries_removes.gff
    """
}
