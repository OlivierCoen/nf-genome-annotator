process AGAT_SPMANAGEIDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    path config

    output:
    tuple val(meta), path("*_ids_fixed.gff"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat_sp_manage_IDs.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    """
    agat_sp_manage_IDs.pl \\
        --gff $gff \\
        --prefix $prefix \\
        ${config_param} \\
       ${args} \\
        --output ${prefix}_ids_fixed.gff
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ids_fixed.gff
    """
}
