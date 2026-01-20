process AGAT_CONVERTSPGFF2GTF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d0194019d3fcefea0875ce4703df12dd2244cf9b8932b935197e9063cbc118ae/data' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.agat.gtf"), emit: output_gtf
    tuple val(meta), path("agat.log")     , emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_convert_sp_gff2gtf.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_convert_sp_gff2gtf.pl \\
        --gff ${gff} \\
        --output ${prefix}.agat.gtf \\
        ${args} \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.agat.gtf
    touch ${gff}.agat.log
    """
}
