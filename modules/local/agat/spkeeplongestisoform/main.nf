process AGAT_SPKEEPLONGESTISOFORM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abfa03eb1d5ee9a8f9aa056751126647577ee62ffac6a6ab84ca7a2184007380/data':
        'community.wave.seqera.io/library/agat:1.7.0--9487e22276dbaaca' }"

    input:
    tuple val(meta), path(gxf)

    output:
    tuple val(meta), path("*.longest_transcript_isoforms.gff"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat_sp_keep_longest_isoform.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff ${gxf} \\
        --out "${prefix}.longest_transcript_isoforms.gff" \\
        ${args}
    """
}
