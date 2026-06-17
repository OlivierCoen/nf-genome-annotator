process AGAT_SPFUNCTIONALSTATISTCS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abfa03eb1d5ee9a8f9aa056751126647577ee62ffac6a6ab84ca7a2184007380/data':
        'community.wave.seqera.io/library/agat:1.7.0--9487e22276dbaaca' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat_sp_functional_statistics.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"), topic: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = meta.final_annotation ? "${meta.id}.final_annotation" : "${gff.simpleName}.intermediate_annotation"
    def genome_size_arg = meta.genome_size ? "--gs ${meta.genome_size}" : ''
    """
    agat_sp_functional_statistics.pl \\
        --gff ${gff} \\
        ${genome_size_arg} \\
        --output ${prefix}.gtf_func_stats.txt \\
        --yaml \\
        ${args}

    # parse yaml file
    parse_gff_stat_file.py \\
        --gff ${prefix}.gtf_stats.txt.yaml \\
        --prefix ${prefix}
    """
}
