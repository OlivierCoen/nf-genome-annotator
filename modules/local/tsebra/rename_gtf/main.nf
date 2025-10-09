process TSEBRA_RENAME_GTF {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/660fcbc43ec4c60cdd3f2ec61f03865b63d8429188b5584b047da06eff4313b2/data':
        'community.wave.seqera.io/library/tsebra:1.1.2.5--8417f53cddae9ef5' }"

    input:
    tuple val(meta), path(gtf)
   
    output:
    tuple val(meta), path("*_renamed.gtf"), emit: renamed_gtf

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('tsebra'), eval("1.1.2.5"),    topic: versions

    script:
    def args        = task.ext.args      ?: ''
    def prefix      = task.ext.prefix    ?: "${meta.id}"
    """
    rename_gtf.py \\
        --gtf $gtf \\
        --prefix $prefix \\
        ${args} \\
        --out ${gtf.baseName}_renamed.gtf
    """

}
