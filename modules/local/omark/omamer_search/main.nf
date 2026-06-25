process OMARK_OMAMERSEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cf/cf01ebcc6512d0cc98453f0556811e1e26dfa008b3e01c107971f47f179dce8a/data'
        : 'community.wave.seqera.io/library/omark:0.4.1--22126bb9fb82176e'}"

    input:
    tuple val(meta), path(fasta)
    path omamer_db

    output:
    tuple val(meta), path("${meta.id}.omamer"), emit: omamer
    // TODO: when done on OMArk's side, add dynamic retrieval of version
    // https://github.com/DessimozLab/OMArk/issues/52
    tuple val("${task.process}"), val('omark'), val('0.4.1'), topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    omamer search \\
        --db ${omamer_db} \\
        --query ${fasta} \\
        --nthreads ${task.cpus} \\
        --out ${meta.id}.omamer \\
        ${args}
    """

   
}
