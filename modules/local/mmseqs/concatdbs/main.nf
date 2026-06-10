process MMSEQS_CONCATDBS {
    tag "${databases.join(' ')}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edfecaaca16ca7fb7b6428dce0ed9c737549b38146360c98fdabf74e6c4cac68/data'
        : 'community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d'}"

    input:
    path databases

    output:
    path "all_mmseqs_dbs", emit: db
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    script:
    def args = task.ext.args ?: ''
    def database_arg = databases.join(' ')
    """
    mmseqs concatdbs \\
        ${database_arg} \\
        all_mmseqs_dbs \\
        --threads 1 \\
        ${args}

    """
}
