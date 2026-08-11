nextflow.enable.types = true

process MMSEQS_CONCATDBS {
    tag "${mmseqs_dbs.join(' ')}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edfecaaca16ca7fb7b6428dce0ed9c737549b38146360c98fdabf74e6c4cac68/data'
        : 'community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d'}"

    input:
        record(
            id: String,
            mmseqs_dbs: Iterable<Path>
        )

    stage:
        stageAs mmseqs_dbs, "mmseqs_dbs/*"

    output:
        record(
            id: id,
            mmseqs_db: file("mmseqs_db")
        )

    topic:
        tuple("${task.process}", 'mmseqs', eval('mmseqs version')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def database_arg = mmseqs_dbs.join(' ')
    """
    mmseqs concatdbs \\
        ${database_arg} \\
        mmseqs_db \\
        --threads 1 \\
        ${args}
    """
}
