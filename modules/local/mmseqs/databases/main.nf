nextflow.enable.types = true

process MMSEQS_DATABASES {
    tag "${database}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edfecaaca16ca7fb7b6428dce0ed9c737549b38146360c98fdabf74e6c4cac68/data'
        : 'community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d'}"

    input:
        database: String

    output:
        record(
            id: database,
            db: file("mmseqs_db/", type: 'dir')
        )


    topic:
        tuple("${task.process}", 'mmseqs', eval('mmseqs version')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    mkdir mmseqs_db/

    mmseqs databases \\
        ${database} \\
        mmseqs_db/${database} \\
        tmp/ \\
        --threads ${task.cpus} \\
        ${args}

    """
}
