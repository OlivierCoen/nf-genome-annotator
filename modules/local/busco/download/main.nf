nextflow.enable.types = true

process BUSCO_DOWNLOAD {
    tag "${lineage}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f67e816ab2f7ccc9cb2d40874dea1e2e1a8e88ef6a44750b66c0ee55fe8de6c/data'
            : 'community.wave.seqera.io/library/busco:6.1.0--0e40710a525d8d44'}"

    input:
        lineage: String

    output:
        record(
            busco_lineage: lineage,
            busco_download_path: file("busco_downloads", type: 'dir')
        )

    topic:
        tuple("${task.process}", 'busco', eval("busco --version 2> /dev/null | sed 's/BUSCO //g'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        --download ${lineage} \\
        ${args}
    """
}
