nextflow.enable.types = true

process BUSCO_DOWNLOAD {
    tag "${lineage}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/963bad66c10646cf0adb1967cc462ad04d02789ddbfae4fbb94182291dbddf8c/data'
        : 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'}"

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
