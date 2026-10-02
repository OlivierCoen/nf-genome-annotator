nextflow.enable.types = true

process MINIMAP2_INDEX {
    label 'process_low'

    // Note: the versions here need to match the versions used in minimap2/align
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['apptainer', 'singularity'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/8412af28f9f12284d4c01d149e84d58e65af80ecf158991c6e45563769267672/data' :
        'community.wave.seqera.io/library/minimap2_samtools:3a608f92b72cacce' }"

    input:
        record(
            id: String,
            fasta: Path
        )

    output:
        record(
            id: id,
            index: file("*.mmi")
        )

    topic:
        tuple("${task.process}", "minimap2", eval("minimap2 --version")) >> 'versions'
        
    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}