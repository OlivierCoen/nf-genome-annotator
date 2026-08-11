nextflow.enable.types = true

process SAMTOOLS_SORT {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
        record(id: String, bam: Path)

    output:
        record(
            id: id,
            bam: file("*.bam"),
            bai: file("*.bai")
        )

    topic:
        tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${id}.sorted"
    """
    samtools sort \\
        ${args} \\
        -T ${prefix} \\
        --threads ${task.cpus} \\
        -o ${prefix}.bam##idx##${prefix}.bam.bai \\
        --write-index \\
        ${bam}
    """
}
