nextflow.enable.types = true

process SAMTOOLS_MERGE {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
        record(
            id: String,
            bams: Set<Path>,
            bais: Set<Path>
        )

    stage:
        stageAs bams, '?/*'
        stageAs bais, '?/*'

    output:
        record(
            id: id,
            bam: file("*.bam")
        )

    topic:
        tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${prefix}.bam \\
        ${bams.join(' ')}
    """
}
