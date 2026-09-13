nextflow.enable.types = true

process FASTPLONG {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['apptainer', 'singularity'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/764f6aec76feeb8e34d30eda96d43f6925f66e0e76fca225bc2a49ec3f5416d0/data':
        'community.wave.seqera.io/library/fastplong:0.7.0--737a477808c8442c' }"

    input:
        record(
            id: String,
            fastq: Path
        )

    output:
        record(
            id: id,
            fastq: file("${prefix}.fq.gz")
        )

    topic:
        tuple("${task.process}", 'fastplong', eval("fastplong --version | cut -d' ' -f2")) >> 'versions'

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${id}.preprocessed"
    """
    fastplong \\
        $args \\
        --in $fastq \\
        --out ${prefix}.fq.gz \\
        --thread ${task.cpus}
    """
}
