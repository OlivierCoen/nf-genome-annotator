nextflow.enable.types = true

process OMARK_OMAMERSEARCH {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7082d2bbe6948025958cd10031f313633a326e5b3c92b122af1224620de0e9c/data'
        : 'community.wave.seqera.io/library/omark:0.5.0--a004ed2a54d1ff9a'}"

    input:
        record(id: String, fasta: Path)
        omamer_db: Path

    output:
        record(
            id: id,
            omamer: file("*.omamer")
        )

    topic:
        tuple("${task.process}", 'omark', eval('omark --version')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    omamer search \\
        --db ${omamer_db} \\
        --query ${fasta} \\
        --nthreads ${task.cpus} \\
        --out ${prefix}.omamer \\
        ${args}
    """

   
}
