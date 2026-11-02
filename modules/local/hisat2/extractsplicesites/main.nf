nextflow.enable.types = true

process HISAT2_EXTRACTSPLICESITES {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a16b8102041e76f25358477471bafe813828313f7a8787066f390354cd7b8b7c/data'
            : 'community.wave.seqera.io/library/hisat2:2.2.3--2616fa83d3b9d8f8'}"

    input:
        record(
            id: String,
            gtf: Path
        )

    output:
        record(
            id: id,
            splice_sites: file("*.splice_sites.txt")
        )

    topic:
        tuple("${task.process}", 'hisat2', eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py \\
        $args \\
        $gtf \\
        > ${id}.splice_sites.txt
    """
}
