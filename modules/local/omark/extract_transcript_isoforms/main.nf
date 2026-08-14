nextflow.enable.types = true

process OMARK_EXTRACT_TRANSCRIPT_ISOFORMS {

    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6bd87171e0392bac597e053166816c63b0ba7917cbb8985b465f8144c7e7013c/data':
        'community.wave.seqera.io/library/polars:1.42.0--968231e32f4804f6' }"

    input:
        record(id: String, gff: Path)

    output:
        record(
            id: id,
            transcript_isoforms: file("*.transcript_isoforms")
        )

    topic:
        tuple("${task.process}", 'python', eval("python3 --version | sed 's/Python //'")) >> 'versions'
        tuple("${task.process}", 'polars', eval('python3 -c "import polars; print(polars.__version__)"')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "$id"
    """
    extract_transcript_isoforms_per_gene.py \\
        --gff $gff \\
        --out ${prefix}.transcript_isoforms
    """

}
