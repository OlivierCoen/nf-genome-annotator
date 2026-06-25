process OMARK_EXTRACT_TRANSCRIPT_ISOFORMS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6bd87171e0392bac597e053166816c63b0ba7917cbb8985b465f8144c7e7013c/data':
        'community.wave.seqera.io/library/polars:1.42.0--968231e32f4804f6' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.transcript_isoforms"), emit: isoforms
    tuple val("${task.process}"), val('python'),eval("python3 --version | sed 's/Python //'"),            topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    extract_transcript_isoforms_per_gene.py \\
        --gff $gff3 \\
        --out ${meta.id}.transcript_isoforms
    """

}
