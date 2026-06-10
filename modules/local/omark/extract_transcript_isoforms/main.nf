process OMARK_EXTRACT_TRANSCRIPT_ISOFORMS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '':
        '' }"

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
