process COMPLEMENT_GFF3_WITH_INTERPROSCAN {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '':
        '' }"

    input:
    tuple val(meta), path(gff3), path(iprscan)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    tuple val("${task.process}"), val('python'),eval("python3 --version | sed 's/Python //'"),            topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.complemented"
    """
    complement_gff3_with_interproscan.py \\
        --annot $gff3 \\
        --iprscan $iprscan \\
        --out ${prefix}.gff3
    """

}
