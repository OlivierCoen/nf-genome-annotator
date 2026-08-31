nextflow.enable.types = true

process COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF {

    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c07a5b09c5b7388edfd0104ee8b33ff07ba74cf59ae10552e2f85a4c8ae97dfd/data':
        'community.wave.seqera.io/library/python_polars_pandas_pyarrow:9d37415130a0bb7d' }"

    input:
        record(
            id: String,
            gff: Path,
            interproscan_gff: Path
        )

    output:
        record(
            id: id,
            gff: file("*.gff3")
        )

    topic:
        tuple("${task.process}", 'python', eval("python3 --version | sed 's/Python //'"))                 >> 'versions'
        tuple("${task.process}", 'polars', eval('python3 -c "import polars; print(polars.__version__)"')) >> 'versions'
        tuple("${task.process}", 'pandas', eval('python3 -c "import pandas; print(pandas.__version__)"')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${id}.complemented"
    """
    complement_gff_with_interproscan_gff.py \\
        --annot ${gff} \\
        --iprscan ${interproscan_gff} \\
        --out ${prefix}.gff3
    """

}
