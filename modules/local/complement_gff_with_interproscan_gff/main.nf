nextflow.enable.types = true

process COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF {

    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/456dbad6a07daee9c47503388d128762d334d96e0574ea149bb3a2abcc0337ee/data':
        'community.wave.seqera.io/library/pandas_polars_pyarrow_python:603d6d02549d908e' }"

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
