nextflow.enable.types = true

process INFER_STRAND_FROM_XS {

    tag "$read_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/bea85d60dc407473ffc48213631e18323a4a46b862c81da99b304f7ba418e33a/data':
        'community.wave.seqera.io/library/python_pysam:9d1352b5cf515a6f' }"

    input:
        record(
            read_id: String,
            bam: Path
        )

    output:
        record(
            read_id: read_id,
            strandedness: env('STRANDEDNESS')
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'
        tuple( "${task.process}", 'pysam',  eval('python3 -c "import pysam; print(pysam.__version__)"') ) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    infer_strand_from_xs.py \\
        $bam \\
        $args
    STRANDEDNESS=\$(cat strandedness.txt)
    """

}
