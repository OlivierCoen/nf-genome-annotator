nextflow.enable.types = true

process INFER_STRAND_FROM_XS {

    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '':
        '' }"

    input:
        record( id: String, bam: Path )

    output:
        record( id: id, strandedness: env('STRANDEDNESS') )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'
        tuple( "${task.process}", 'pysam',  eval('python3 -c "import pysam; print(pysam.__version__)"') ) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    bin/infer_strand_from_xs.py \\
        $bam \\
        $args
    STRANDEDNESS=\$(cat strandedness.txt)
    """

}
