nextflow.enable.types = true

process INFER_STRAND_FROM_RSEQC {

    tag "$read_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7977dd87c418c76b502f06c515e51fd8eedd1aeee9c76378d0ec0bf1a8e4c3e/data':
        'community.wave.seqera.io/library/python:3.14.7--629fbce4f44dac86' }"

    input:
        record(
            read_id: String,
            report: Path
        )

    output:
        record(
            read_id: read_id,
            strandedness: env('STRANDEDNESS')
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'

    script:
    """
    parse_rseqc_infer_experiment_report.py \\
        --report $report
    STRANDEDNESS=\$(cat strandedness.txt)
    """

}
