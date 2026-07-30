nextflow.enable.types = true

process ORTHODB_GETCLADES {

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ef/efdbe0b1f00aeb72d474fd6dc2bab29103586eb97decbd010092f36d58837cd6/data':
        'community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94' }"

    output:
        clades: Path = file("odb12v2_levels.tab.gz")

    topic:
        tuple( "${task.process}", 'wget', eval("wget -h 2>&1 | head -1 | cut -d' ' -f3 | sed 's/,//g'") ) >> 'versions'

    script:
    """
    wget --tries 10 https://data.orthodb.org/current/download/odb_data_dump/odb12v2_levels.tab.gz
    """

}
