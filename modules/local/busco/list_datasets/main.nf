process BUSCO_LISTDATASETS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4137d65ab5b90d2ae4fa9d3e0e8294ddccc287e53ca653bb3c63b8fdb03e882f/data'
        : 'community.wave.seqera.io/library/busco:6.0.0--a9a1426105f81165'}"

    input:

    output:
    path("busco_datasets.yaml"), emit: datasets
    tuple val("${task.process}"), val('busco'), eval('busco --version | sed "s/^BUSCO //"'),    topic: versions


    script:
    """
    busco --list-datasets > busco_datasets.yaml
    """
}
