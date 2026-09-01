nextflow.enable.types = true

process BUSCO_LISTDATASETS {

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f67e816ab2f7ccc9cb2d40874dea1e2e1a8e88ef6a44750b66c0ee55fe8de6c/data'
            : 'community.wave.seqera.io/library/busco:6.1.0--0e40710a525d8d44'}"

    input:
        _s: List // only used to enable caching

    output:
        datasets: Path = file("busco_datasets.yaml")

    topic:
        tuple( "${task.process}", 'busco', eval('busco --version | sed "s/^BUSCO //"') ) >> 'versions'

    script:
    """
    busco --list-datasets > busco_datasets.yaml
    """
}
