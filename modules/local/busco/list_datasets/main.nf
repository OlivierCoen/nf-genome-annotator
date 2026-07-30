nextflow.enable.types = true

process BUSCO_LISTDATASETS {

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/963bad66c10646cf0adb1967cc462ad04d02789ddbfae4fbb94182291dbddf8c/data'
            : 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'}"

    output:
        datasets: Path = file("busco_datasets.yaml")

    topic:
        tuple( "${task.process}", 'busco', eval('busco --version | sed "s/^BUSCO //"') ) >> 'versions'

    script:
    """
    busco --list-datasets > busco_datasets.yaml
    """
}
