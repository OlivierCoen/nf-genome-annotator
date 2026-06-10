process FETCH_ENA_FASTQ_URLS {

    label 'process_medium'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e08817d4ad945790377716f555bb30a1ff4a11caf213f84d5b14ab9b08ba0cd/data':
        'community.wave.seqera.io/library/requests_tenacity:4bfd6dae88df84fc' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path('ena_ftp_links.txt'),                                                                                                 emit: ftp_urls
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                               topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),                           topic: versions
    tuple val("${task.process}"), val('tenacity'), eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"'),   topic: versions

    script:
    """
    fetch_ena_fastq_urls.py \\
        --accession $accession
    """

}
