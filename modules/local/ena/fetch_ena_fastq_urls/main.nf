nextflow.enable.types = true

process FETCH_ENA_FASTQ_URLS {

    label 'process_medium'

    tag "$id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e08817d4ad945790377716f555bb30a1ff4a11caf213f84d5b14ab9b08ba0cd/data':
        'community.wave.seqera.io/library/requests_tenacity:4bfd6dae88df84fc' }"

    input:
        id: String

    output:
        record(
            id: id, 
            ftp_urls: file('ena_ftp_links.txt')
        ) 
    
    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )                                               >> 'versions'
        tuple( "${task.process}", 'requests', eval('python3 -c "import requests; print(requests.__version__)"') )                         >> 'versions'
        tuple( "${task.process}", 'tenacity', eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"') ) >> 'versions'
        
    script:
    """
    fetch_ena_fastq_urls.py \\
        --accession $id 
    """

}
