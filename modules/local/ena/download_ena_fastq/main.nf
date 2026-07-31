nextflow.enable.types = true

process DOWNLOAD_ENA_FASTQ {

    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
        record(
            id: String, 
            ftp_urls: Path
        )

    output:
        record(
            id: id, 
            reads: files('*.fastq.gz')
        )
        
    topic:
        tuple( "${task.process}", 'aria2', eval('aria2c --version | head -1 | cut -d" " -f3') ) >> 'versions'

    script:
    """
    for url in \$(cat ${ftp_urls}); do
        echo "Downloading \${url}"
        aria2c \\
            -x ${task.cpus} \\
            -s ${task.cpus} \\
            -o \${url##*/} \\
            \${url}
    done
    """

}
