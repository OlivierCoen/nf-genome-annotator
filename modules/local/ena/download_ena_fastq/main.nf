process DOWNLOAD_ENA_FASTQ {

    label 'process_single'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
    tuple val(meta), path(ena_ftp_url_file)

    output:
    tuple val(meta), path('*.fastq.gz'),                                                              emit: fastq
    tuple val("${task.process}"), val('aria2'), eval('aria2c --version | head -1 | cut -d" " -f3'),   topic: versions

    script:
    """
    for url in \$(cat ${ena_ftp_url_file}); do
        echo "Downloading \${url}"
        aria2c \\
            -x ${task.cpus} \\
            -s ${task.cpus} \\
            -o \${url##*/} \\
            \${url}
    done
    """

}
