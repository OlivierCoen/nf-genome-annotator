process OMARK_DOWNLOADDB {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
    val db_url

    output:
    path("omamer_db/*"), emit: db
    tuple val("${task.process}"), val('aria2'), eval("aria2c -v | head -1 | sed 's/aria2 version //g'"), topic: versions

    script:
    """
    echo "Downloading ${db_url}"
    aria2c \\
        -s ${task.cpus} \\
        -x ${task.cpus} \\
        --max-tries=10 \\
        --retry-wait=30 \\
        --timeout=60 \\
        --optimize-concurrent-downloads \\
        --check-integrity \\
        --dir omamer_db \\
        ${db_url}
    """

}
