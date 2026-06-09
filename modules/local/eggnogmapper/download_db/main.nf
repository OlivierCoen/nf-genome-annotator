process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    storeDir "${workflow.projectDir}/.nextflow/cache/eggnogmapper"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b514ff051f837a2f06ab62a9a16592b0a2171917f1c7b11bbfef84c632e99f3/data':
        'community.wave.seqera.io/library/aria2_httpx_pigz_tenacity:0cf226d5365c27ae' }"

    input:

    output:
    path "data",          emit: eggnog_data_dir
    tuple val("${task.process}"), val('aria2'), eval("aria2c -v | head -1 | sed 's/aria2 version //g'"), topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),           topic: versions
    tuple val("${task.process}"), val('httpx'), eval('python3 -c "import httpx; print(httpx.__version__)"'), topic: versions

    script:
    """
    mkdir data

    download_eggnog_data.modified.py \\
        --db diamond \\
        --out data
    """
}
