nextflow.enable.types = true

process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    storeDir "${workflow.projectDir}/.nextflow/cache/eggnogmapper_db"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b514ff051f837a2f06ab62a9a16592b0a2171917f1c7b11bbfef84c632e99f3/data':
        'community.wave.seqera.io/library/aria2_httpx_pigz_tenacity:0cf226d5365c27ae' }"

    output:
        file("data", type: 'dir')

    topic:
        tuple("${task.process}", 'eggnog-mapper', eval('emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    mkdir data
    
    download_eggnog_data.py \\
        --data_dir data \\
        ${args}
    """
}
