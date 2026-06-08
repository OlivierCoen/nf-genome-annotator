process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d12124094d33e21ac770e9447e8e3be6f9208e8bc49f84af236996eb61b243dc/data':
        'community.wave.seqera.io/library/eggnog-mapper:2.1.13--c99d97a9121734e6' }"

    input:

    output:
    path "data",          emit: eggnog_data_dir
    tuple val("${task.process}"), val('eggnog-mapper'), eval('emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//"'), topic: versions


    script:
    """
    DOWNLOAD_SCRIPT=\$(dirname \$(which emapper.py))/download_eggnog_data.py

    mkdir data
    
    \$DOWNLOAD_SCRIPT \\
        --data_dir data \\
        -y
    """
}
