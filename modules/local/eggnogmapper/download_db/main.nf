process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:

    output:
    path "data", emit: eggnog_data_dir
    tuple val("${task.process}"), val('eggnog-mapper'), eval('emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//"'), topic: versions

    script:
    """
    mkdir data

    download_eggnog_data.modified.py \\
        -y \\
        --data_dir data
    """
}
