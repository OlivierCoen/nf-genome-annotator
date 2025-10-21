process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    storeDir "${workflow.projectDir}/.nextflow/cache/eggnogmapper"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:

    output:
    path "data",          emit: eggnog_data_dir
    path "versions.yml",  emit: versions

    script:
    """
    mkdir data

    download_eggnog_data.modified.py \\
        -y \\
        --data_dir data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS
    """
}
