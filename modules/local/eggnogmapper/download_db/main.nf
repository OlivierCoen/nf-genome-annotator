process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    storeDir "${workflow.projectDir}/.nextflow/cache/eggnogmapper"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d12124094d33e21ac770e9447e8e3be6f9208e8bc49f84af236996eb61b243dc/data':
        'community.wave.seqera.io/library/eggnog-mapper:2.1.13--c99d97a9121734e6' }"

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
