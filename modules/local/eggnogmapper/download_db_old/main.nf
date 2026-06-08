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
    path "versions.yml",  emit: versions

    script:
    """
    mkdir data

    download_eggnog_data.modified.py \\
        --db diamond \\
        --db-version 7.0.0 \\
        --out data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(aria2c -v | head -1 | sed 's/aria2 version //g')
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """
}
