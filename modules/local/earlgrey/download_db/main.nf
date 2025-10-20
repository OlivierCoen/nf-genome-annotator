process EARLGREY_DOWNLOADDB {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy = {
        if (task.exitStatus == 100) {
            log.warn("md5 checksum failed for FamDB URL ${db_url}. Please delete the local file and relaunch the pipeline.")
            return 'retry'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/facc3740411bc51eaa97d2e0208a7dd571bf5b650b16173e747cc3b6de2b0c3c/data':
        'community.wave.seqera.io/library/aria2_pigz:6b38092500fd4da6' }"

    input:
    val partition

    output:
    path("*/data"), emit: db
    tuple val("${task.process}"), val('aria2'), eval("aria2c -v | head -1 | sed 's/aria2 version //g'"), topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),           topic: versions

    script:
    def filename = "dfam39_full.${partition}.h5.gz"
    def url = "https://dfam.org/releases/current/families/FamDB/${filename}"
    """
    aria2c \\
        -s ${task.cpus} \\
        -x ${task.cpus} \\
        -c \\
        --max-tries=10 \\
        --retry-wait=30 \\
        --timeout=60 \\
        "${url}"

    echo "Checking md5"
    aria2c -c "${url}.md5"
    md5sum -c --status ${filename}.md5 && echo "Chechsum ok" || exit 100

    echo "Extracting archive"
    pigz -d ${filename}

    echo "Deleting archive"
    rm ${filename} ${filename}.md5
    """

}
