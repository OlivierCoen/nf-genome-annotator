process INTERPROSCAN_DOWNLOADDB {

    label 'process_medium'
    tag "${meta.id}"

    storeDir "${workflow.projectDir}/.nextflow/cache/interproscan"

    errorStrategy {
        if (task.exitStatus == 100) {
            log.warn("md5 checksum failed for Interproscan URL ${db_url}. Please delete the local file and relaunch the pipeline.")
            return 'terminate'
        }
    }

    input:
    tuple val(meta), val(db_url)

    output:
    path("*/data"),            emit: db
    path "versions.yml",       emit: versions

    script:
    def filename = db_url.tokenize("/")[-1]
    """
    aria2c \\
        -s ${task.cpus} \\
        -x ${task.cpus} \\
        -c \\
        --max-tries=10 \\
        --retry-wait=30 \\
        --timeout=60 \\
        "${db_url}"

    echo "Checking md5"
    aria2c -c "${db_url}.md5"
    md5sum -c --status ${filename}.md5 && echo "ok" || exit 100

    echo "Extracting archive"
    tar -pxzf ${filename}

    echo "Deleting archive"
    rm ${filename} ${filename}.md5


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(aria2c -v | head -1 | sed 's/aria2 version //g')
    END_VERSIONS
    """


}
