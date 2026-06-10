process SRATOOLS_FASTERQDUMP {
    tag "${meta.taxid} :: ${sra.name}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41cb64e962a993ad12270759922b4c4d6cc03dd977f2d32e044552c3a8a91c98/data' :
        'community.wave.seqera.io/library/seqkit_sra-tools_awk_pigz:885a0bc209e70207' }"

    input:
    tuple val(meta), path(sra)
    path ncbi_settings

    output:
    tuple val(meta), path('*.fastq.gz'),                                                                      emit: reads
    tuple val("${task.process}"), val('sratools'), eval("fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+'"),  topic: versions
    tuple val("${task.process}"), val('pigz'),     eval("pigz --version 2>&1 | sed 's/pigz //g'"),            topic: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${sra.name}"
    """
    export NCBI_SETTINGS="\$PWD/${ncbi_settings}"

    fasterq-dump \\
        $args \\
        --split-files \\
        --threads $task.cpus \\
        --outfile ${prefix}.fastq \\
        ${sra}

    # compressing
    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${sra.name}"
    """
    touch ${prefix}.fastq.gz
    """
}
