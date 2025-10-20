process EARLGREY_EARLGREY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/16c8c9a119c6aaf7d52786859fe29f57480848c13a5a34d3127e392c1b366570/data':
        'community.wave.seqera.io/library/earlgrey:6.3.3--4a2200a6b48c86ec' }"

    input:
    tuple val(meta), path(fasta), path(lib)

    output:
    tuple val(meta), path("${prefix}.masked"),                                                                        emit: masked
    tuple val("${task.process}"), val('earlgrey'), eval("earlGrey | grep version | sed 's/earlGrey version //g'"),    topic: versions

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    ${fasta_gz_cmd}
    RepeatMasker \\
        $lib_arg \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        ${args} \\
        ${out_fasta}
    """

}
