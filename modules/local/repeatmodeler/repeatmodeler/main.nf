process REPEATMODELER_REPEATMODELER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.5--pl5321hdfd78af_0':
        'biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.fa"),                                                                                          emit: fasta, optional: true
    tuple val(meta), path("*.stk"),                                                                                         emit: stk, optional: true
    tuple val(meta), path("*.log"),                                                                                         emit: log, optional: true
    tuple val("${task.process}"), val('repeatmodeler'), eval("RepeatModeler --version | sed 's/RepeatModeler version //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def db_name = file(db[0]).getBaseName()
    """
    RepeatModeler \\
        -database $db_name \\
        $args \\
        -threads $task.cpus

    mv ${db_name}-families.fa   ${prefix}.fa || echo "Could not find families fasta file"
    mv ${db_name}-families.stk  ${prefix}.stk || echo "Could not find stk file"
    mv ${db_name}-rmod.log      ${prefix}.log || echo "Could not find log file"
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.stk
    touch ${prefix}.log
    """
}
