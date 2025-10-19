process REPEATMODELER_REPEATMODELER {
    tag "$meta.id"
    label 'process_high'

    errorStrategy = {
        if (task.exitStatus == 100) {
            log.warn("RepeatModeler did not find families for database ${meta.id}")
            return 'ignore'
        } else {
            log.error("RepeatModeler failed with exit status ${task.exitStatus} for ${meta.id}")
            if (task.ext.args.contains('-LTRStruct')) {
                log.info("RepeatModeler on ${meta.id}: retrying without -LTRStruct argument")
                task.ext.args = task.ext.args.replace('-LTRStruct', ' ')
                return 'retry'
            } else  {
                return 'ignore'
            }
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b9637a9d4d72993f80c5014ece53d39161871c94845f420a5313a31a7ae5d2a/data':
        'community.wave.seqera.io/library/repeatmodeler:2.0.7--136d59ab97ab30de' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.fa"),  emit: fasta, optional: true
    tuple val(meta), path("*.stk"), emit: stk, optional: true
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('repeatmodeler'), eval("RepeatModeler --version | sed 's/RepeatModeler version //'"), topic: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def db_name = file(db[0]).getBaseName()
    """
    RepeatModeler \\
        -database $db_name \\
        $args \\
        -threads $task.cpus

    mv ${db_name}-families.stk  ${prefix}.stk || echo "Could not find stk file"
    mv ${db_name}-rmod.log      ${prefix}.log || echo "Could not find log file"
    mv ${db_name}-families.fa   ${prefix}.fa  || (echo "Could not find families fasta file" && exit 100)
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.stk
    touch ${prefix}.log
    """
}
