nextflow.enable.types = true

process REPEATMODELER_REPEATMODELER {
    tag "$id"
    label 'process_high'

    errorStrategy {
        log.warn("RepeatModeler failed with exit status ${task.exitStatus} for ${id}")
        if (task.ext.args && task.ext.args.contains('-LTRStruct')) {
            log.info("RepeatModeler on ${id}: retrying without -LTRStruct argument")
            task.ext.args = task.ext.args.replace('-LTRStruct', ' ')
            return 'retry'
        } else  {
            return 'ignore'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4618e2da2755ead8794cee7e6fe99de4ef739a97364e08ba9fa8b2988cdba123/data':
        'community.wave.seqera.io/library/repeatmodeler:2.0.9--7529329ebd736619' }"

    input:
        record(id: String, db: Iterable<Path>)

    output:
        record(id: id, lib: file("*.fa", optional: true))
            
    topic:
        tuple('repeatmodeler', file("*.stk", optional: true)) >> 'additional_results'
        tuple('repeatmodeler', file("*.log", optional: true)) >> 'logs'
        tuple("${task.process}", 'repeatmodeler', eval("RepeatModeler --version | sed 's/RepeatModeler version //'")) >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${id}"
    def db_name = file(db[0]).getBaseName()
    def args_cli = args ? "--args ${args}" : ""
    """
    run_repeatmodeler.sh \\
        --database ${db_name} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${args_cli}

    mv ${db_name}-families.stk  ${prefix}.stk || echo "Could not find stk file"
    mv ${db_name}-rmod.log      ${prefix}.log || echo "Could not find log file"
    mv ${db_name}-families.fa   ${prefix}.fa  || (echo "Could not find families fasta file" && exit 0)
    """

    stub:
    def prefix  = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.stk
    touch ${prefix}.log
    """
}
