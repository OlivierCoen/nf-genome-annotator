nextflow.enable.types = true

process STAR_GENOMEGENERATE {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b03f5c57d28f4975bbbda74a56202f192c69744e3f4533463cc2dfc1bde2bba/data' :
        'community.wave.seqera.io/library/star:2.7.11b--5300af0cf0d14492' }"

    input:
        record(
            id: String,
            fasta: Path,
            gtf: Path
        )
        ignore_existing_gtf: Boolean

    output:
        record(
            id: id,
            index: file('star', type: 'dir')
        )

    topic:
        tuple("${task.process}", 'star', eval('STAR --version | sed "s/STAR_//g"')) >> 'versions'

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def gtf_arg     = ignore_existing_gtf ? "" : gtf ? "--sjdbGTFfile $gtf" : ''
    """
    NUM_BASES=\$(grep -v '^>' $fasta | tr -d '\n' | wc -c)

    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $gtf_arg \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases \$NUM_BASES \\
        $memory \\
        $args
    """

}
