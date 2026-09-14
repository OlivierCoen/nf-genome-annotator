nextflow.enable.types = true

process STAR_ALIGN {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b03f5c57d28f4975bbbda74a56202f192c69744e3f4533463cc2dfc1bde2bba/data' :
        'community.wave.seqera.io/library/star:2.7.11b--5300af0cf0d14492' }"

    input:
        record(
            id: String,
            reads: Iterable<Path>,
            index: Path
        )

    stage:
        stageAs reads, "input*/*"

    output:
        record(
            id: id,
            bam: file('*.Aligned.out.bam')
        )

    topic:
        tuple(id, file('*Log.final.out'))         >> 'star_multiqc'
        tuple('star', id, file('*Log.final.out')) >> 'logs'
        tuple("${task.process}", 'star', eval('STAR --version | sed "s/STAR_//g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    def read_file_command_arg = reads[0].extension == 'gz' ? "--readFilesCommand zcat": ''
    """
    # Note: '--outSAMstrandField intronMotif' is required for BRAKER
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads.join(",")} \\
        $read_file_command_arg \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --outSAMstrandField intronMotif \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes All \\
        $args
    """
}
