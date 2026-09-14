nextflow.enable.types = true

process STAR_ALIGN {
    tag "${id} :: ${read_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db9b7883c4db989e21dc1d702f4e502219e7dc28ec9fb7322e9a3ad84b40e94f/data' :
        'community.wave.seqera.io/library/star_gzip:07738e810fc0bb60' }"

    input:
        record(
            id: String,
            reads: Set<Path>,
            read_id: String,
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
    def prefix = task.ext.prefix ?: "$read_id"
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
