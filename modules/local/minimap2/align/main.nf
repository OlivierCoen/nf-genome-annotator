nextflow.enable.types = true

process MINIMAP2_ALIGN {
    tag "${read_id} on ${fasta.baseName}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['apptainer', 'singularity'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/8412af28f9f12284d4c01d149e84d58e65af80ecf158991c6e45563769267672/data' :
        'community.wave.seqera.io/library/minimap2_samtools:3a608f92b72cacce' }"

    input:
        record(
            id: String,
            read_id: String,
            fastq: Path,
            fasta: Path
        )

    output:
        record(
            id: id,
            read_id: read_id,
            bam: file("*.bam"),
            bai: file("*.bai")
        )

    topic:
        tuple("${task.process}", 'minimap2', eval('minimap2 --version'))                                  >> 'versions'
        tuple("${task.process}", 'samtools', eval('samtools --version | head -1 | sed "s/samtools //g"')) >> 'versions'

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "$read_id"
    // samtools flags: 
    // -F 4 : exclude read unmapped
    // -F 256: exclude if not primary alignment
    """
    minimap2 \\
        -ax splice \\
        -uf \\
        -t $task.cpus \\
        $fasta \\
        $fastq \\
        $args \\
        | samtools view \\
            -bS \\
            -F 4 -F 256 - \\
            | samtools sort \\
                -@ ${task.cpus-1} \\
                -o ${prefix}.bam##idx##${prefix}.bam.bai \\
                --write-index \\
                ${args2}
    """

    stub:
    def prefix = task.ext.prefix ?: "$read_id"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
