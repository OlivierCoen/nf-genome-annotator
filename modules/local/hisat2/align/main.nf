process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c36472269e8898f63b7b65dd40433462d541f9e75f9401f0bf8488021275d006/data' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.log")                   , topic: hisat2_summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness_arg = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness_arg = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    } 
    
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    if (meta.single_end) {
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        hisat2 \\
            -x \$INDEX \\
            -U $reads \\
            $strandedness_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $rg \\
            $unaligned \\
            $args \\
            | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $strandedness_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $rg \\
            $unaligned \\
            --no-mixed \\
            --no-discordant \\
            $args \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
        """
    }

}
