nextflow.enable.types = true

process HISAT2_ALIGN {
    tag "${id} :: ${read_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca541d2e69c03b3cb2a10da5dc6fbe21d4548f90bbdc08d5e4f8e13f0fea9a75/data' :
            'community.wave.seqera.io/library/hisat2_samtools:add4b555d95c067d' }"

    input:
        record(
            id: String,
            reads: Set<Path>,
            read_id: String,
            index: Path,
            strandedness: String?
        )

    output:
        record(
            id: id,
            read_id: read_id,
            bam: file("*.bam")
        )

    topic:
        tuple(id, file('*.hisat2.summary.log'))           >> 'hisat2_multiqc'
        tuple('hisat2', id, file('*.hisat2.summary.log')) >> 'logs'
        tuple("${task.process}", 'hisat2', eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2')) >> 'versions'
        tuple("${task.process}", 'samtools', eval("samtools --version | sed -n '1s/samtools //p'"))              >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$read_id"
    def single_end = reads.size() == 1

    def strandedness_arg = ''
    if (strandedness == 'forward') {
        strandedness_arg = single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (strandedness == 'reverse') {
        strandedness_arg = single_end ? '--rna-strandness R' : '--rna-strandness RF'
    } else {
        strandedness_arg = ''
    }

    // NOTE: --dta is required for braker (makes HISAT2 write the XS:A:+/- tag)
    // samtools flags: 
    // -F 4 : exclude read unmapped
    // -F 8 : exclude mate unmapped
    // -F 256: exclude if not primary alignment
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    if ( single_end ) {
        """
        # find is not included in the Docker image, so use ls instead
        INDEX=\$(ls -1 hisat2/*.1.ht2 | sed 's/\\.1.ht2.*\$//')

        hisat2 \\
            -x \$INDEX \\
            -U ${reads[0]} \\
            $strandedness_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            --dta \\
            $rg \\
            $args \\
            | samtools view \\
                -bS \\
                -F 4 -F 256 - \\
                | samtools sort \\
                    -@ ${task.cpus-1} \\
                    -o ${prefix}.bam
        """
    } else {
        """
        # find is not included in the Docker image, so use ls instead
        INDEX=\$(ls -1 hisat2/*.1.ht2 | sed 's/\\.1.ht2.*\$//')

        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $strandedness_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            --dta \\
            $rg \\
            --no-mixed \\
            --no-discordant \\
            $args \\
            | samtools view \\
                -bS \\
                -F 4 -F 8 -F 256 - \\
               | samtools sort \\
                   -@ ${task.cpus-1} \\
                   -o ${prefix}.bam
        """
    }

}
