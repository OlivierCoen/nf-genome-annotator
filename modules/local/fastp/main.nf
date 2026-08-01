nextflow.enable.types = true

process FASTP {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    ?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013aad5427d824afe472e6607ea47685ff0181f1fb09e52a179e0ec39e43e88/data'
    :         'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde' }"

    input:
        record(id: String, reads: List<Path>)

    output:
        record(
            id: id,
            reads: file('*.fastp.fastq.gz')
        )

    topic:
        tuple(id, file('*.json'))                                                              >> 'fastp_multiqc'
        tuple('fastp', id, file('*.json'))                                                     >> 'logs'
        tuple("${task.process}", 'fastp', eval('fastp --version 2>&1 | sed -e "s/fastp //g"')) >> 'versions'


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    def single_end = reads.size() == 1 ? true : false
    def out_fq1 = single_end ? "--out1 ${prefix}.fastp.fastq.gz" : "--out1 ${prefix}_R1.fastp.fastq.gz"
    def out_fq2 = single_end ? "" : "--out2 ${prefix}_R2.fastp.fastq.gz"
    if (single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            $out_fq1 \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    } else {
        """
        [ ! -f  ${prefix}_R1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f  ${prefix}_R2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_R2.fastq.gz
        fastp \\
            --in1 ${prefix}_R1.fastq.gz \\
            --in2 ${prefix}_R2.fastq.gz \\
            $out_fq1 \\
            $out_fq2 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    }
}
