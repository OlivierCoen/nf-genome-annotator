nextflow.enable.types = true


process FASTP {
    tag "${id} :: ${read_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    ?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013aad5427d824afe472e6607ea47685ff0181f1fb09e52a179e0ec39e43e88/data'
    :         'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde' }"

    input:
        record(
            id: String, 
            read_id: String,
            reads: Iterable<Path>
        )

    output:
        record(
            id: id,
            read_id: read_id,
            reads: files('*.fastp.fastq.gz')
        )

    topic:
        tuple(id, file('*.json'))                                                              >> 'fastp_multiqc'
        tuple('fastp', id, file('*.fastp.log'))                                                >> 'logs'
        tuple("${task.process}", 'fastp', eval('fastp --version 2>&1 | sed -e "s/fastp //g"')) >> 'versions'


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$read_id"
    if ( reads.size() == 1 ) { // single-end
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    } else { // paired-end
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_R1.fastp.fastq.gz \\
            --out2 ${prefix}_R2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    }
}
