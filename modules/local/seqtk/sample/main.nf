nextflow.enable.types = true

process SEQTK_SAMPLE {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1aa6084536813ab32b6373b04407a898ab5de44692763b59ee0705e44974e0de/data' :
        'community.wave.seqera.io/library/seqtk_pigz:aa99a20f06d8e9a8' }"

    input:
        record(
            id: String,
            reads: Set<Path>
        )
        sample_size: Number

    output:
        record(
            id: id,
            reads: files("*.sampled.fastq.gz")
        )

    topic:
        tuple("${task.process}", 'seqtk', eval("seqtk 2>&1 | awk 'NR==3' | sed 's/Version: //g'")) >> 'versions'
        tuple("${task.process}", 'pigz', eval("pigz --version 2>&1 | sed 's/pigz //g'"))           >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    if (!(args ==~ /.*\ -s\ ?[0-9]+.*/)) {
        args += " -s100"
    }
    // ensuring that sample_size is an integer when it is >= 1
    sample_size = sample_size >= 1 ? sample_size.toInteger() : sample_size
    """
    for file in ${reads.join(' ')}
    do
        FILE_STEM=\$(echo "\$file" | sed -E 's/\\.(fastq|fq)(\\.gz)?\$//')
        seqtk \\
            sample \\
            $args \\
            \$file \\
            $sample_size \\
            | pigz --no-name -p ${task.cpus} > \${FILE_STEM}.sampled.fastq.gz
    done
    """
}
