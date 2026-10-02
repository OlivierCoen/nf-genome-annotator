nextflow.enable.types = true

process BBMAP_BBNORM {
    tag "$read_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
        record(
            read_id: String,
            fastq: Path
        )

    output:
        record(
            read_id: read_id,
            fastq: file("*.fastq.gz")
        )

    topic:
        tuple('bbnorm', id, file('*.log')) >> 'logs'
        tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$read_id"

    input  = meta.single_end ? "in=${fastq.join(',')}" : "in=${fastq[0]} in2=${fastq[1]}"
    output = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.nm.fastq.gz out2=${prefix}_2.nm.fastq.gz"

    memory = '-Xmx3g'
    if ( ! task.memory ) {
        log.info '[BBNorm]: Available memory not known, defaulting to 3 GB. Specify process memory requirements to change this.'
    } else {
        memory = "-Xmx${Math.round(Math.max(1, Math.floor(task.memory.toGiga() * 0.95)))}g"
    }

    // Nextflow forwards the host's \$TMPDIR into Singularity/Apptainer containers without
    // bind-mounting it, so multipass temp files must stay under the task's own (bound) work dir.
    """
    mkdir -p tmp

    bbnorm.sh \\
        $input \\
        $output \\
        $args \\
        threads=$task.cpus \\
        $memory \\
        tmpdir=tmp \\
        &> ${prefix}.bbnorm.log

    rm -rf tmp
    """
}