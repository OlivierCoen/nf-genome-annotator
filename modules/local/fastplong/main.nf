nextflow.enable.types = true

process FASTPLONG {
    tag "${id} :: ${read_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['apptainer', 'singularity'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b4f7df1210b7d0545a6c62f3f914248aea0e5c5d2e18652af971adf21d0be17/data':
        'community.wave.seqera.io/library/fastplong_gzip:6a68cfd1567defb8' }"

    input:
        record(
            id: String,
            read_id: String,
            fastq: Path
        )

    output:
        record(
            id: id,
            read_id: read_id,
            fastq: file("*.fq.gz", optional: true)
        )

    topic:
        tuple('fastplong', id, file('fastplong.log',  optional: true)) >> 'logs'
        tuple('fastplong', id, file('fastplong.html', optional: true)) >> 'additional_results'
        tuple('fastplong', id, file('fastplong.json', optional: true)) >> 'additional_results'
        tuple("${task.process}", 'fastplong', eval("fastplong --version | cut -d' ' -f2")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${read_id}.fastplong"
    """
    fastplong \\
        $args \\
        --in $fastq \\
        --out ${prefix}.fq.gz \\
        --thread ${task.cpus} \\
        | tee fastplong.log
        
    # checking if filtered file is empty
    if [ \$(zcat ${prefix}.fq.gz | head -c 1 | wc -c) -eq "0" ]; then 
        echo "${prefix}.fq.gz is empty... Removing file to avoid downstream issues"
        rm ${prefix}.fq.gz
    fi
    """
}
