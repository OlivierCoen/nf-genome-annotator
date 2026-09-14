nextflow.enable.types = true

process FASTQC {
    tag "${id} :: ${read_id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
            : 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'}"

    input:
        record(
            id: String, 
            read_id: String,
            reads: Set<Path>
        )

    stage:
        stageAs reads, '?/*'

    //topic:
        //tuple(id, files("*.zip"))                                                                  >> 'fastqc_multiqc'
        //tuple('fastqc', id, files("*.html"))                                                       >> 'additional_results'
        //tuple("${task.process}", 'fastqc', eval('fastqc --version | sed "/FastQC v/!d; s/.*v//"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$read_id"

    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpus allows to stick to requested amount of RAM in the label
    def memory_in_mb = task.memory
        ? (task.memory.toUnit('MB') / task.cpus).intValue()
        : null
    // FastQC memory value allowed range (100 - 10000)
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    def fastqc_memory_arg = fastqc_memory ? "--memory ${fastqc_memory}" : ''

    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${fastqc_memory_arg} \\
        ${reads.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "$read_id"
    """
    touch ${prefix}.html
    touch ${prefix}.zip
    """
}
