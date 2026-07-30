nextflow.enable.types = true

process REPEATMASKER_REPEATMASKER {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f55b93073469a6e7aea3b0db14cf8c7181528ccecfca7fc06a10ae16818cf47/data':
        'community.wave.seqera.io/library/repeatmasker:4.2.4--eb5ce0b75f3001ca' }"

    input:
        record(
            id: String, 
            fasta: Path, 
            lib: Path
        )

    output:
        record(
            id: id, 
            softmasked: file("*.masked"), 
            repeats_gff: file("*.gff", optional: true)
        )
   
    topic:
        tuple('repeatmasker', file("*.tbl"))                 >> 'additional_results'
        tuple('repeatmasker', file("*.out"))                 >> 'logs'
        tuple("${task.process}", 'repeatmasker', eval("RepeatMasker -v | sed 's/RepeatMasker version //1'")) >>  'versions'

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${id}"
    def lib_arg = lib               ? "-lib $lib"   : ''
    
    def is_compressed  = fasta.getExtension() == "gz" ? true : false
    def fasta_name     = is_compressed                ? fasta.getBaseName() : fasta.name
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi
    
    RepeatMasker \\
        $lib_arg \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        ${args} \\
        ${fasta_name}

    mv $prefix/${fasta_name}.masked  ${prefix}.masked
    mv $prefix/${fasta_name}.out     ${prefix}.out
    mv $prefix/${fasta_name}.tbl     ${prefix}.tbl
    mv $prefix/${fasta_name}.out.gff ${prefix}.gff       || echo "GFF is not produced"
    """

    stub:
    def prefix      = task.ext.prefix       ?: "${id}"
    def args        = task.ext.args         ?: ''
    def touch_gff   = args.contains('-gff') ? "touch ${prefix}.gff" : ''

    """
    touch ${prefix}.masked
    touch ${prefix}.out
    touch ${prefix}.tbl
    $touch_gff

    """
}
