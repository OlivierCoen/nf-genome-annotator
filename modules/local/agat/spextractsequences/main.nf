process AGAT_SPEXTRACTSEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gxf), path(genome)
    val codon_usage_id
    path config

    output:
    tuple val(meta), path("*.prot.faa"), emit: proteins
    tuple val(meta), path("agat.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_sp_extract_sequences.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = meta.main_annotation ? "${meta.id}" : "${gxf.baseName}"
    def config_arg  = config ? "-c ${config}" : ''
    def is_compressed = genome.getExtension() == "gz" ? true : false
    def genome_fasta = is_compressed ? genome.getBaseName() : genome
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${genome} > ${genome_fasta}
    fi

    agat_sp_extract_sequences.pl \\
        ${args} \\
        --gff ${gxf} \\
        --fasta ${genome_fasta} \\
        --protein \\
        --clean_final_stop \\
        --codon $codon_usage_id \\
        ${config_arg} \\
        --output ${prefix}.prot.faa \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "${genome}" == "${genome}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta
    """
}
