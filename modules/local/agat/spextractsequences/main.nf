process AGAT_SPEXTRACTSEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gxf), path(fasta)
    val codon_usage_id
    path config

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val("${task.process}"), val('agat'), eval("agat_sp_extract_sequences.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def config_arg  = config ? "-c ${config}" : ''
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_extract_sequences.pl \\
        ${args} \\
        --gff ${gxf} \\
        --fasta ${fasta} \\
        --protein \\
        --codon $codon_usage_id \\
        ${config_arg} \\
        --output ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "${fasta}" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta
    """
}
