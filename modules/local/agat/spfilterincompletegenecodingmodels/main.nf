process AGAT_SPFILTERINCOMPLETEGENECODINGMODELS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0':
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff), path(genome)
    path config

    output:
    tuple val(meta), path("*_incomplete_gene_models_fixed.gff"), emit: gff
    tuple val(meta), path("agat.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_sp_filter_incomplete_gene_coding_models.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    def is_compressed = genome.getExtension() == "gz" ? true : false
    def genome_fasta = is_compressed ? genome.getBaseName() : genome
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${genome} > ${genome_fasta}
    fi

    agat_sp_filter_incomplete_gene_coding_models.pl \\
        --gff $gff \\
        --fasta $genome_fasta \\
        --add_flag \\
        ${config_param} \\
       ${args} \\
        --output ${prefix}_incomplete_gene_models_fixed.gff \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_incomplete_gene_models_fixed.gff
    """
}
