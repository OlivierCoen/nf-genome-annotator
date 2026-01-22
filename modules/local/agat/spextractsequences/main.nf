process AGAT_SPEXTRACTSEQUENCES {
    tag "$meta.id"
    label 'process_single'

    // for now, the version of AGAT is 1.4.2 for this module
    // version 1.6.1 gives issues
    // TODO: update when issues are resolved
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7fd8135f4654d5e5827ef55f7a9eb17995b2e7d2a30633d672c127dde340072/data':
        'community.wave.seqera.io/library/agat:1.4.2--f0c60073d54a9afe' }"

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
    def prefix      = meta.final_annotation ? "${meta.id}" : "${gxf.baseName}"
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
