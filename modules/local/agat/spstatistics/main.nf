process AGAT_SPSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91396f03e6b1ac316141776c8782c8cbc085e53c6fc390f67aa272e5d4337813/data' :
        'community.wave.seqera.io/library/agat_pyyaml:b4d19f33ad15b73b' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.txt"),                       emit: stats_txt
    tuple val(meta), path("*.yaml"),                      emit: stats_yaml
    path("*_with_isoforms_mrna_gff_stats.csv"),          topic: mqc_mrna_with_isoforms_gff_stats,          optional: true
    path("*_with_isoforms_rna_gff_stats.csv"),           topic: mqc_rna_with_isoforms_gff_stats,           optional: true
    path("*_with_isoforms_transcript_gff_stats.csv"),    topic: mqc_transcript_with_isoforms_gff_stats,    optional: true
    path("*_without_isoforms_mrna_gff_stats.csv"),       topic: mqc_mrna_without_isoforms_gff_stats,       optional: true
    path("*_without_isoforms_rna_gff_stats.csv"),        topic: mqc_rna_without_isoforms_gff_stats,        optional: true
    path("*_without_isoforms_transcript_gff_stats.csv"), topic: mqc_transcript_without_isoforms_gff_stats, optional: true
    tuple val("${task.process}"), val('agat'), eval("agat_sp_statistics.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"), topic: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = meta.final_annotation ? "${meta.id}.final_annotation" : "${gff.simpleName}.intermediate_annotation"
    def genome_size_arg = meta.genome_size ? "--gs ${meta.genome_size}" : ''
    """
    agat_sp_statistics.pl \\
        --gff ${gff} \\
        ${genome_size_arg} \\
        --output ${prefix}.gtf_stats.txt \\
        --yaml \\
        ${args}

    # parse yaml file
    parse_gff_stat_file.py \\
        --gff ${prefix}.gtf_stats.txt.yaml \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats.txt
    touch ${prefix}.stats.yaml
    """
}
