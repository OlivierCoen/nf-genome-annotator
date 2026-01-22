process AGAT_SPSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e7e31998b8d46049a71d795c8b63c59cc38823b0f06eca2869b98b1a1515cd9/data' :
        'community.wave.seqera.io/library/agat_python_pyyaml:bd0b5d997a84b24b' }"

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

    when:
    task.ext.when == null || task.ext.when

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
        ${args} \\
        > agat.log 2>&1

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
