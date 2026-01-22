process AGAT_SPFIXLONGESTORF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d0194019d3fcefea0875ce4703df12dd2244cf9b8932b935197e9063cbc118ae/data':
        'community.wave.seqera.io/library/agat:1.6.1--d39db4f54af12afb' }"

    input:
    tuple val(meta), path(gff), path(genome_fasta)
    path config

    output:
    tuple val(meta), path("*_longest_orfs_fixed.gff"), emit: gff
    tuple val(meta), path("agat.log"), emit: log
    tuple val("${task.process}"), val('agat'), eval("agat_sp_fix_longest_ORF.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'"),    topic: versions

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    """
    agat_sp_fix_longest_ORF.pl \\
        --gff $gff \\
        --fasta $genome_fasta \\
        --model 1,2,3,4,6 \\
        ${config_param} \\
       ${args} \\
        --output ${prefix}_longest_orfs_fixed.gff \\
        > agat.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_longest_orfs_fixed.gff
    """
}
