process HELIXER {
    tag "$meta.id"
    label 'process_gpu'

    container 'docker://gglyptodon/helixer-docker:helixer_v0.3.3_cuda_11.2.0-cudnn8'

    input:
    tuple val(meta), path(fasta)
    val  lineage

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    path "versions.yml",             emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Helixer.py \\
        --fasta-path $fasta \\
        --lineage $lineage \\
        --gff-output-path ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        helixer: \$(Helixer.py --version 2>&1 | head -1)
    END_VERSIONS
    """
}
