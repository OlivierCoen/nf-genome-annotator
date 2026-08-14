nextflow.enable.types = true

process AGAT_SPFUNCTIONALSTATISTICS {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91396f03e6b1ac316141776c8782c8cbc085e53c6fc390f67aa272e5d4337813/data' :
        'community.wave.seqera.io/library/agat_pyyaml:b4d19f33ad15b73b' }"

    input:
        record(
            id: String, 
            gff: Path, 
            genome_size: Integer?
        )

    output:
        record(
            id: id,
            gff_stats: file("*.yaml")
        )

    topic:
        tuple(id, files("*_gff_stats.csv", optional: true)) >> 'multiqc'
        tuple("${task.process}", 'agat', eval("agat_sp_functional_statistics.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "$id"
    def genome_size_arg = genome_size ? "--gs ${genome_size}" : ''
    """
    agat_sp_functional_statistics.pl \\
        --gff ${gff} \\
        ${genome_size_arg} \\
        --output ${prefix}.gtf_func_stats.txt \\
        --yaml \\
        ${args}

    # parse yaml file
    parse_gff_stat_file.py \\
        --gff ${prefix}.gtf_stats.txt.yaml \\
        --prefix ${prefix}
    """
}
