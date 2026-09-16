nextflow.enable.types = true 

process GFFREAD {

    tag "${id} :: ${gff.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f2e0989f34dd166393d87c5265fe5860224ab0cf92061e45a8b21602605b134/data':
        'community.wave.seqera.io/library/gffread_gzip:3a40774e2f2a9631' }"

    input:
        record(
            id: String,
            gff: Path,
            fasta: Path
        )

    output:
        record(
            id: id,
            proteins: file("*.proteins.faa")
        )

    //topic:
    // tuple('gffread', id, file(".mrna.fna")) >> 'additional_results'
    // tuple('gffread', id, file(".cds.fna")) >> 'additional_results'
    //    tuple("${task.process}", 'agat', eval("agat_sp_extract_sequences.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = "${gff.baseName}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def genome_fasta = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${genome_fasta}
    fi

    gffread \\
    $gff \\
      -g $genome_fasta \\
      -w ${prefix}.mrna.fna \\
      -x ${prefix}.cds.fna \\
      -y ${prefix}.proteins.faa \\
      -W \\
      $args
    """
}
