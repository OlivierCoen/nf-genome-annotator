nextflow.enable.types = true 

process AGAT_SPEXTRACTSEQUENCES {

    tag "${id} :: ${gff.baseName}"
    label 'process_single'

    // for now, the version of AGAT is 1.4.2 for this module
    // version 1.6.1 gives issues
    // TODO: update when issues are resolved
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7fd8135f4654d5e5827ef55f7a9eb17995b2e7d2a30633d672c127dde340072/data':
        'community.wave.seqera.io/library/agat:1.4.2--f0c60073d54a9afe' }"

    input:
        record(
            id: String,
            gff: Path,
            fasta: Path
        )
        codon_usage_id: Integer?

    output:
        record(
            id: id,
            extracted_fasta: file("*.{faa,fna}")
        )

    //topic:
    //    tuple("${task.process}", 'agat', eval("agat_sp_extract_sequences.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = "${gff.baseName}"
    
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def genome_fasta = is_compressed ? fasta.getBaseName() : fasta

    def extract_proteins = args.contains("--protein") ? true : false
    def codon_usage_arg = extract_proteins ? "--codon $codon_usage_id" : ""
    def suffix          = extract_proteins ? "prot.faa" : "cds.fna"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${genome_fasta}
    fi

    agat_sp_extract_sequences.pl \\
        ${args} \\
        --gff ${gff} \\
        --fasta ${genome_fasta} \\
        ${codon_usage_arg} \\
        --output ${prefix}.${suffix}
    """
}
