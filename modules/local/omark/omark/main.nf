nextflow.enable.types = true

process OMARK_OMARK {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7082d2bbe6948025958cd10031f313633a326e5b3c92b122af1224620de0e9c/data'
        : 'community.wave.seqera.io/library/omark:0.5.0--a004ed2a54d1ff9a'}"

    input:
        record(
            id: String,
            fasta: Path,
            omamer: Path,
            transcript_isoforms: Path
        )
        omamer_db: Path
    
    topic:
        tuple('omark', id, files("*_omark_out/*")) >> 'additional_results'
        tuple("${task.process}", 'omark', eval('omark --version')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    # for matplotlib
    export MPLCONFIGDIR=\${PWD}
    
    omark \\
        --file ${omamer} \\
        --database ${omamer_db} \\
        --ete_ncbi_db \${PWD}/.etetoolkit \\
        --isoform_file ${transcript_isoforms} \\
        --og_fasta ${fasta} \\
        --outputFolder ${prefix}_omark_out \\
        ${args}
    """

   
}
