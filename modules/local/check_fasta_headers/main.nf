process CHECK_FASTA_HEADERS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386577baf029b1bd12bc58855e52a8b73575c62ccdb4c09dfcfa158dcf1bbcee/data':
        'community.wave.seqera.io/library/biopython:1.87--7dac4c636fcf94ae' }"

    input:
    tuple val(meta), path(fasta)
    val fix

    output:
    tuple val(meta), path("*headers_cleaned.{fasta,fa,fas,fna,faa}.gz"), emit: fasta, optional: true
    tuple val("${task.process}"), val('python'),eval("python3 --version | sed 's/Python //'"),            topic: versions
    //tuple val("${task.process}"), val('Bio'),    eval('python3 -c "import Bio; print(Bio.__version__)"'), topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.headers_cleaned"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta.name
    def fasta_ext = fasta_name.tokenize('.')[-1]
    def fix_arg = fix ? "--fix": ""
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    outfile=${prefix}.${fasta_ext}

    check_fasta_headers.py \\
        --in $fasta_name \\
        --out \$outfile \\
        $fix_arg

    if [ "${fix}" == "true" && ${is_compressed}" == "true" ]; then
        gzip \$outfile
    fi
    """

}
