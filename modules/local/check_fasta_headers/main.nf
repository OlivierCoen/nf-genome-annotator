process CHECK_FASTA_HEADERS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/039340bbfffa2261c6ab74f8f66c68151b9116f7784d82cdfe167b2dc90eca1e/data':
        'community.wave.seqera.io/library/biopython_python:f180d02b12dd489c' }"

    input:
    tuple val(meta), path(fasta)
    val fix

    output:
    tuple val(meta), path("*headers_cleaned.{fasta,fa,fas,fna,faa}.gz"), emit: fasta, optional: true
    tuple val("${task.process}"), val('python'),eval("python3 --version | sed 's/Python //'"),            topic: versions
    tuple val("${task.process}"), val('Bio'),    eval('python3 -c "import Bio; print(Bio.__version__)"'), topic: versions

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

    if [ "${is_compressed}" == "true" ]; then
        echo "Removing ${fasta_name}"
        rm ${fasta_name}
        if [ "${fix}" == "true" ]; then
            echo "Compressing \$outfile"
            gzip \$outfile
        fi
    fi
    """

}
