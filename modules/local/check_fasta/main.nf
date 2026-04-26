process CHECK_FASTA {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/039340bbfffa2261c6ab74f8f66c68151b9116f7784d82cdfe167b2dc90eca1e/data':
        'community.wave.seqera.io/library/biopython_python:f180d02b12dd489c' }"

    input:
    tuple val(meta), path(fasta)
    val type
    val fix_headers
    val fix_sequences

    output:
    tuple val(meta), path("*.cleaned.{fasta,fa,fas,fna,faa}*"), emit: fasta, optional: true
    tuple val("${task.process}"), val('python'),eval("python3 --version | sed 's/Python //'"),            topic: versions
    tuple val("${task.process}"), val('Bio'),    eval('python3 -c "import Bio; print(Bio.__version__)"'), topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.cleaned"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta.name
    def fasta_ext = fasta_name.tokenize('.')[-1]
    def fix_headers_arg = fix_headers ? "--fix-headers": ""
    def fix_sequences_arg = fix_sequences ? "--fix-sequences": ""
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    outfile=${prefix}.${fasta_ext}

    check_fasta.py \\
        --in $fasta_name \\
        --out \$outfile \\
        --type $type \\
        $fix_headers_arg \\
        $fix_sequences_arg

    if [ "${is_compressed}" == "true" ]; then
        echo "Removing ${fasta_name}"
        rm ${fasta_name}
        if [ -f \$outfile ]; then
            echo "Compressing \$outfile"
            gzip \$outfile
        fi
    fi
    """

}
