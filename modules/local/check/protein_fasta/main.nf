nextflow.enable.types = true

process CHECK_PROTEIN_FASTA {

    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/039340bbfffa2261c6ab74f8f66c68151b9116f7784d82cdfe167b2dc90eca1e/data':
        'community.wave.seqera.io/library/biopython_python:f180d02b12dd489c' }"

    input:
        record(
            id: String,
            fasta: Path
        )
        minlen: Integer

    output:
        record(
            id: id,
            fasta: file("*.cleaned.{fasta,fa,fas,fna,faa}*")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )           >> 'versions'
        tuple( "${task.process}", 'Bio',    eval('python3 -c "import Bio; print(Bio.__version__)"') ) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${id}.cleaned"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta.name
    def fasta_ext = fasta_name.tokenize('.')[-1]
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    outfile=${prefix}.${fasta_ext}

    check_protein_fasta.py \\
        --in $fasta_name \\
        --out \$outfile \\
        --minlen $minlen

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
