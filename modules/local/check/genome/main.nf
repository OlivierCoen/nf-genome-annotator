nextflow.enable.types = true

process CHECK_GENOME {

    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/039340bbfffa2261c6ab74f8f66c68151b9116f7784d82cdfe167b2dc90eca1e/data':
        'community.wave.seqera.io/library/biopython_python:f180d02b12dd489c' }"

    input:
        record( id: String, genome_fasta: Path )

    output:
        checked = record( id: id, genome_fasta: file("output/*.cleaned.*") )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )           >> 'versions'
        tuple( "${task.process}", 'Bio',    eval('python3 -c "import Bio; print(Bio.__version__)"') ) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${id}.cleaned"
    def is_compressed = genome_fasta.getExtension() == "gz" ? true : false
    def genome_name = is_compressed ? genome_fasta.getBaseName() : genome_fasta.name
    def genome_ext = genome_name.tokenize('.')[-1]
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${genome_fasta} > ${genome_name}
    fi

    mkdir output
    outfile=output/${prefix}.${genome_ext}

    check_genome.py \\
        --in $genome_name \\
        --out \$outfile

    if [ "${is_compressed}" == "true" ]; then
        echo "Removing ${genome_name}"
        rm ${genome_name}
        if [ -f \$outfile ]; then
            echo "Compressing \$outfile"
            gzip \$outfile
        fi
    fi
    """

}
