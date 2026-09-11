nextflow.enable.types = true

process RED {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2f4bbe923930c05cb260401b16ef62605e509270a7242dafa618291a57c1fb2/data':
        'community.wave.seqera.io/library/red_gzip:dd4c51f9a4d75d0d' }"

    input:
        record(id: String, fasta: Path)

    output:
        record(id: id, softmasked: file("${fasta.baseName}.softmasked.fa.gz"))

    topic:
        tuple('red', id, files("repeats/*.rpt", optional: true)) >> 'additional_results'
        tuple("${task.process}", 'red', eval("Red 2>&1 | grep Version | cut -d' ' -f2")) >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${id}"
    def is_compressed      = fasta.getExtension() == "gz"    ? true : false
    """
    # adapted from https://github.com/Gaius-Augustus/BRAKER4/blob/main/rules/preprocessing/run_red_masking.smk
    # credits: KatharinaHoff

    # Red expects a directory of FASTA files with .fa extension
    mkdir -p input/ masked/ repeats/

    ################################################
    # DECOMPRESSING INPUT FASTA IF NEEDED
    ################################################
    
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > input/genome.fa
    else
        ln -s ${fasta} input/genome.fa
    fi

    ################################################
    # RUNNING RED
    ################################################

    Red \\
        -gnm input \\
        -msk masked \\
        -rpt repeats

    ################################################
    # HANDLE RED OUTPUT
    ################################################

    # Red produces .msk files with the same basename as input
    if [ ! -f masked/genome.msk ]; then
        echo "Red failed to produce masked output"
        exit 1
    fi

    softmasked_genome=${fasta.baseName}.softmasked.fa
    mv masked/genome.msk \$softmasked_genome

    # Count masked bases for the log
    TOTAL=\$(grep -v '^>' \$softmasked_genome | tr -d '\n' | wc -c)
    MASKED=\$(grep -v '^>' \$softmasked_genome | tr -d '\n' | tr -cd 'a-z' | wc -c)
    PCT=\$(awk "BEGIN {printf \\"%.1f\\", 100.0*\$MASKED/\$TOTAL}")
    echo "Masked \$MASKED / \$TOTAL bp (\$PCT%)"
    echo "Red masking complete for ${id}"

    echo "Compressing softmasked genome"
    gzip \$softmasked_genome
    """

}
