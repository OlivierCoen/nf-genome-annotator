nextflow.enable.types = true

process RED_RED {
    tag "$id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/78b7b35e8581533249f184933bc11db9f8c895d3550bdb6d2ae94853f088db65/data':
        'community.wave.seqera.io/library/red:2018.09.10--e81556edfad56018' }"

    input:
        record(id: String, fasta: Path)

    output:
        record(id: id, softmasked: file("${fasta.baseName}.softmasked.fa"))

    topic:
        tuple('red', id, files("repeats/*.rpt", optional: true)) >> 'additional_results'
        tuple("${task.process}", 'red', eval("Red 2>&1 | grep Version | cut -d' ' -f2")) >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${id}"
    """
    # adapted from https://github.com/Gaius-Augustus/BRAKER4/blob/main/rules/preprocessing/run_red_masking.smk
    # credits: KatharinaHoff

    # Red expects a directory of FASTA files with .fa extension
    mkdir -p input/ masked/ repeats/

    ln -s \$PWD/${fasta} input/genome.fa

    Red \\
        -gnm input \\
        -msk masked \\
        -rpt repeats

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
    """

}
