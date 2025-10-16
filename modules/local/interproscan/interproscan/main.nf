process INTERPROSCAN_INTERPROSCAN {

    tag "$meta.id"
    label 'process_low'
    label 'process_long'

    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    path(interproscan_database, stageAs: 'data')

    output:
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    tuple val(meta), path('*.xml') , optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    tuple val("${task.process}"), val('interproscan'), eval("interproscan.sh --version | sed '1!d; s/.*version //'"),    topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    """
     if [ -d 'data' ]; then
        # Find interproscan.properties to link data/ from work directory
        INTERPROSCAN_DIR="\$( dirname "\$( dirname "\$( which interproscan.sh )" )" )"
        echo \$INTERPROSCAN_DIR
        INTERPROSCAN_PROPERTIES="\$( find "\$INTERPROSCAN_DIR/share" -name "interproscan.properties" )"
        echo \$INTERPROSCAN_PROPERTIES
        cp "\$INTERPROSCAN_PROPERTIES" .
        sed -i "/^bin\\.directory=/ s|.*|bin.directory=\$INTERPROSCAN_DIR/bin|" interproscan.properties
        cat interproscan.properties
        export INTERPROSCAN_CONF=interproscan.properties
    fi # else use sample DB included with conda ( testing only! )

    if ${is_compressed} ; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if grep -q '\\*' ${fasta_name}; then
      echo "Found * — removing trailing stop codons..."
      sed -i '/^>/! s/\\*\$//' ${fasta_name}
    fi

    interproscan.sh \\
        --cpu ${task.cpus} \\
        --input ${fasta_name} \\
        ${args} \\
        --output-file-base ${prefix}

    rm -rf data
    """
}
