nextflow.enable.types = true

process INTERPROSCAN5_INTERPROSCAN {

    tag "$id"
    label 'process_high'

    // there are issues with the interproscan db whith symlinks... 
    // see if it's really not possible to use simlinks
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
        record(
            id: String,     
            fasta: Path
        )
        interproscan5_db: Path

    stage:
        stageAs interproscan5_db, 'data'

    output:
        record(
            id: id,
            gff: file('*.gff3', optional: true)
        )
    
    topic:
        tuple('interproscan', id, file('*.tsv', optional: true))  >> 'additional_results'
        tuple('interproscan', id, file('*.xml', optional: true))  >> 'additional_results'
        tuple('interproscan', id, file('*.json', optional: true)) >> 'additional_results'

        tuple("${task.process}", 'interproscan', eval("interproscan.sh --version | sed '1!d; s/.*version //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
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
