nextflow.enable.types = true

process INTERPROSCAN5_INTERPROSCAN {

    tag "$id"
    label 'process_high'

    // no conda package for the latest releases, and anyway it's really tough to have interproscan work with conda...
    //conda "${moduleDir}/environment.yml"
    container "docker.io/interpro/interproscan:5.78-109.0"

    input:
        record(
            id: String,     
            fasta: Path
        )
        interproscan5_db: Path

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
    # Note: since /opt/interproscan is read-only, we need to provide the path to the data through the conf 
    # No parameter is exposed to do that...
    cp /opt/interproscan/interproscan.properties .
    sed -i "/^data\\.directory=/ s|.*|data.directory=\${PWD}/data|" interproscan.properties
    export INTERPROSCAN_CONF=\${PWD}/interproscan.properties
        
    # uncompress proteome
    if ${is_compressed} ; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # remove trailing stars in sequences
    if grep -q '\\*' ${fasta_name}; then
      echo "Found * — removing trailing stop codons..."
      sed -i '/^>/! s/\\*\$//' ${fasta_name}
    fi

    mkdir tmp
    
    /opt/interproscan/interproscan.sh \\
        --cpu ${task.cpus} \\
        --input ${fasta_name} \\
        ${args} \\
        --output-file-base ${prefix} \\
        --tempdir tmp
    """
}

