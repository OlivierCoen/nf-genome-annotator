nextflow.enable.types = true

process ORTHODB_MAKECLADEDB {

    label 'process_medium'
    tag "$orthodb_clade"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b4860b3d61daf2a3c47bf82ff4fe5d69a058ffd25ed4212fb50c34c5baa1ddd/data':
        'TODO' }"

    input:
        record(
            orthodb_clade: String,
            orthodb_excluded_clades: List<String>,
            orthodb_excluded_species: List<String>
        )

    output:
        record(
            orthodb_clade: orthodb_clade,
            orthodb_excluded_clades: orthodb_excluded_clades,
            orthodb_excluded_species: orthodb_excluded_species,
            orthodb_proteins: file("${orthodb_clade}.orthodb_proteins.faa.gz")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )           >> 'versions'
        tuple( "${task.process}", 'aria2',  eval("aria2c -v | head -1 | sed 's/aria2 version //g'") ) >> 'versions'
        tuple( "${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'") )          >> 'versions'
        tuple( "${task.process}", 'seqkit', eval("seqkit version | sed 's/^.*v//'") )                 >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def odb_levels_species_urls = [
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_level2species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_levels.tab.gz"
    ].join(' ').trim()
    def odb_protein_fasta_url = "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_aa_fasta.gz"
    def excluded_clades_arg = orthodb_excluded_clades ? "--exclude ${orthodb_excluded_clades.join(',')}" : ""
    def excluded_species_arg = orthodb_excluded_species ? "--excludeSpecies ${orthodb_excluded_species.join(',')}" : ""
    def nb_splits = Math.min(16, task.cpus.toInteger())
    def nb_max_connections = Math.min(16, task.cpus.toInteger())
    """
    for url in ${odb_levels_species_urls}
    do
        outfile=\$(basename \$url)

        echo "Downloading \$url to \$outfile"
        aria2c \\
            -x ${nb_splits} \\
            -s ${nb_max_connections} \\
            -o \$outfile \\
            \$url

        echo "Decompressing \$outfile"
        pigz -d \$outfile
    done

    PROTEIN_FASTA_FILE=odb12v2.faa.gz
    SPECIES_IDS_FILE=species_ids.txt
    
    echo "Downloading ${odb_protein_fasta_url}"
    aria2c \\
        -x ${nb_splits} \\
        -s ${nb_max_connections} \\
        -o \$PROTEIN_FASTA_FILE \\
        $args \\
        $odb_protein_fasta_url

    echo "Done downloading"

    echo "Getting clade specific species ids"
    get_orthodb_clade_specific_species_ids.py \\
        --clade "$orthodb_clade" \\
        --levels odb12v2_levels.tab \\
        --levels2species odb12v2_level2species.tab \\
        --species odb12v2_species.tab \\
        --outfile \$SPECIES_IDS_FILE \\
        $excluded_clades_arg \\
        $excluded_species_arg

    echo "Filtering \$PROTEIN_FASTA_FILE"
    seqkit grep -f \$SPECIES_IDS_FILE \$PROTEIN_FASTA_FILE | pigz > ${orthodb_clade}.orthodb_proteins.faa.gz
    
    echo "Removing intermediate files"
    rm \$PROTEIN_FASTA_FILE odb12v2_levels.tab odb12v2_level2species.tab odb12v2_species.tab
    """

}
