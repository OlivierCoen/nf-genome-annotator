nextflow.enable.types = true

process ORTHODB_MAKECLADEDB {

    label 'process_download_db'
    tag "$clade"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/7563efb3dae4e088733d162f0fa33a00c5467e3e448eb263c9148032d43abecc/data':
        'community.wave.seqera.io/library/aria2_pigz_python:13735d8e32a1c063' }"

    input:
        record(
            clade: String,
            excluded_clades: Iterable<String>,
            excluded_species: Iterable<String>
        )

    output:
        record(
            clade: clade,
            excluded_clades: excluded_clades,
            excluded_species: excluded_species,
            orthodb_proteins: file("${clade}.orthodb_proteins.faa.gz")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )           >> 'versions'
        tuple( "${task.process}", 'aria2',  eval("aria2c -v | head -1 | sed 's/aria2 version //g'") ) >> 'versions'
        tuple( "${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'") )          >> 'versions'

    script:
    def orthodb_file_urls = [
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_aa_fasta.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_level2species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_levels.tab.gz"
    ].join(' ').trim()
    def excluded_clades_arg = excluded_clades ? "--exclude ${excluded_clades.join(',')}" : ""
    def excluded_species_arg = excluded_species ? "--excludeSpecies ${excluded_species.join(',')}" : ""
    def nb_splits = Math.min(16, task.cpus.toInteger())
    def nb_max_connections = Math.min(16, task.cpus.toInteger())
    """
    for url in ${orthodb_file_urls}
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

    echo "Done downloading and decomprissing"
    # Renaming fasta file
    mv odb12v2_aa_fasta odb12v2_all.faa

    echo "Filtering odb12v2_all.faa"
    select_clade_from_orthodb.py \\
        odb12v2_all.faa \\
        odb12v2_levels.tab \\
        odb12v2_level2species.tab \\
        odb12v2_species.tab \\
        --clade "$clade" \\
        $excluded_clades_arg \\
        $excluded_species_arg \\
        > ${clade}.orthodb_proteins.faa

    echo "Compressing ${clade}.orthodb_proteins.faa"
    pigz ${clade}.orthodb_proteins.faa

    echo "Removing intermediate files"
    rm odb12v2_all.faa odb12v2_levels.tab odb12v2_level2species.tab odb12v2_species.tab
    """

}
