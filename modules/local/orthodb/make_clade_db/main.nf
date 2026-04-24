process ORTHODB_MAKECLADEDB {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/7563efb3dae4e088733d162f0fa33a00c5467e3e448eb263c9148032d43abecc/data':
        'community.wave.seqera.io/library/aria2_pigz_python:13735d8e32a1c063' }"

    input:
    val clade
    val excluded_clades
    val excluded_species

    output:
    path("${clade}.orthodb_proteins.faa.gz"), emit: db
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),           topic: versions
    tuple val("${task.process}"), val('aria2'),  eval("aria2c -v | head -1 | sed 's/aria2 version //g'"), topic: versions
    tuple val("${task.process}"), val('pigz'),   eval("pigz --version 2>&1 | sed 's/pigz //g'"),          topic: versions

    script:
    def orthodb_file_urls = [
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_aa_fasta.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_level2species.tab.gz",
        "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_levels.tab.gz"
    ].join(' ').trim()
    def excluded_clades_arg = excluded_clades != "" ? "--exclude $excluded_clades": ""
    def excluded_species_arg = excluded_species != "" ? "--excludeSpecies $excluded_species" : ""
    """
    for url in ${orthodb_file_urls}
    do
        outfile=\$(basename \$url)

        echo "Downloading \$url to \$outfile"
        aria2c \\
            --split ${task.cpus} \\
            --max-connection-per-server ${task.cpus} \\
            --optimize-concurrent-downloads \\
            --check-integrity=true \\
            --max-tries=10 \\
            --retry-wait=30 \\
            --timeout=60 \\
            \$url \\
            -o \$outfile

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
