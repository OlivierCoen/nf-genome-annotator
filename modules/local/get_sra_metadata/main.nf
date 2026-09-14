nextflow.enable.types = true

process GET_SRA_METADATA {

    maxForks 1 // to avoid issues with the NCBI E-Utilities API

    label 'process_high'

    tag "$taxid"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/7140956576e6779ff2bea610b8e41bde95d4909d67b31634d0ddb6dba50aef5a/data':
        'community.wave.seqera.io/library/requests_tenacity_xmltodict:9e74a2aeeb88aab9' }"

    input:
        taxid: String
        nb_short_read_sra_datasets: Integer
        nb_long_read_sra_datasets: Integer
        sra_max_size: String
        sra_allow_single_end: Boolean

    output:
        record(
            taxid: taxid,
            short_read_sra_ids_file: file("sra_ids.short_read.txt", optional: true),
            long_read_sra_ids_file: file("sra_ids.long_read.txt", optional: true)
        )

    topic:
        tuple('sra_metadata', taxid, files("sra_metadata.*.json", optional: true)) >> 'additional_results'
        tuple("${task.process}", 'python',   eval("python3 --version | sed 's/Python //'"))                                             >> 'versions'                  
        tuple("${task.process}", 'requests', eval('python3 -c "import requests; print(requests.__version__)"'))                         >> 'versions'
        tuple("${task.process}", 'tenacity', eval('python3 -c "from importlib.metadata import version; print(version(\'tenacity\'))"')) >> 'versions'
        tuple("${task.process}", 'xmltodict', eval('python3 -c "import xmltodict; print(xmltodict.__version__)"'))                      >> 'versions'    
    
    script:
    def paired_only_arg = sra_allow_single_end ? "" : "--paired-only"
    """
    get_sra_metadata.py \\
        --taxid $taxid \\
        --max-short-reads $nb_short_read_sra_datasets \\
        --max-long-reads $nb_long_read_sra_datasets \\
        --max-size $sra_max_size \\
        $paired_only_arg
    """

}}

