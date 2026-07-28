process GET_TAXONOMY_INFO {

    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2ec166d51ecd375d538623ea698a7d1003aea274329d2c5623c203755d25c839/data':
        'community.wave.seqera.io/library/python_httpx_tenacity_pigz:aa5cc8a91450ff8d' }"

    input:
    tuple val(meta), path(species)
    path busco_datasets
    path orthodb_clades

    output:
    tuple val(meta), path("*.cleaned.{fasta,fa,fas,fna,faa}*"), emit: fasta
    tuple val(meta), env("TAXID"),                              emit: taxid
    tuple val(meta), env("BUSCO_DATASET"),                      emit: busco_dataset
    tuple val(meta), env("ORTHODB_CLADE"),                     emit: orthodb_clade
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),               topic: versions
    tuple val("${task.process}"), val('httpx'),  eval('python3 -c "import httpx; print(httpx.__version__)"'), topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"),                topic: versions

    script:
    def orthodb_clades_is_compressed = orthodb_clades.getExtension() == "gz" ? true : false
    def orthodb_clades_file = orthodb_clades_is_compressed ? orthodb_clades.baseName : orthodb_clades
    """
    # removing header from busco dataset files
    grep -vP '^(\\d|#|D|\$)' ${busco_datasets} > busco_datasets.cleaned.txt

    if [ "${orthodb_clades_is_compressed}" == "true" ]; then
        pigz -c -d ${orthodb_clades_is_compressed} > ${orthodb_clades_file}
    fi

    get_taxonomy_info.py \\
        --species ${species} \\
        --busco-datasets busco_datasets.cleaned.txt \\
        --orthodb-clades ${orthodb_clades_file}

    TAXID=\$(cat found_taxid.txt)
    BUSCO_DATASET=\$(cat found_busco_dataset.txt)
    ORTHODB_CLADE=\$(cat found_orthodb_clade.txt)
    """

}
