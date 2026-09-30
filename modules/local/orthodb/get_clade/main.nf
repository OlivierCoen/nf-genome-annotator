nextflow.enable.types = true

process GET_ORTHODB_CLADE {

    tag "$species"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2ec166d51ecd375d538623ea698a7d1003aea274329d2c5623c203755d25c839/data':
        'community.wave.seqera.io/library/python_httpx_tenacity_pigz:aa5cc8a91450ff8d' }"

    input:
        species: String
        orthodb_clades_list: List<Path> // because we used collect()

    output:
        record(
            species: species,
            orthodb_clade: env("ORTHODB_CLADE") == 'None' ? null : env("ORTHODB_CLADE")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'
        tuple( "${task.process}", 'httpx',  eval('python3 -c "import httpx; print(httpx.__version__)"') ) >> 'versions'
        tuple( "${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'") )              >> 'versions'

    script:
    def orthodb_clades = orthodb_clades_list[0]
    def orthodb_clades_is_compressed = orthodb_clades.getExtension() == "gz" ? true : false
    def orthodb_clades_file = orthodb_clades_is_compressed ? orthodb_clades.baseName : orthodb_clades
    """
    # if the orthodb clades file is compressed, decompress it
    if [ "${orthodb_clades_is_compressed}" == "true" ]; then
        pigz -c -d ${orthodb_clades} > ${orthodb_clades_file}
    fi

    get_orthodb_clade.py \\
        --species "${species}" \\
        --orthodb-clades ${orthodb_clades_file}

    ORTHODB_CLADE=\$(cat orthodb_clade.txt)
    """

}
