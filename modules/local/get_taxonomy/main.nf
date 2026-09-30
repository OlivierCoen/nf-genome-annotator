nextflow.enable.types = true

process GET_TAXONOMY {

    tag "$species"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2ec166d51ecd375d538623ea698a7d1003aea274329d2c5623c203755d25c839/data':
        'community.wave.seqera.io/library/python_httpx_tenacity_pigz:aa5cc8a91450ff8d' }"

    input:
        species: String

    output:
        record(
            species: species,
            taxid: env("TAXID"),
            taxid_lineage: env("TAXID_LINEAGE")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'
        tuple( "${task.process}", 'httpx',  eval('python3 -c "import httpx; print(httpx.__version__)"') ) >> 'versions'
        tuple( "${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'") )              >> 'versions'

    script:
    """
    get_taxonomy_info.py \\
        --species "${species}"

    TAXID=\$(cat taxid.txt)
    TAXID_LINEAGE=\$(cat taxid_lineage.txt)
    """

}
