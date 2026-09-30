nextflow.enable.types = true

process BUSCO_GET_LINEAGE {

    tag "$species"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2ec166d51ecd375d538623ea698a7d1003aea274329d2c5623c203755d25c839/data':
        'community.wave.seqera.io/library/python_httpx_tenacity_pigz:aa5cc8a91450ff8d' }"

    input:
        species: String
        busco_datasets_list: List<Path> // because we used collect()

    output:
        record(
            species: species,
            busco_lineage: env("BUSCO_LINEAGE") == 'None' ? null : env("BUSCO_LINEAGE")
        )

    topic:
        tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )               >> 'versions'
        tuple( "${task.process}", 'httpx',  eval('python3 -c "import httpx; print(httpx.__version__)"') ) >> 'versions'
        tuple( "${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'") )              >> 'versions'

    script:
    def busco_datasets = busco_datasets_list[0]
    """
    # removing header from busco dataset files
    grep -vP '^(\\d|#|D|\$)' ${busco_datasets} > busco_datasets.cleaned.txt

    get_busco_lineage.py \\
        --species "${species}" \\
        --busco-datasets busco_datasets.cleaned.txt

    BUSCO_LINEAGE=\$(cat busco_lineage.txt)
    """

}
