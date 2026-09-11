nextflow.enable.types = true

process GET_RANDOM_SAMPLE {

    tag "${taxid} :: ${type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7977dd87c418c76b502f06c515e51fd8eedd1aeee9c76378d0ec0bf1a8e4c3e/data':
        'community.wave.seqera.io/library/python:3.14.7--629fbce4f44dac86' }"

    input:
        record( 
            taxid: String, 
            type: String,
            input_file: Path,
            nb_to_sample: Integer
        )
        seed: Integer

    output:
        record( taxid: taxid, type: type, sampled: file("*.sampled.txt") )

    //topic:
    //    tuple( "${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )           >> 'versions'

    script:
    def prefix = "${input_file.baseName}.sampled"
    """
    get_random_sample.py \\
        --in $input_file \\
        --seed $seed \\
        --nb $nb_to_sample \\
        --out ${prefix}.txt
    """

}
