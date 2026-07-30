nextflow.enable.types = true

process REPEATMODELER_BUILDDATABASE {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4618e2da2755ead8794cee7e6fe99de4ef739a97364e08ba9fa8b2988cdba123/data':
        'community.wave.seqera.io/library/repeatmodeler:2.0.9--7529329ebd736619' }"

    input:
        record(id: String, fasta: Path)

    output:
        record(id: id, db: files("${prefix}.*"))
    
    topic:
        tuple( "${task.process}", 'repeatmodeler', eval("RepeatModeler --version | sed 's/RepeatModeler version //'") ) >> 'versions'
        
    script:
    prefix = task.ext.prefix ?: "${id}"
    """
    BuildDatabase \\
        -name $prefix \\
        $fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.nhr
    touch ${prefix}.nin
    touch ${prefix}.njs
    touch ${prefix}.nnd
    touch ${prefix}.nni
    touch ${prefix}.nog
    touch ${prefix}.nsq
    touch ${prefix}.translation
    """
}
