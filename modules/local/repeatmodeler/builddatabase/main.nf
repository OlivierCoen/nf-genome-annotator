process REPEATMODELER_BUILDDATABASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b9637a9d4d72993f80c5014ece53d39161871c94845f420a5313a31a7ae5d2a/data':
        'community.wave.seqera.io/library/repeatmodeler:2.0.7--136d59ab97ab30de' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.*")    , emit: db
    tuple val("${task.process}"), val('repeatmodeler'), eval("RepeatModeler --version | sed 's/RepeatModeler version //'"), topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    BuildDatabase \\
        -name $prefix \\
        $fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
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
