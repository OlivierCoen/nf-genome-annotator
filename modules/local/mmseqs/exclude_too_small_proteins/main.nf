process MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(db)
    val(min_prot_db_seq_length)

    output:
    tuple val(meta), path("*.${suffix}"), emit: filtered_db
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    script:
    def args  = task.ext.args  ?: '' // args is used for the main arguments of the tool
    """
    awk '\$3 > 102 {print \$1}' seqDb.index > ids.gt100
    """

}
