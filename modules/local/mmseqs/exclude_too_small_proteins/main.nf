nextflow.enable.types = true

process MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7c4fe261ef11dff32838cf078e2965a48866ca957c798d488c3d4a61437f9e50/data' :
        'community.wave.seqera.io/library/mmseqs2_gawk:544919ce21ad9400' }"

    input:
        record(
            id: String,
            db: Path
        )
        min_len: Integer

    output:
        record(
            id: id,
            db: file("${db}.gt${min_len}")
        )

    topic:
        tuple("${task.process}", 'mmseqs', eval('mmseqs version'))                                 >> 'versions'
        tuple("${task.process}", 'gawk',   eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'

    script:
    // see https://github.com/soedinglab/MMseqs2/wiki#manipulating-databases
    def actual_min_len = min_len + 2
    """
    awk '\$3 > ${actual_min_len} {print \$1}' ${db}.index > ids.gt${min_len}
    mmseqs createsubdb ids.gt${min_len} ${db} ${db}.gt${min_len}
    """

}
