nextflow.enable.types = true

process EARLGREY_EARLGREY {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f798eb5fb7d14db38bde92881b8c76e6843f22b88e61d7b2224be2b7a3fb60a4/data':
        'community.wave.seqera.io/library/earlgrey:7.3.1--739efdc54cd5d7d9' }"

    input:
        record(
            id: String,
            fasta: Path,
            dfam_db: Path
        )

    output:
        record(
            id: id,
            softmasked: file("${prefix}.masked")
        )
  
    topic:
        tuple("${task.process}", 'earlgrey', eval("earlGrey | grep version | sed 's/earlGrey version //g'")) >> 'versions'

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "$id"
    """
    # configuring RepeatMasker
    BIN_DIR=\$(dirname \$(which RepeatMasker))
    REPEATMASKER_SHARE_DIR=\$(dirname \$BIN_DIR)/share/RepeatMasker
	perl \${REPEATMASKER_SHARE_DIR}/configure -libdir $dfam_db

    earlGrey \\
        -g $fasta \\
        -o results \\
        -s ${prefix} \\
        -t ${task.cpus} \\
        -d \\
        ${args}
    """

}
