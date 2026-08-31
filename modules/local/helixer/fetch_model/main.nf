nextflow.enable.types = true

process HELIXER_FETCHMODEL {
    tag "$meta.id"
    label 'process_single'

    // Helixer does not provide a conda package
    container "docker.io/gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1"

    input:
        lineage: String

    output:
        record(
            lineage: lineage,
            model: file("models/${lineage}/*")
        )

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (workflow.containerEngine == 'singularity') {
        log.warn("Running Helixer with Singularity is not recommended since you may encounter issues with permissions. " +
                 "Consider using Apptainer instead. See https://github.com/gglyptodon/helixer-docker for more information.")
    }
    """
    Helixer/scripts/fetch_helixer_models.py \\
        --lineage $lineage \\
        --custom-path models
    """

}
