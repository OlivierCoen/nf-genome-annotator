nextflow.enable.types = true

process HELIXER_FETCHMODEL {
    tag "$lineage"
    label 'process_single'

    // Helixer does not provide a conda package
    container "docker.io/gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1"

    input:
        lineage: String

    output:
        record(
            lineage: lineage,
            models_path: file("models")
        )

    script:
    if (workflow.containerEngine == 'singularity') {
        log.warn("Running Helixer with Singularity is not recommended since you may encounter issues with permissions. " +
                 "Consider using Apptainer instead. See https://github.com/gglyptodon/helixer-docker for more information.")
    }
    """
    mkdir models 
    
    fetch_helixer_models.py \\
        --lineage $lineage \\
        --custom-path models
    """

}
