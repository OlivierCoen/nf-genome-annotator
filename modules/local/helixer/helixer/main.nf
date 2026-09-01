nextflow.enable.types = true

process HELIXER_HELIXER {

    tag "${id} :: ${lineage}"
    label 'process_high'
    label 'process_gpu'

    // Helixer does not provide a conda package
    // Note: when using with docker, Nvidia docker toolkit should be installed
    // https://github.com/usadellab/helixer-docker#installing-docker
    container "docker.io/gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8_1"

    input:
        record(
            id: String,
            fasta: Path,
            lineage: String,
            models_path: Path
        )

    output:
        record(
            id: id,
            gff: file("*.gff3")
        )

    topic:
        tuple("${task.process}", 'helixer', eval("Helixer.py --version 2>&1 | grep 'Helixer.py' | cut -d' ' -f 2")) >> 'versions'

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${id}.helixer"
    // Warning if Singularity is being used
    if (workflow.containerEngine == 'singularity') {
        log.warn("Running Helixer with Singularity is not recommended since you may encounter issues with permissions. " +
                 "Consider using Apptainer instead. See https://github.com/gglyptodon/helixer-docker for more information.")
    }
    // Warning if not using GPUs
    if ( !workflow.profile.contains('gpu') ) {
        log.warn("Running Helixer without GPU(s) may take a long time. Running with GPU(s) can be activated by adding 'gpu' to profiles.")
    }
    """
    mkdir tmp
    
    Helixer.py \\
        --fasta-path $fasta \\
        --lineage $lineage \\
        --downloaded-model-path $models_path \\
        --gff-output-path ${prefix}.uncleaned.gff3 \\
        --temporary-dir tmp \\
        --deterministic \\
        ${args}

    # all gene and transcript names start with '_X' because no species name is provided here
    clean_helixer_gene_ids.py \\
        --gff ${prefix}.uncleaned.gff3 \\
        --out ${prefix}.gff3

    rm ${prefix}.uncleaned.gff3
    """

}
