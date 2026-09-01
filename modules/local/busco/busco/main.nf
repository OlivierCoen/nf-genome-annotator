nextflow.enable.types = true

process BUSCO_BUSCO {
    
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f67e816ab2f7ccc9cb2d40874dea1e2e1a8e88ef6a44750b66c0ee55fe8de6c/data'
            : 'community.wave.seqera.io/library/busco:6.1.0--0e40710a525d8d44'}"
    // Note: one test had to be disabled when switching to Busco 6.0.0, cf https://github.com/nf-core/modules/pull/8781/files
    // Try to restore it when upgrading Busco to a later version

    input:
        record(
            id: String,
            fasta: Iterable<Path>,
            lineage: String,
            download_path: Path
        )
        mode: String

    stage:
        stageAs fasta, 'tmp_input/*'

    topic:
        tuple(id, files("short_summaries/*.txt"))                                      >> 'busco_multiqc'
        tuple('busco', id, file("*-busco.batch_summary.txt"))                          >> 'additional_results'
        tuple('busco', id, file('*-busco.log'))                                        >> 'logs'
        tuple("${task.process}", 'busco', eval('busco --version | sed "s/^BUSCO //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}.${mode}"
    def intermediate_files = [
        './*-busco/*/auto_lineage',
        './*-busco/*/**/{miniprot,hmmer,.bbtools}_output',
        './*-busco/*/prodigal_output/predicted_genes/tmp/',
    ]
    def bbtools_memory_preferred = task.memory * 0.25
    def bbtools_memory_minimum = 120.Mb
    def bbtools_memory = bbtools_memory_preferred > bbtools_memory_minimum ? "${bbtools_memory_preferred.toGiga()}g" : "${bbtools_memory_minimum.toMega()}m"
    """
    export BUSCO_BBTOOLS_MEMORY=${bbtools_memory}
    
    # Fix Augustus for Apptainer
    ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
    set +u
    if [ -z "\${AUGUSTUS_CONFIG_PATH}" ] && [ -f "\${ENV_AUGUSTUS}" ]; then
        source "\${ENV_AUGUSTUS}"
    fi
    set -u

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix}-busco \\
        --mode ${mode} \\
        --lineage_dataset $lineage \\
        --download_path ${download_path} \\
        ${args}

    # clean up
    rm -rf "\$INPUT_SEQS"
    rm -fr ${intermediate_files.join(' ')}

    # Move files to avoid staging/publishing issues
    mv ${prefix}-busco/batch_summary.txt ${prefix}-busco.batch_summary.txt

    mkdir short_summaries
    mv ${prefix}-busco/*/short_summary.*.txt short_summaries/${prefix}.busco.txt

    mv ${prefix}-busco/logs/busco.log ${prefix}-busco.log

    if grep 'Run failed; check logs' ${prefix}-busco.batch_summary.txt > /dev/null
    then
        echo "Busco run failed"
        exit 1
    fi
    """
}
