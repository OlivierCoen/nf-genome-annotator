process BUSCO_BUSCO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/963bad66c10646cf0adb1967cc462ad04d02789ddbfae4fbb94182291dbddf8c/data'
            : 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'}"
    // Note: one test had to be disabled when switching to Busco 6.0.0, cf https://github.com/nf-core/modules/pull/8781/files
    // Try to restore it when upgrading Busco to a later version

    input:
    tuple val(meta), path(fasta, stageAs: 'tmp_input/*'), path(busco_lineages_path)
    // One of genome, proteins, or transcriptome
    val mode

    output:
    tuple val(meta), path("*-busco.batch_summary.txt"),                                      emit: batch_summary
    tuple val(meta), path("*-busco.log"),                                                    emit: log, optional: true
    tuple val(meta),path("short_summaries/*.txt"),                                           topic: mqc_busco_short_summaries_txt
    tuple val("${task.process}"), val('busco'), eval('busco --version | sed "s/^BUSCO //"'), topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = mode == 'genome' ? "${meta.id}.genome" : meta.final_annotation ? "${meta.id}.final_proteome" : "${fasta.baseName}.intermediate_proteome"
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
        --download_path ${busco_lineages_path} \\
        ${args}

    # clean up
    rm -rf "\$INPUT_SEQS"
    rm -fr ${intermediate_files.join(' ')}
    
    # find and remove broken symlinks from the cleanup
    find . -xtype l -delete

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
