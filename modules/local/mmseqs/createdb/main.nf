nextflow.enable.types = true

process MMSEQS_CREATEDB {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
        record(
            id: String,
            sequences: Set<Path>
        )

    stage:
        stageAs sequences, "tmp_input/*"

    output:
        record(
            id: id,
            custom_mmseqs_db: file("mmseqs_db/")
        )

    topic:
        tuple("${task.process}", 'mmseqs', eval('mmseqs version')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    # Ensure the input is uncompressed
    mkdir input_seqs
    cd input_seqs
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    prepared_sequences=\$(ls -1 input_seqs | tr '\n' ' ')

    mkdir -p mmseqs_db

    mmseqs \\
        createdb \\
        \${prepared_sequences} \\
        mmseqs_db/${prefix} \\
        ${args}

    """
}
