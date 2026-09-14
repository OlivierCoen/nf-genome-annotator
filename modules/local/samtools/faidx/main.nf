nextflow.enable.types = true

process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
        record(id: String, fasta: Path)

    output:
        record(
            id: id,
            fai: file("*.fai")
        )

    topic:
        tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta.name
    """
    # uncompressing, because samtools does not like gzipped fasta (only bgzipped)
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi
    
    samtools \\
        faidx \\
         ${fasta_name} \\
        ${args}

    if [ "${is_compressed}" == "true" ]; then
        echo "Removing ${fasta_name}"
        rm ${fasta_name}
    fi
    """
}
