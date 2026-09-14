nextflow.enable.types = true

process HISAT2_BUILD {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5bee187a0639f17702fc686a0244cfd32df6b2ad5786b97befdbacadc8ff03d/data' :
        'community.wave.seqera.io/library/hisat2:2.2.2--3dea1097582b823a'}"

    input:
        record(
            id: String,
            fasta: Path,
            splice_sites: Path?,
            exons: Path?
        )

    output:
        record(
            id: id,
            index: file("hisat2", type: 'dir')
        )

    topic:
        tuple("${task.process}", 'hisat2', eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def splice_site_arg = splice_sites ? "--ss ${splice_sites}" : ""
    def exon_arg = exons ? "--exon ${exons}" : ""
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta.name
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi
    
    mkdir hisat2
    hisat2-build \\
        -p ${task.cpus} \\
        ${splice_site_arg} \\
        ${exon_arg} \\
        ${args} \\
        ${fasta_name} \\
        hisat2/${id}

    if [ "${is_compressed}" == "true" ]; then
        rm ${fasta_name}
    fi
    """
}
