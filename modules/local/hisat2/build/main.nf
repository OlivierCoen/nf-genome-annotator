process HISAT2_BUILD {
    tag "${fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5bee187a0639f17702fc686a0244cfd32df6b2ad5786b97befdbacadc8ff03d/data' :
        'community.wave.seqera.io/library/hisat2:2.2.2--3dea1097582b823a'}"

    input:
    tuple val(meta), path(fasta), path(splicesites), path(exons)

    output:
    tuple val(meta), path("hisat2"), emit: index
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n 's/.*version \\([^ ]*\\).*/\\1/p'"), topic: versions


    script:
    def args = task.ext.args ?: ''
    def splice_site_arg = splicesites ? "--ss ${splicesites}" : ""
    def exon_arg = exons ? "--exon ${exons}" : ""
    """
    mkdir hisat2
    hisat2-build \\
        -p ${task.cpus} \\
        ${splice_site_arg} \\
        ${exon_arg} \\
        ${args} \\
        ${fasta} \\
        hisat2/${fasta.baseName}
    """
}
