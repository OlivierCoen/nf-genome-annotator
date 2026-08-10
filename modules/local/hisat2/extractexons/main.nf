nextflow.enable.types = true

process HISAT2_EXTRACTEXONS {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5bee187a0639f17702fc686a0244cfd32df6b2ad5786b97befdbacadc8ff03d/data' :
        'community.wave.seqera.io/library/hisat2:2.2.2--3dea1097582b823a'}"

    input:
        record(
            sample_id: String,
            gtf: Path
        )

    output:
        record(
            sample_id: sample_id,
            exons: file("*.exons.txt")
        )

    topic:
        tuple("${task.process}", 'hisat2', eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2')) >> 'versions'
        tuple("${task.process}", 'samtools', eval("samtools --version | sed -n '1s/samtools //p'"))              >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_exons.py \\
        $args \\
        $gtf \\
        > ${sample_id}.exons.txt
    """
}
