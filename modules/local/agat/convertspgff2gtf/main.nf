nextflow.enable.types = true

process AGAT_CONVERTSPGFF2GTF {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abfa03eb1d5ee9a8f9aa056751126647577ee62ffac6a6ab84ca7a2184007380/data' :
        'community.wave.seqera.io/library/agat:1.7.0--9487e22276dbaaca' }"

    input:
        record(
            id: String,
            gff: Path
        )

    output:
        record(
            id: id,
            gtf: file("converted/*.gtf")
        )

    topic:
        tuple("${task.process}", 'agat', eval("agat_convert_sp_gff2gtf.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'

    script:
    def args = task.ext.args   ?: ''
    def is_compressed = gff.getExtension() == "gz" ? true : false
    def gff_name = is_compressed ? gff.getBaseName() : gff.name
    def gff_stem = gff_name.tokenize('.')[0..-2].join('.')
    """
    mkdir converted
    agat_convert_sp_gff2gtf.pl \\
        --gff ${gff} \\
        --output converted/${gff_stem}.gtf \\
        ${args}
    """
}
