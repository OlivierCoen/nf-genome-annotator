nextflow.enable.types = true

process AGAT_CONVERTSPGFF2BED {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abfa03eb1d5ee9a8f9aa056751126647577ee62ffac6a6ab84ca7a2184007380/data' :
        'community.wave.seqera.io/library/agat:1.7.0--9487e22276dbaaca' }"

    input:
        record(
            id: String,
            gxf: Path
        )

    output:
        record(
            id: id,
            bed: file("*.bed")
        )

    topic:
        tuple("${task.process}", 'agat', eval("agat_convert_sp_gff2bed.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'

    script:
    def args   = task.ext.args   ?: ''
    def is_compressed = gxf.getExtension() == "gz" ? true : false
    def gxf_name = is_compressed ? gxf.getBaseName() : gxf.name
    def gxf_stem = gxf_name.tokenize('.')[0..-2].join('.')
    """
    agat_convert_sp_gff2bed.pl \\
        --gff ${gxf} \\
        --output ${gxf_stem}.bed \\
        ${args}
    """
}
