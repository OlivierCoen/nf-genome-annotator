nextflow.enable.types = true

process AGAT_SPCOMPLEMENTANNOTATIONS {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abfa03eb1d5ee9a8f9aa056751126647577ee62ffac6a6ab84ca7a2184007380/data':
        'community.wave.seqera.io/library/agat:1.7.0--9487e22276dbaaca' }"

    input:
        record(
            id: String, 
            ref_gff: Path,
            other_gff: Path
        )

    output:
        record(
            id: id,
            gff: file("*_complemented.gff")
        )

    topic:
        tuple("${task.process}", 'agat', eval("sp_complement_annotations.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p'")) >> 'versions'
        
    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    agat_sp_complement_annotations.pl \\
        --ref $ref_gff \\
        --add $other_gff \\
       ${args} \\
        --output ${prefix}_complemented.gff
    """
}
