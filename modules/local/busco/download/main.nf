process BUSCO_DOWNLOADPROTEINS {
    tag "$lineage"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4137d65ab5b90d2ae4fa9d3e0e8294ddccc287e53ca653bb3c63b8fdb03e882f/data'
        : 'community.wave.seqera.io/library/busco:6.0.0--a9a1426105f81165'}"

    input:
    val lineage

    output:
    path "*.faa", emit: proteins
    tuple val("${task.process}"), val('busco'), eval('busco --version | sed "s/^BUSCO //"'),    topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        --download $lineage \\
        $args

    gzip -dc busco_downloads/lineages/${lineage}/refseq_db.faa.gz > ${lineage}.faa
    """
}
