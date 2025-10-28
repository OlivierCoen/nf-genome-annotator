process REPEATMASKER_REPEATMASKER {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb2c806e7a34fc52eabb8809fa9f4e0e34117e397162a516cbb80454e50b2e72/data':
        'community.wave.seqera.io/library/repeatmasker:4.2.1--94dd45a91fa85e83' }"

    input:
    tuple val(meta), path(fasta), path(lib)

    output:
    tuple val(meta), path("${prefix}.masked"),                                                                        emit: masked
    tuple val(meta), path("${prefix}.out"),                                                                           emit: out
    tuple val(meta), path("${prefix}.tbl"),                                                                           emit: tbl
    tuple val(meta), path("${prefix}.gff"),                                                                           emit: gff, optional: true
    tuple val("${task.process}"), val('repeatmasker'), eval("RepeatMasker -v | sed 's/RepeatMasker version //1'"),    topic: versions

    script:
    def args    = task.ext.args     ?: ''
    prefix      = task.ext.prefix   ?: "${meta.id}"
    def lib_arg = lib               ? "-lib $lib"   : ''

    def out_fasta    = fasta.getBaseName(fasta.name.endsWith('.gz') ? 1 : 0)
    def fasta_gz_cmd = fasta.name.endsWith('.gz') ? "gunzip -c ${fasta} > ${out_fasta}" : ""

    """
    ${fasta_gz_cmd}
    RepeatMasker \\
        $lib_arg \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        ${args} \\
        ${out_fasta}

    mv $prefix/${out_fasta}.masked  ${prefix}.masked
    mv $prefix/${out_fasta}.out     ${prefix}.out
    mv $prefix/${out_fasta}.tbl     ${prefix}.tbl
    mv $prefix/${out_fasta}.out.gff ${prefix}.gff       || echo "GFF is not produced"
    """

    stub:
    prefix          = task.ext.prefix       ?: "${meta.id}"
    def args        = task.ext.args         ?: ''
    def touch_gff   = args.contains('-gff') ? "touch ${prefix}.gff" : ''

    """
    touch ${prefix}.masked
    touch ${prefix}.out
    touch ${prefix}.tbl
    $touch_gff

    """
}
