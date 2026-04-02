process EGGNOGMAPPER_EMAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d12124094d33e21ac770e9447e8e3be6f9208e8bc49f84af236996eb61b243dc/data':
        'community.wave.seqera.io/library/eggnog-mapper:2.1.13--c99d97a9121734e6' }"

    input:
    tuple val(meta), path(fasta), path(gff)
    path(eggnog_data_dir)

    output:
    tuple val(meta), path("*.emapper.decorated.gff"),  emit: decorated_gff
    tuple val(meta), path("*.emapper.annotations"),    emit: annotations
    tuple val(meta), path("*.emapper.orthologs"),      emit: orthologs
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: seed_orthologs
    tuple val(meta), path("*.emapper.hits"),           emit: hits
    tuple val("${task.process}"), val('eggnog-mapper'), eval('emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//"'), topic: versions


    script:
    def args            = task.ext.args                 ?: ''
    def prefix          = task.ext.prefix               ?: "${meta.id}"
    def is_compressed   = fasta.extension == '.gz'      ? true                              : false
    def fasta_name      = is_compressed                 ? fasta.baseName                    : "$fasta"
    def dbmem           = task.memory.toMega() > 40000  ? '--dbmem'                         : ''
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    emapper.py \\
        --cpu ${task.cpus} \\
        -i ${fasta_name} \\
        --data_dir ${eggnog_data_dir} \\
        -m diamond \\
        --report_orthologs \\
        --decorate_gff ${gff} \\
        --output ${prefix} \\
        ${dbmem} \\
        $args
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.annotations
    touch ${prefix}.emapper.seed_orthologs
    touch ${prefix}.emapper.hits
    """
}
