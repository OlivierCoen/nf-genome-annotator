nextflow.enable.types = true

process EGGNOGMAPPER_EMAPPER {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d1/d12124094d33e21ac770e9447e8e3be6f9208e8bc49f84af236996eb61b243dc/data':
        'community.wave.seqera.io/library/eggnog-mapper:2.1.13--c99d97a9121734e6' }"

    input:
        record(
            id: String, 
            fasta: Path,
            gff: Path
        )
        eggnog_data_dir: Path

    output:
        record(
            id: id,
            gff: file("*.emapper.decorated.gff")
        )

    topic:
        tuple('eggnog-mapper', id, file("*.emapper.annotations"))    >> 'additional_results'
        tuple('eggnog-mapper', id, file("*.emapper.orthologs"))      >> 'additional_results'
        tuple('eggnog-mapper', id, file("*.emapper.seed_orthologs")) >> 'additional_results'
        tuple('eggnog-mapper', id, file("*.emapper.hits"))           >> 'additional_results'
        tuple("${task.process}", 'eggnog-mapper', eval('emapper.py --version | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//"')) >> 'versions'

    script:
    def args            = task.ext.args                 ?: ''
    def prefix          = task.ext.prefix               ?: "$id"
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
        --report_orthologs \\
        --decorate_gff ${gff} \\
        --output ${prefix} \\
        ${dbmem} \\
        ${args}
    """
}
