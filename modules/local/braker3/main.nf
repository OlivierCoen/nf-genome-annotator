process BRAKER3 {
    tag "${meta.id}"
    label 'process_high'

    // Re. Conda from the BRAKER team:
    // Warning: installing GeneMark-ETP for BRAKER in conda environments has lead to multiple problems reported by users (Issues!).
    // We can not offer support for conda installations. Please use the singularity image instead.
    container "docker.io/teambraker/braker3:v3.0.7.5"

    input:
    tuple val(meta), path(fasta), path(proteins), path(bam)
    val species
    val busco_lineage

    output:
    tuple val(meta), path("$prefix/braker.gtf")         , emit: gtf
    tuple val(meta), path("$prefix/braker.codingseq")   , emit: cds
    tuple val(meta), path("$prefix/braker.aa")          , emit: aa
    tuple val(meta), path("$prefix/braker.log")         , emit: log
    tuple val(meta), path("$prefix/hintsfile.gff")      , emit: hintsfile   , optional: true
    tuple val(meta), path("$prefix/braker.gff3")        , emit: gff3        , optional: true
    tuple val(meta), path("$prefix/what-to-cite.txt")   , emit: citations

    tuple val("${task.process}"), val('braker3'), eval("braker.pl --version 2>/dev/null | sed 's/braker.pl version //'"),           topic: versions
    tuple val("${task.process}"), val('augustus'), eval("augustus --version |& sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p'"), topic: versions
    tuple val("${task.process}"), val('augustus-etp'), eval("gmetp.pl | sed -n 's/ETP version \\(.*\\)/\\1/p'"),                    topic: versions
    tuple val("${task.process}"), val('prothint'), eval("prothint.py --version | sed 's/prothint.py //1'"),                         topic: versions

    script:
    def args               = task.ext.args                   ?: ''
    prefix                 = task.ext.prefix                 ?: "${meta.id}"
    // The number of CPUs cannot exceed 48, otherwise BRAKER warns that it could create problems with GeneMark
    def nb_threads         = Math.min(48, task.cpus)
    def is_compressed      = fasta.getExtension() == "gz"    ? true : false
    def fasta_name         = is_compressed                   ? fasta.getBaseName() : fasta.name
    //def rna_ids     = rnaseq_sets_ids           ? "--rnaseq_sets_ids=$rnaseq_sets_ids"      : ''
    //def rna_dirs    = rnaseq_sets_dirs          ? "--rnaseq_sets_dirs=$rnaseq_sets_dirs"    : ''
    def bam_arg            = bam                             ? "--bam=$bam" : ''

    def prot_is_compressed = proteins && proteins.getExtension() == "gz" ? true : false
    def prot_fasta_name    = proteins && prot_is_compressed              ? proteins.getBaseName() : "fake"
    def prot_arg           = proteins ? "--prot_seq=$prot_fasta_name": ""

    //def hints       = hintsfile                 ? "--hints=$hintsfile"                      : ''
    def new_species        = args.contains('--species')      ? '' : '--species new_species'
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ -f $proteins -a "${prot_is_compressed}" == "true" ]; then
        gzip -c -d ${proteins} > ${prot_fasta_name}
    fi

    cp -r \$AUGUSTUS_CONFIG_PATH \\
        augustus_config

    chmod -R a+w \\
        augustus_config

    # keep only IDs in genome fasta header (and remove description)
    perl -p -e 's/^(>\\S+).*\$/\$1/' \\
        $fasta_name \\
        > ${prefix}.genome.masked.fasta

    braker.pl \\
        --genome ${prefix}.genome.masked.fasta \\
        --workingdir $prefix \\
        --AUGUSTUS_CONFIG_PATH "\$(pwd)/augustus_config" \\
        --AUGUSTUS_ab_initio \\
        --threads $task.cpus \\
        --busco_lineage=$busco_lineage
        $new_species \\
        $bam_arg \\
        $prot_arg \\
        $args
    """
}
