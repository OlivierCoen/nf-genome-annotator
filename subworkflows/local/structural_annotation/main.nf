nextflow.enable.types = true

include { BRAKER           } from '../braker'
include { HELIXER          } from '../helixer'
include { METAEUK          } from '../metaeuk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    species: String
    fasta: Path
    orthodb_clade: String
    excluded_clades: List<String>
    excluded_species: List<String>
    mmseqs_db: String
    training_proteins: Set<Path>
    short_read_sorted_bams_bais: Set<Record>
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_input: Channel<Input>
    structural_annotator: String
    mmseqs_db: String
    skip_orthodb_download: Boolean
    skip_mmseqs_db_download: Boolean
    min_prot_db_seq_length: Integer

    main:

    if ( structural_annotator == "braker3" ) {

        BRAKER(
            ch_input,
            skip_orthodb_download,
            min_prot_db_seq_length
        )
        ch_annotated = BRAKER.out.annotated

    } else if ( structural_annotator == "helixer" ) {

        HELIXER( ch_input )
        ch_annotated = HELIXER.out.annotated
        
    } else if ( structural_annotator == "metaeuk" ) {

        METAEUK(
            ch_input,
            mmseqs_db,
            skip_mmseqs_db_download,
            min_prot_db_seq_length
        )
        ch_annotated = METAEUK.out.annotated
        
    }

    emit:
    ch_annotated

}
