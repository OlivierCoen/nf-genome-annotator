nextflow.enable.types = true

include { BRAKER3          } from '../braker3'
include { HELIXER          } from '../helixer'
include { METAEUK          } from '../metaeuk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Mapping {
    bam: Path
    bai: Path
}

record Input {
    id: String
    species: String
    fasta: Path
    orthodb_excluded_clades: List<String>
    orthodb_excluded_species: List<String>
    mmseqs_db: String
    proteins_fasta: Path
    short_read_mappings: List<Mapping>
    long_read_mappings: List<Mapping>
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_input: Channel<Input>
    structural_annotator: String
    mmseqs_db: String
    skip_mmseqs_db_download: Boolean

    main:

    if ( structural_annotator == "braker3" ) {

        ch_annotated = BRAKER3( ch_input )

    } else if ( structural_annotator == "helixer" ) {

        ch_annotated = HELIXER( ch_input )
        
    } else if ( structural_annotator == "metaeuk" ) {

        ch_annotated = METAEUK(
            ch_input,
            mmseqs_db,
            skip_mmseqs_db_download
        )
        
    }

    emit:
    ch_annotated

}
