#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
# Adapted from script originally written by Tomas Bruna:
# https://github.com/tomasbruna/orthodb-clades/blob/3163c23f98417e20cddc5da42ce6e92ba0d18dc3/selectClade.py

"""
Get the species IDs corresponding to a species clade, optionnally excluding species / clades from given lists.
"""
                                 
import argparse
import csv
import sys
import logging

from pathlib import Path
from dataclasses import dataclass, field


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


SPECIES_IDS_OUTFILE = "species_ids.txt"


# ---------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clade', dest="clade", type=str, required=True, help='Clade to select proteins from')
    parser.add_argument('--levels', dest="levels_file", type=Path, required=True, help='Path to the levels OrthoDB file')
    parser.add_argument('--levels2species', dest="levels2species_file", type=Path, required=True, help='Path to the levels2species OrthoDB file')
    parser.add_argument('--species', dest="species_file", type=str, required=True, help='Path to the species OrthoDB file')
    parser.add_argument('--outfile', dest="outfile", type=str, required=True, help='Output file path')
    parser.add_argument('--exclude-clades', dest="excluded_clades", type=str, help='A comma separated set of clade names to exclude')
    parser.add_argument(
        '--exclude-species', dest="excluded_species",  type=str,
        help='A comma separated set of species names to exclude. When used, --species must be supplied.'
    )
    return parser.parse_args()


# ---------------------------------------------------
# SPECIES SELECTOR
# ---------------------------------------------------

@dataclass
class SpeciesSelector:

    clade: str
    levels_file: Path
    levels2species_file: Path
    species_file: Path
    excluded_clades: str | None
    excluded_species: str | None

    def get_species_ids(self) -> set[str]:

        clade_id = self.get_clade_id(self.clade)
        species_ids: set = self.parse_species_ids(clade_id)
    
        if self.excluded_clades:
            for clade_to_exclude in self.excluded_clades.split(','):
                excluded_clade_id = self.get_clade_id(clade_to_exclude)
                excluded_species = self.parse_species_ids(excluded_clade_id)
                species_ids.remove(excluded_species)
    
        if self.excluded_species:
            excluded_species = self.get_species_ids_from_names(self.excluded_species)
            species_ids -= excluded_species

        return species_ids
    

    def get_clade_id(self, clade: str):
        with open(self.levels_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if row[1] == clade:
                    return row[0]
        sys.exit(f'Error: Clade {clade} not found in {self.levels_file}')
    
    
    def parse_species_ids(self, clade_id: str) -> set:
        species_ids = set()
        with open(self.levels2species_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if clade_id in row[3].strip("{}").split(","):
                    species_ids.add(row[1])
        return species_ids
    
    
    def get_species_ids_from_names(self, species_names: str) -> set:
        species_names = {s.strip() for s in species_names.split(",")}
        excluded = set()
        with open(self.species_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if row[2] in species_names:
                    excluded.add(row[1])
        return excluded


# ---------------------------------------------------
# MAIN
# ---------------------------------------------------

def main():
    args = parse_args()

    # GET SPECIES IDS
    species_selector = SpeciesSelector(
        args.clade, 
        args.levels_file,
        args.levels2species_file,
        args.species_file,
        args.excluded_clades,
        args.excluded_species
    )
    
    species_ids = species_selector.get_species_ids()
    logger.info(f"Found {len(species_ids)} species for clade {args.clade}")
    
    logger.info(f"Writing species ids to {SPECIES_IDS_OUTFILE}")
    with open(SPECIES_IDS_OUTFILE, "w") as fout:
        # adding '^' for seqkit grep (as a regex pattern)
        fout.writelines([f"^{species_id}\n" for species_id in species_ids])


if __name__ == '__main__':
    main()
