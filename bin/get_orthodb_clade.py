#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
get_orthodb_clade.py - Get OrthoDB clade corresponding to a specific species, among a list of available clades.
Write OrthoDB clade to orthodb_clade.txt 

Usage: get_orthodb_clade.py --species <species name> --orthodb-clades <file containing list of clades, one per line>
"""

import argparse
import logging
from pathlib import Path

from taxonomy import Taxonomy

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

ORTHODB_CLADES_OUTFILE = "orthodb_clade.txt"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--species", 
        type=str, 
        required=True,
        help="Species name (or Taxon ID)"
    )
    parser.add_argument(
        "--orthodb-clades",
        dest="orthodb_clades",
        type=Path,
        required=True,
        help="Path to file containing OrthoDB clades"
    )
    return parser.parse_args()


def parse_orthodb_clades(file_path: Path) -> dict[str, str]:
    with open(file_path, "r") as fin:
        lines = [line.strip() for line in fin.readlines()]
    clades = [line.split("\t")[1] for line in lines]
    return {clade.lower(): clade for clade in clades}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    species_taxonomy = Taxonomy.get_taxonomy(args.species)
    
    lineages_to_orthodb_clades = parse_orthodb_clades(args.orthodb_clades)
    logger.info(f"Parsed {len(lineages_to_orthodb_clades)} orthodb clades")

    orthodb_clade = species_taxonomy.find_matching_lineage(lineages_to_orthodb_clades)

    with open(ORTHODB_CLADES_OUTFILE, "w") as fout:
        fout.write(orthodb_clade)

    logger.info("Done")
