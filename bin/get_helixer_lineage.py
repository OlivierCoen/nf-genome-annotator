#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
get_helixer_lineage.py - Get Helixer lineage corresponding to a specific species. The available lineages are hardcoded in the script.
Write Helixer lineage to helixer_lineage.txt

Usage: get_helixer_lineage.py --species <species name>
"""

import argparse
import logging

from taxonomy import Taxonomy

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

HELIXER_LINEAGE_OUTFILE = "helixer_lineage.txt"

HELIXER_LINEAGES = {
    "land_plant": {
        "included": 3193, #Embryophyta
        "excluded": None
    },
    "vertebrate": {
        "included": 7742,
        "excluded": None
    },
    "invertebrate": { # everything in metazoans except for vertebrates
        "included": 33208,
        "excluded": 7742
    },
    "fungi":  {
        "included": 4751,
        "excluded": None
    }
}

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
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    
    species_taxonomy = Taxonomy.get_taxonomy(args.species)
    
    helixer_lineage = species_taxonomy.get_tool_lineage('Helixer', HELIXER_LINEAGES)

    with open(HELIXER_LINEAGE_OUTFILE, "w") as fout:
        fout.write(str(helixer_lineage))

    logger.info("Done")
