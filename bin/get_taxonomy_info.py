#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
get_taxonomy_info.py - Get taxonomy ID and lineage for a given species name or Taxon ID
"""

import argparse
import logging

from taxonomy import Taxonomy

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

TAXID_OUTFILE = "taxid.txt"
TAXID_LINEAGE_OUTFILE = "taxid_lineage.txt"


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
    
    with open(TAXID_OUTFILE, "w") as fout:
        fout.write(str(species_taxonomy.taxid))

    with open(TAXID_LINEAGE_OUTFILE, "w") as fout:
        fout.write(','.join([str(taxid) for taxid in species_taxonomy.taxid_lineage]))

    logger.info("Done")
