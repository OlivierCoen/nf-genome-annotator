#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
get_busco_lineage.py - Get BUSCO lineage corresponding to a specific species, among a list of available lineages.
Write BUSCO lineage to busco_lineage.txt

Usage: get_busco_lineage.py --species <species name> --busco-datasets <output of busco --list-datasets>
"""

import argparse
import logging
from pathlib import Path
import re

from taxonomy import Taxonomy

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

BUSCO_LINEAGE_OUTFILE = "busco_lineage.txt"

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
        "--busco-datasets",
        dest="busco_datasets",
        type=Path,
        required=True,
        help="Path to file containing BUSCO datasets"
    )
    return parser.parse_args()


def parse_busco_datasets(file_path: Path) -> dict[str, str]:
    with open(file_path, "r") as fin:
        lines = [line.strip() for line in fin.readlines()]
    busco_datasets = {}
    for line in lines:
        if match := re.findall(r"\w+_odb\d+\.?\d+?", line):
            dataset = match[0]
            lineage = dataset.split("_")[0]
            busco_datasets[lineage] = dataset
        else:
            logger.warning(f"Could not parse line: {line}")
    return busco_datasets


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    species_taxonomy = Taxonomy.get_taxonomy(args.species)
    
    lineages_to_busco_datasets = parse_busco_datasets(args.busco_datasets)
    logger.info(f"Parsed {len(lineages_to_busco_datasets)} busco lineages")

    busco_lineage = species_taxonomy.find_matching_lineage(lineages_to_busco_datasets)

    with open(BUSCO_LINEAGE_OUTFILE, "w") as fout:
        fout.write(str(busco_lineage))

    logger.info("Done")
