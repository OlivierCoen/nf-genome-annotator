#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from multiprocessing import Value
from pathlib import Path
import re


import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)
from dataclasses import dataclass

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_TAXONOMY_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}
STOP_RETRY_AFTER_DELAY = 120

TAXID_OUTFILE = "found_taxid.txt"
BUSCO_LINEAGE_OUTFILE = "found_busco_lineage.txt"
ORTHODB_CLADES_OUTFILE = "found_orthodb_clade.txt"
HELIXER_LINEAGE_OUTFILE = "found_helixer_lineage.txt"

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
# CLASSES
#####################################################
#####################################################

@dataclass
class Taxonomy:
    species: str
    taxid: str | int
    common_name: str
    lineage: list[int]
    busco_lineage: str | None
    orthodb_clade: str | None
    helixer_lineage: str | None

    def __str__(self) -> str:
        return f"{self.species} [{self.common_name}] (taxid: {self.taxid})"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get taxonomy ID for a given species name or Taxon ID"
    )
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
    parser.add_argument(
        "--orthodb-clades",
        dest="orthodb_clades",
        type=Path,
        required=True,
        help="Path to file containing OrthoDB clades"
    )
    return parser.parse_args()


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str | int):
    logger.info(f"Sending POST request to {NCBI_TAXONOMY_API_URL}")
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    response = httpx.post(NCBI_TAXONOMY_API_URL, headers=NCBI_API_HEADERS, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_esearch_query(query: str, database: str, ncbi_api_key: str | None):
    """
    Query NCBI's db with Esearch API
    """
    params = dict(db=database, term=query, retmax=ESEARCH_RETMAX)
    if ncbi_api_key is not None:
        params["api_key"] = ncbi_api_key
    response = httpx.get(ESEARCH_BASE_URL, params=params)
    response.raise_for_status()
    return response.text


def get_taxon_metadata(taxid: str | int) -> dict:
    result = send_request_to_ncbi_taxonomy(taxid)
    
    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for {taxid}")
        
    try:
        taxonomy_node = result["taxonomy_nodes"][0]
    except (IndexError, KeyError):
        msg = f"Could not find taxonomy results for taxid {taxid}"
        if "errors" in result:
            for error in result["errors"]:
                msg += f"\n{error['reason']}"
        raise RuntimeError(msg)
        
    if "taxonomy" not in taxonomy_node:
        msg = f"Could not find taxonomy results for species {species}"
        if "errors" in taxonomy_node:
            for error in taxonomy_node["errors"]:
                msg += f"\n{error['reason']}"
        raise RuntimeError(msg)
        
    return taxonomy_node["taxonomy"]


def get_taxonomy(
    species: str,
    lineages_to_busco_datasets: dict[str, str],
    lineages_to_orthodb_clades: dict[str, str],
) -> Taxonomy:
    
    metadata = get_taxon_metadata(species)

    # trying to find a lineage match in the busco datasets and orthodb clades
    busco_lineage = None
    orthodb_clade = None

    if "lineage" not in metadata:
        raise KeyError(f"No lineage found for taxid {species}!")
    
    lineage = metadata["lineage"]
    # looping from the most specific (species level or below) to the most general
    for parent_taxid in lineage[::-1]:
        parent_metadata = get_taxon_metadata(parent_taxid)
        parent_organism_name = parent_metadata["organism_name"].lower()
        
        if parent_organism_name in lineages_to_busco_datasets:
            logger.info(f"Found lineage match for BUSCO: {parent_organism_name} -> {lineages_to_busco_datasets[parent_organism_name]}")
            busco_lineage = lineages_to_busco_datasets[parent_organism_name]

        if parent_organism_name in lineages_to_orthodb_clades:
            logger.info(f"Found lineage match for OrthoDB: {parent_organism_name} -> {lineages_to_orthodb_clades[parent_organism_name]}")
            orthodb_clade = lineages_to_orthodb_clades[parent_organism_name]

        # if both were found, break out of the loop
        if busco_lineage and orthodb_clade:
            break

    if busco_lineage is None:
        logger.warning(f"No lineage match found for species {species} and lineage {' '.join([str(l) for l in lineage[::-1]])}.")

    if orthodb_clade is None:
        logger.warning(f"No orthodb clade match found for species {species} and lineage {' '.join([str(l) for l in lineage[::-1]])}.")

    helixer_lineage = get_helixer_lineage(lineage)

    common_name = metadata.get("common_name", metadata["organism_name"])

    return Taxonomy(
        species=metadata["organism_name"],
        taxid=int(metadata["tax_id"]),
        common_name=common_name,
        lineage=lineage,
        busco_lineage=busco_lineage,
        orthodb_clade=orthodb_clade,
        helixer_lineage=helixer_lineage
    )


def format_species_name(species: str) -> str:
    return species.replace("_", " ").capitalize().strip()


def parse_busco_datasets(file_path: Path) -> dict[str, str]:
    with open(file_path, "r") as fin:
        lines = [line.strip() for line in fin.readlines()]
    busco_datasets = {}
    for line in lines:
        if match := re.findall(r"\w+_odb\d+\.?\d+?", line):
            dataset = match[0]
            lineage = dataset.split("_")[0]
            busco_datasets[lineage] = lineage
        else:
            logger.warning(f"Could not parse line: {line}")
    return busco_datasets


def parse_orthodb_clades(file_path: Path) -> dict[str, str]:
    with open(file_path, "r") as fin:
        lines = [line.strip() for line in fin.readlines()]
    clades = [line.split("\t")[1] for line in lines]
    return {clade.lower(): clade for clade in clades}


def get_helixer_lineage(taxid_lineage: list[int]) -> str | None:
    # searching for a candidate lineage
    # looping from the most general to the most specific (species level or below)
    candidate_lineages = []
    for taxid in taxid_lineage:
        for lineage, lineage_dict in HELIXER_LINEAGES.items():
            if taxid == lineage_dict["included"]:
                candidate_lineages.append(lineage)
                
    if candidate_lineages:
        logger.info(f"Helixer candidate lineages: {candidate_lineages}")
    else:
        logger.info("No candidate lineage found for Helixer")
        return None
        
    # among the candidate lineages, searching for the ones that are excluded
    for taxid in taxid_lineage:
        for lineage in candidate_lineages:
            if taxid == HELIXER_LINEAGES[lineage]["excluded"]:
                del candidate_lineages[lineage]

    if candidate_lineages:
        if len(candidate_lineages) > 1:
            raise ValueError(f"Multiple Helixer lineages found for taxid lineage: {taxid_lineage}")
        found_lineage = candidate_lineages[0]
        logger.info(f"Found Helixer lineage: {found_lineage}")
        return found_lineage[0]
    else:
        logger.info("All candidate lineages were excluded")
        return None

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    species = format_species_name(args.species)

    lineages_to_busco_datasets = parse_busco_datasets(args.busco_datasets)
    logger.info(f"Parsed {len(lineages_to_busco_datasets)} busco lineages")

    lineages_to_orthodb_clades = parse_orthodb_clades(args.orthodb_clades)
    logger.info(f"Parsed {len(lineages_to_orthodb_clades)} orthodb clades")

    logger.info(f"Querying taxid for species: {species}")
    taxonomy = get_taxonomy(
        species,
        lineages_to_busco_datasets,
        lineages_to_orthodb_clades
    )
    logger.info(f"Obtained taxonomy for: {taxonomy}")

    with open(TAXID_OUTFILE, "w") as fout:
        fout.write(str(taxonomy.taxid))

    with open(BUSCO_LINEAGE_OUTFILE, "w") as fout:
        fout.write(str(taxonomy.busco_lineage))

    with open(ORTHODB_CLADES_OUTFILE, "w") as fout:
        fout.write(str(taxonomy.orthodb_clade))

    with open(HELIXER_LINEAGE_OUTFILE, "w") as fout:
        fout.write(str(taxonomy.helixer_lineage))

    logger.info("Done")
