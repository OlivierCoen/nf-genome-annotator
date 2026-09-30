#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
from dataclasses import dataclass

import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logger = logging.getLogger(__name__)


# Modern NCBI API
NCBI_TAXONOMY_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}
STOP_RETRY_AFTER_DELAY = 120


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

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

#####################################################
#####################################################
# CLASS
#####################################################
#####################################################

@dataclass
class Taxonomy:
    
    organism: str
    taxid: int
    common_name: str
    taxid_lineage: list[int]
    organism_lineage: list[str]
        
    def __str__(self) -> str:
        return f"{self.organism} [{self.common_name}] (taxid: {self.taxid})"

    @staticmethod
    def format_species_name(species: str) -> str:
        return species.replace("_", " ").capitalize().strip()


    @classmethod
    def get_taxonomy(cls, species: str) -> 'Taxonomy':
        logger.info(f"Querying taxid for species: {species}")
        formated_species = cls.format_species_name(species)
        metadata = cls.get_taxon_metadata(formated_species)
        if "lineage" not in metadata:
            raise KeyError(f"No lineage found for taxon {species}!")

        organism = metadata["organism_name"]
        taxid = int(metadata["tax_id"])
        lineage = metadata["lineage"]
        
        common_name = metadata.get("common_name", metadata["organism_name"])
        
        # from the most specific (species level or below) to the most general
        taxid_lineage = [taxid] + [int(taxid) for taxid in lineage[::-1]]

        organism_lineage = []
        for parent_taxid in lineage[::-1]:
            parent_metadata = cls.get_taxon_metadata(parent_taxid)
            organism_name = parent_metadata["organism_name"].lower()
            organism_lineage.append(organism_name)

        logger.info(f"Obtained taxonomy for: {organism}")
        return Taxonomy(
            organism=organism,
            taxid=taxid,
            taxid_lineage=taxid_lineage,
            common_name=common_name,
            organism_lineage=organism_lineage
        )


    @staticmethod
    def get_taxon_metadata(taxon: str | int) -> dict:
        result = send_request_to_ncbi_taxonomy(taxon)
        
        if len(result["taxonomy_nodes"]) > 1:
            raise ValueError(f"Multiple taxids for {taxon}")
            
        try:
            taxonomy_node = result["taxonomy_nodes"][0]
        except (IndexError, KeyError):
            msg = f"Could not find taxonomy results for taxon {taxon}"
            if "errors" in result:
                for error in result["errors"]:
                    msg += f"\n{error['reason']}"
            raise RuntimeError(msg)
            
        if "taxonomy" not in taxonomy_node:
            msg = f"Could not find taxonomy results for taxon {taxon}"
            if "errors" in taxonomy_node:
                for error in taxonomy_node["errors"]:
                    msg += f"\n{error['reason']}"
            raise RuntimeError(msg)
            
        return taxonomy_node["taxonomy"]


    def find_matching_lineage(self, organism_names: dict[str, str]) -> str | None:
        """
        Loops through the lineage of organism names (names of the current species and all successive parents, upward)
        and tries to find a matching name in the keys of organism_names. The first match is returned.

        Args:
            organism_names (dict[str, str]): A dictionary mapping organism names to lineage identifiers.

        Returns:
            str | None: The lineage identifier of the first matching organism name.
        """
        logger.info(f"Finding matching lineage for {self.organism} among {len(organism_names)} names")
        matching_lineage = None
    
        # looping from the most specific (species level or below) to the most general
        for parent_organism in self.organism_lineage:
            if parent_organism in organism_names:
                logger.info(f"Found lineage match: {parent_organism} -> {organism_names[parent_organism]}")
                matching_lineage = organism_names[parent_organism]
                break
    
        if matching_lineage is None:
            logger.warning(
                f"No lineage match found for species {self.organism} "
                f"and lineage {' '.join([str(l) for l in self.organism_lineage])}."
            )
    
        return matching_lineage


    def get_tool_lineage(self, tool: str, model_lineages: dict) -> str | None:
        """
        Wrapper function that finds the best model given a lineage
        """
        # searching for a candidate lineage
        # looping from the most general to the most specific (species level or below)
        candidate_lineages = []
        for taxid in self.taxid_lineage:
            for lineage, lineage_dict in model_lineages.items():
                if taxid == lineage_dict["included"]:
                    candidate_lineages.append(lineage)
                    
        if candidate_lineages:
            logger.info(f"{tool} candidate lineages: {candidate_lineages}")
        else:
            logger.warning(f"No candidate lineage found for {tool}")
            return None
            
        # among the candidate lineages, searching for the ones that are excluded
        for taxid in self.taxid_lineage:
            for lineage in candidate_lineages:
                if taxid == model_lineages[lineage]["excluded"]:
                    del candidate_lineages[lineage]
    
        if candidate_lineages:
            if len(candidate_lineages) > 1:
                logger.warning(f"Multiple {tool} lineages found for taxid lineage: {self.taxid_lineage}. Taking the first one.")
                return None
            found_lineage = candidate_lineages[0]
            logger.info(f"Found {tool} lineage: {found_lineage}")
            return found_lineage
        else:
            logger.warning(f"All candidate lineages were excluded for {tool}")
            return None