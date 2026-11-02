#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import asyncio
import argparse

import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = 'orthodb.faa'

ORTHODB_API_URL = "https://data.orthodb.org/v13"
STOP_RETRY_AFTER_DELAY = 120
HTTPX_CLIENT_KWARGS = dict(
    follow_redirects=True,
    timeout=httpx.Timeout(600, connect=15)
)
MAX_PARALLEL_REQUESTS = 10

sem = asyncio.Semaphore(MAX_PARALLEL_REQUESTS) # semaphore, to limit the number of parallel requests

# ---------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--taxid-lineage', dest="taxid_lineage", type=str, required=True, 
        help='List of NCBI taxids representing the organism and its parents, from the most specific to the most general'
    )
    return parser.parse_args()
    

@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_orthodb(endpoint: str, params: dict | None = None):
    with httpx.Client(**HTTPX_CLIENT_KWARGS) as client:
        response = client.get(endpoint, params=params)
        response.raise_for_status()
    return response


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
async def fetch_species_fasta(client: httpx.AsyncClient, endpoint: str, params: dict):
    async with sem:
        response = await client.get(endpoint, params=params)
        response.raise_for_status()
        return response.text


def get_species_tree() -> dict:
    endpoint = f"{ORTHODB_API_URL}/tree"
    response = send_request_to_orthodb(endpoint)
    return response.json()['data']


def get_clade_species(clade_taxid: str) -> list[str]:
    endpoint = f"{ORTHODB_API_URL}/species"
    params = {'clade': clade_taxid}
    response = send_request_to_orthodb(endpoint, params)
    return [species_dict['taxon_id'] for species_dict in response.json()['data']]


def get_most_specific_taxid(taxids: list[str], tree: dict | list) -> dict | None:
    
    result = None
    if isinstance(tree, list):
        for subtree in tree:
           result = get_most_specific_taxid(taxids, subtree)
           if result is not None:
               break
    elif isinstance(tree, dict):
        # first, searching in children
        if 'children' in tree:
            result = get_most_specific_taxid(taxids, tree['children'])
        # if nothing in children, searching in the current dict
        if result is None:
            for taxid in taxids:
                if taxid == (tree_taxid := str(tree.get('key'))):
                    result = {
                        'taxid': tree_taxid,
                        'name': tree.get('name', 'Unknown')
                    }
                    break
                    
    return result


async def get_clade_fasta(species_taxids: list[str]):
    endpoint = f"{ORTHODB_API_URL}/fasta"
    async with httpx.AsyncClient(**HTTPX_CLIENT_KWARGS) as client:
        # launches all requests in parallel
        results = await asyncio.gather(
            *(fetch_species_fasta(client, endpoint, {'species': taxid}) for taxid in species_taxids),
            return_exceptions=True,
        )
    return results


# ---------------------------------------------------
# MAIN
# ---------------------------------------------------

def main():
    args = parse_args()

    lineage_taxids = args.taxid_lineage.split(',')

    logger.info("Fetching OrthoDB species tree")
    species_tree = get_species_tree()

    logger.info("Get the lowest level clade related to the provided taxid lineage")
    clade_result = get_most_specific_taxid(lineage_taxids, species_tree)
    
    if clade_result is None:
        raise ValueError(f"Could not find any matching taxid from lineage: {lineage_taxids}")

    clade_taxid = clade_result['taxid']
    clade_name = clade_result['name']

    logger.info(f"Found clade: {clade_name} (taxid {clade_taxid})")

    species_taxids = get_clade_species(clade_taxid)
    logger.info(f"Found {len(species_taxids)} species related to clade {clade_name}")

    outfile = f"{clade_name}.{OUTFILE_SUFFIX}"
    with open(outfile, 'w') as fout:
        logger.info(f"Getting fasta sequences related to {clade_name} species")
        results = asyncio.run(get_clade_fasta(species_taxids))
        for result, taxid in zip(results, species_taxids):
            if isinstance(result, Exception):
                logger.warning(f"Failed for species {taxid}: {result!r}")
                continue
            fout.write(result)

    logger.info("Done")

if __name__ == '__main__':
    main()
