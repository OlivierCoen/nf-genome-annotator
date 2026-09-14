#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import json
import logging
from enum import Enum
import xml.etree.ElementTree as ET

import requests
import urllib3
import xmltodict
from tenacity import (
    before_sleep_log,
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# E-UTILITIES OL API
ESEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_BASE_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db}&id={id}"
)
ESEARCH_RETMAX = 1000000000  # max retmax that worked
CHUNKSIZE = 400 # 500 is already too much, but 400 seems to work fine

SHORT_READ_SRA_IDS_OUTFILE = "sra_ids.short_read.txt"
LONG_READ_SRA_IDS_OUTFILE = "sra_ids.long_read.txt"
SHORT_READ_METADATA_OUTFILE = "sra_metadata.short_read.json"
LONG_READ_METADATA_OUTFILE = "sra_metadata.long_read.json"

#####################################################
# ESEARCH QUERY TERMS
#####################################################

RNASEQ_SELECTION_TEMRS = [
    '"cdna oligo dt"[Selection]',
    '"cdna randompriming"[Selection]',
    '"oligo dt"[Selection]',
    '"polya"[Selection]',
    '"random"[Selection]',
    '"inverse rrna"[Selection]'
]

BASE_TAXID_TERM = 'txid{taxid}[Taxonomy]'

BASE_QUERY_ITEMS = [
    '"biomol rna"[Properties]',
    '"filetype fastq"[Properties]',
    f'("rna seq"[Strategy] OR ("other"[Strategy] AND {" OR ".join(RNASEQ_SELECTION_TEMRS)}))'
]

class Platform(Enum):
    SHORT_READ = "short_reads"
    SHORT_READ_FALLBACK = "short_reads_fallback"
    LONG_READ = "long_reads"


PLATFORM_TERM = {
    Platform.SHORT_READ: '"illumina"[Platform]',
    Platform.SHORT_READ_FALLBACK: '("bgiseq"[Platform] OR "ion_torrent"[Platform])',
    Platform.LONG_READ: '("pacbio_smrt"[Platform] OR "oxford_nanopore"[Platform])'
}

PAIRED_LAYOUT_TERM = '"paired"[Layout]'


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument("--taxid", type=int, required=True, help="NCBI taxon ID")
    parser.add_argument("--max-short-reads", required=True,dest="max_short_reads_data", type=int, help="Maximum number of short read sequencing data to fetch")
    parser.add_argument("--max-long-reads", required=True, dest="max_long_reads_data", type=int, help="Maximum number of short read sequencing data to fetch")
    parser.add_argument("--paired-only", dest="paired_only", action='store_true', help="For short reads, fetch only paired-end sequencing data")
    parser.add_argument("--max-size", dest="max_size", type=str, help="Maximum size (in Mb / Gb) for an experiment. Must end with 'Mb' or 'Gb'. Example: 10Mb / 2Gb")
    
    return parser.parse_args()
    


class RateLimitException(Exception):
    pass


@retry(
    retry=retry_if_exception_type(
        (
            RateLimitException,
            urllib3.exceptions.ReadTimeoutError,
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
        )
    ),
    stop=stop_after_delay(3600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_esearch_query(query: str, database: str):
    """
    Query NCBI's db with Esearch API
    """
    params = {
        'db': database,
        'term': query,
        'retmax': ESEARCH_RETMAX,
    }
    response = requests.get(ESEARCH_BASE_URL, params=params)
    if response.status_code == 429:
        raise RateLimitException("Rate limit exceeded")
    response.raise_for_status()
    return response.text


@retry(
    retry=retry_if_exception_type(
        (
            RateLimitException,
            urllib3.exceptions.ReadTimeoutError,
            requests.exceptions.ConnectionError,
        )
    ),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.DEBUG),
)
def send_efetch_query(query_ids: list[str], database: str) -> str:
    """
    Query NCBI's db with Efetch API
    """
    query_ids_str = ",".join(query_ids)
    url = EFETCH_BASE_URL.format(db=database, id=query_ids_str)
    response = requests.get(url)
    if response.status_code == 429:
        raise RateLimitException("Rate limit exceeded")
    response.raise_for_status()
    return response.text


def parse_ids_from_xml(xml_string: str) -> list[str]:
    """
    Parse XML string and get text in <Id>...</Id> blocks
    :param xml_string:
    :return: list of Ids
    """
    # parsing XML returned by API
    root = ET.fromstring(xml_string)
    # species taxon IDs are contained in <Id>...</Id> blocks
    return [
        id_element.text 
        for id_element in root.findall(".//Id") 
        if id_element.text is not None
    ]


def parse_sra_accessions_from_xml(xml_string: str) -> list[dict]:
    try:
        xml_dict = xmltodict.parse(xml_string)
    except Exception as e:
        logger.error(f"Failed to parse XML: {e}")
        return []
    experiments = xml_dict["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]
    # if type is dict, turn into a list of dict (of one element)
    if isinstance(experiments, dict):
        experiments = [experiments]
    return [
        experiment_dict
        for experiment_dict in experiments
        if isinstance(experiment_dict, dict)
    ]


def get_esearch_query(taxid: int, platform: Platform, paired_only: bool) -> str:
    """
    Constructs an Esearch query string for the given taxid, platform, and paired_only flag.
    """
    taxid_term = BASE_TAXID_TERM.format(taxid=taxid)
    query_items = [taxid_term] + BASE_QUERY_ITEMS
    query_items.append(PLATFORM_TERM[platform])
    if platform in [Platform.SHORT_READ, Platform.SHORT_READ_FALLBACK] and paired_only:
        query_items.append(PAIRED_LAYOUT_TERM)
    return " AND ".join(query_items)


def get_sra_ids(taxid: int, platform: Platform, paired_only: bool = False) -> list[str]:
    """
    Get list of SRA experiment IDs given a NCBI taxonomy ID.
    Esearch returns experiments corresponding to all the NCBI TaxID underneath 
    the provided organism name
    :param organism_name: str
    :return: list of SRA experiment IDs
    """
    query = get_esearch_query(taxid, platform, paired_only)
    logger.info(f"Sending Esearch query: {query}")
    xml_string = send_esearch_query(query=query, database="sra")
    return parse_ids_from_xml(xml_string)


def fetch_sra_experiments_for_ids(sra_uids: list[str]) -> list[dict]:
    """
    Get list of SRA experiment IDs given a NCBI taxonomy ID
    :param taxid:
    :return: list of SRA experiment IDs
    """
    xml_string = send_efetch_query(query_ids=sra_uids, database="sra")
    return parse_sra_accessions_from_xml(xml_string)


def fetch_sra_experiments_chunks(sra_uids: list[str]) -> list[dict]:
    """
    Fetch SRA experiment metadata for a list of SRA experiment IDs
    :param sra_uids: list of SRA experiment IDs
    :return: list of dictionaries, each one containing experiment metadata for a specific SRA ID
    """
    try:
        return fetch_sra_experiments_for_ids(sra_uids)
    except requests.exceptions.HTTPError as e:
        # handling known errors
        if e.response.status_code == 414:  # request too long
            logger.warning(
                f"Too long request URI with {len(sra_uids)} SRA IDs: dividing in 2 and sending 2 separates requests"
            )
            # dividing into 2 chunks and launching new requests recursively
            sra_uids_part_1 = sra_uids[: len(sra_uids) // 2]
            sra_uids_part_2 = sra_uids[len(sra_uids) // 2 :]
            return fetch_sra_experiments_chunks(sra_uids_part_1) + fetch_sra_experiments_chunks(
                sra_uids_part_2
            )
        else:
            raise

def fetch_sra_experiments(sra_uids: list[str], read_type: str):
    experiments = []
    sra_ids_chunks = [
        short_read_sra_raw_ids[i : i + CHUNKSIZE]
        for i in range(0, len(short_read_sra_raw_ids), CHUNKSIZE)
    ]
    logger.info(
        f"Fetching sra experiment metadata for each {read_type} SRA experiment ID ({len(sra_ids_chunks)} chunks in total)"
    )
    for i, sra_ids_chunk in enumerate(sra_ids_chunks):
        logger.info(f"Fetching chunk {i + 1}")
        result = fetch_sra_experiments_chunks(sra_ids_chunk)
        experiments += result
    return experiments


def parse_max_size(max_size: str) -> float:
    error_msg = f"Invalid max size: {max_size}. Must end with 'Mb' or 'Gb'. Example: 10Mb / 2.1Gb"
    try:
        if max_size.endswith('Mb'):
            return 1e6 * float(max_size.replace('Mb', '')) 
        elif max_size.endswith('Gb'):
            return 1e9 * float(max_size.replace('Gb', ''))
        else:
            raise ValueError(error_msg)
    except ValueError:
        raise ValueError(error_msg)


def get_srrs(experiments: list[dict], max_size) -> list[dict]:
    filtered_experiments = []
    if max_size:
        max_size_bytes = parse_max_size(max_size)
        filtered_experiments = []
        for exp in short_read_experiments:
            try:
                size_in_bytes = float(exp["RUN_SET"]["@bytes"])
            except (ValueError, KeyError):
                srr = exp["EXPERIMENT"]["@accession"]
                logger.warning(f"Could not get size in bytes for {srr}")
                continue
            if size_in_bytes and size_in_bytes <= max_size_bytes:
                filtered_experiments.append(exp)
    else:
        filtered_experiments = experiments
    return [exp["EXPERIMENT"]["@accession"] for exp in filtered_experiments]
    

#####################################################
#####################################################
# MAIN
#####################################################srrs = get_srrs(long_read_experiments, args.max_size)
#####################################################

if __name__ == "__main__":
    args = parse_args()

    logger.info(f"Fetching SRA experiment UIDs for NCBI taxon ID {args.taxid}")

    #####################################################
    # SHORT READS
    #####################################################

    if args.max_short_reads_data > 0:
        logger.info(f"Fetching short read sequencing data (paired_only: {args.paired_only})")
        short_read_sra_raw_ids = get_sra_ids(args.taxid, Platform.SHORT_READ, args.paired_only)
        logger.info(f"Got {len(short_read_sra_raw_ids)} short read SRA experiment IDs")
    else:
        logger.info("Skipping short read data fetching")
        short_read_sra_raw_ids = []

    if len(short_read_sra_raw_ids) < args.max_short_reads_data:
        logger.warning(f"Could not find {args.max_short_reads_data} from Illumina platform alone. Trying with other platforms")
        fallback_short_read_sra_ids = get_sra_ids(args.taxid, Platform.SHORT_READ_FALLBACK, args.paired_only)
        short_read_sra_raw_ids += fallback_short_read_sra_ids
        logger.info(f"Got {len(fallback_short_read_sra_ids)} additional SRA experiment UIDs. Total: {short_read_sra_raw_ids}")

    #####################################################
    # LONG READS
    #####################################################

    if args.max_long_reads_data > 0:
        logger.info("Fetching long read sequencing data")
        long_read_sra_raw_ids = get_sra_ids(args.taxid, Platform.LONG_READ)
        logger.info(f"Got {len(long_read_sra_raw_ids)} long read SRA experiment IDs")
    else:
        logger.info("Skipping long read data fetching")
        long_read_sra_raw_ids = []
    
    #####################################################
    # GETTING ADDITIONAL INFORMATION FOR FETCHED SRA EXPERIMENTS 
    #####################################################

    short_read_experiments = []
    if short_read_sra_raw_ids:
        short_read_experiments = fetch_sra_experiments(short_read_sra_raw_ids, read_type='short read')
    
    long_read_experiments = []
    if long_read_sra_raw_ids:
        long_read_experiments = fetch_sra_experiments(long_read_sra_raw_ids, read_type='long read')
    
    #####################################################
    # EXPORTING DATA
    #####################################################

    if short_read_experiments:
        logger.info(f"Writing short read SRA IDs to {SHORT_READ_SRA_IDS_OUTFILE}")
        srrs = get_srrs(short_read_experiments, args.max_size)
        with open(SHORT_READ_SRA_IDS_OUTFILE, "w") as fout:
            fout.writelines([f"{srr}\n" for srr in srrs])
    
        logger.info(f"Writing metadata of short read SRAs to {SHORT_READ_METADATA_OUTFILE}")
        with open(SHORT_READ_METADATA_OUTFILE, "w") as fout:
            json.dump(short_read_experiments, fout)

    if long_read_experiments:
        logger.info(f"Writing long read SRA IDs to {LONG_READ_SRA_IDS_OUTFILE}")
        srrs = get_srrs(long_read_experiments, args.max_size)
        with open(LONG_READ_SRA_IDS_OUTFILE, "w") as fout:
            fout.writelines([f"{srr}\n" for srr in srrs])
    
        logger.info(f"Writing metadata of long read SRAs to {LONG_READ_METADATA_OUTFILE}")
        with open(LONG_READ_METADATA_OUTFILE, "w") as fout:
            json.dump(long_read_experiments, fout)

    logger.info("Done")
