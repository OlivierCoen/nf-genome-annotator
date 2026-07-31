#!/usr/bin/env python

# Adapted from from nf-core/fetchngs. Released under the MIT license.


import argparse
import csv
import logging
import re
import sys

import requests
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


# Example ids supported by this script
SRA_IDS = (
    "PRJNA63463",
    "SAMN00765663",
    "SRA023522",
    "SRP003255",
    "SRR390278",
    "SRS282569",
    "SRX111814",
)
ENA_IDS = (
    "PRJEB7743",
    "SAMEA3121481",
    "ERA2421642",
    "ERP120836",
    "ERR674736",
    "ERS4399631",
    "ERX629702",
)
DDBJ_IDS = (
    "PRJDB4176",
    "SAMD00114846",
    "DRA008156",
    "DRP004793",
    "DRR171822",
    "DRS090921",
    "DRX162434",
)
GEO_IDS = ("GSE18729", "GSM465244")
ID_REGEX = re.compile(r"^([A-Z]+)([0-9]+)$")
PREFIX_LIST = sorted(
    {ID_REGEX.match(id).group(1) for id in SRA_IDS + ENA_IDS + DDBJ_IDS + GEO_IDS}
)


# List of metadata fields fetched from the ENA API - can be overriden by options
# `-ef` or `--ena_metadata_fields`.
# Full list of accepted fields can be obtained here:
# https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&format=tsv&result=read_run
ENA_METADATA_FIELDS = (
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "secondary_study_accession",
    "submission_accession",
    "run_alias",
    "experiment_alias",
    "sample_alias",
    "study_alias",
    "library_layout",
    "library_selection",
    "library_source",
    "library_strategy",
    "library_name",
    "instrument_model",
    "instrument_platform",
    "base_count",
    "read_count",
    "tax_id",
    "scientific_name",
    "sample_title",
    "experiment_title",
    "study_title",
    "sample_description",
    "fastq_md5",
    "fastq_bytes",
    "fastq_ftp",
    "fastq_galaxy",
    "fastq_aspera",
)

OUTFILE = "ena_ftp_links.txt"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ENA_API_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

#####################################################
#####################################################
# CLASSES
#####################################################
#####################################################


class DatabaseIdentifierChecker:
    """Define a service class for validating database identifiers."""

    _VALID_PREFIXES = frozenset(PREFIX_LIST)

    @classmethod
    def is_valid(cls, identifier):
        """
        Check the validity of the given database identifier.

        Args:
            identifier (str): A short identifier presumably belonging to one of the
                supported databases.

        Returns:
            bool: Whether or not the identifier is valid.

        """
        match = ID_REGEX.match(identifier)
        if match is None:
            return False
        return match.group(1) in cls._VALID_PREFIXES


class DatabaseResolver:
    """Define a service class for resolving various identifiers to experiments."""

    _GEO_GSM_PREFIXES = {"GSM"}
    _GEO_GSE_PREFIXES = {"GDS", "GSE"}
    _SRA_PREFIXES = {
        "DRA",
        "DRP",
        "DRS",
        "DRX",
        "PRJDB",
        "SAMD",
    }
    _ENA_PREFIXES = {"ERR", "SRR", "SAMN", "DRR"}

    @classmethod
    def expand_identifier(cls, identifier):
        """
        Expand the given identifier to potentially multiple experiment identifiers.

        Args:
            identifier (str): A short identifier presumably belonging to one of the
                supported databases.

        Returns:
            list: A list of one or more SRA/ENA experiment identifiers.

        """
        prefix = ID_REGEX.match(identifier).group(1)
        if prefix in cls._GEO_GSM_PREFIXES:
            return cls._gsm_to_srx(identifier)
        elif prefix in cls._GEO_GSE_PREFIXES:
            return cls._gse_to_srx(identifier)
        elif prefix in cls._SRA_PREFIXES:
            return cls._id_to_srx(identifier)
        elif prefix in cls._ENA_PREFIXES:
            return cls._id_to_erx(identifier)
        else:
            return [identifier]

    @classmethod
    def _id_to_srx(cls, identifier):
        """Resolve the identifier to SRA experiments."""
        params = {
            "id": identifier,
            "db": "sra",
            "rettype": "runinfo",
            "retmode": "text",
        }
        response = fetch_url_content(EFETCH_URL, params=params)
        return [row["Experiment"] for row in open_table(response, delimiter=",")]

    @classmethod
    def _gsm_to_srx(cls, identifier):
        """Resolve the GEO identifier to SRA experiments."""
        ids = []
        params = {"term": identifier, "db": "sra", "retmode": "json"}
        response = fetch_url_content(ESEARCH_URL, params=params)
        data = response.json()
        gsm_ids = data["esearchresult"]["idlist"]
        for gsm_id in gsm_ids:
            ids += cls._id_to_srx(gsm_id)
        return ids

    @classmethod
    def _gds_to_gsm(cls, identifier):
        """Resolve the GEO UIDs to GSM IDs to then resolve to SRA IDs."""
        ids = []
        params = {"id": identifier, "db": "gds", "retmode": "json", "retmax": 10}
        response = fetch_url_content(
            ESUMMARY_URL,
            params=params,
        )
        data = response.json()
        for each in data["result"][identifier]["samples"][0:]:
            ids += cls._gsm_to_srx(each["accession"])
        return ids

    @classmethod
    def _gse_to_srx(cls, identifier):
        """Resolve the GSE identifier to GEO UIDs."""
        ids = []
        params = {"term": identifier, "db": "gds", "retmode": "json"}
        response = fetch_url_content(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params=params
        )
        data = response.json()
        gds_uids = data["esearchresult"]["idlist"]
        for gds_uid in gds_uids:
            ids += cls._gds_to_gsm(gds_uid)
        return ids

    @classmethod
    def _id_to_erx(cls, identifier):
        """Resolve the identifier to ENA experiments."""
        fields = ["run_accession", "experiment_accession"]
        params = {
            "accession": identifier,
            "result": "read_run",
            "fields": ",".join(fields),
        }
        response = fetch_url_content(ENA_API_FILEREPORT_URL, params=params)
        return [
            row["experiment_accession"]
            for row in open_table(response.text, delimiter="\t")
        ]


class ENAMetadataFetcher:
    """Define a service class for fetching metadata from ENA."""

    def __init__(self, ena_metadata_fields: list[str], **kwargs):
        """
        Initialize the service with the desired metadata fields.

        Args:
            ena_metadata_fields (list[str]): An iterable of the desired fields.
            **kwargs: Passed to parent constructor.
        """
        super().__init__(**kwargs)
        self._params = {
            "result": "read_run",
            "fields": ",".join(ena_metadata_fields),
        }

    def open_experiment_table(self, accession: str) -> list[dict]:
        """
        Open the metadata table belonging to the given experiment accession.

        Args:
            accession (str): An ENA experiment accession.

        Returns:
            csv.DictReader: A CSV reader instance of the metadata.

        """
        params = {**self._params, "accession": accession}
        response = fetch_url_content(ENA_API_FILEREPORT_URL, params=params)
        return open_table(response.text, delimiter="\t")


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def open_table(content: str, delimiter=",") -> list[dict]:
    """
    Return a list derived from a CSV reader instance from the given response.

    Args:
        response (Response): An instance of the local HTTP response class.
        delimiter (str): The delimiter separating the table fields.

    Returns:
            csv.DictReader: A CSV reader instance of the response body.

    """
    reader = csv.DictReader(content.splitlines(), delimiter=delimiter)
    return [row for row in reader]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download and create a run information metadata file from SRA / "
        "ENA / DDBJ / GEO identifiers.",
        epilog="Example usage: python fetch_sra_runinfo.py <FILE_IN> <FILE_OUT>",
    )
    parser.add_argument(
        "--accession", required=True, type=str, help="Database identifier."
    )
    return parser.parse_args()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def fetch_url_content(url: str, params: dict):
    HEADERS = {"accept": "application/json", "content-type": "application/json"}
    response = requests.get(url, params=params, headers=HEADERS)
    response.raise_for_status()
    return response


def check_accession(accession: str):
    if not DatabaseIdentifierChecker.is_valid(accession):
        id_str = ", ".join([x + "*" for x in PREFIX_LIST])
        logger.error(
            f"Invalid accession: {accession}\n"
            f"Please provide a valid database id starting with {id_str}!\n"
        )
        sys.exit(1)
    ids = DatabaseResolver.expand_identifier(accession)
    if not ids:
        logger.error(f"No matches found for database id {accession}")
        sys.exit(1)


def main():
    args = parse_args()

    logger.info("Checking accessions")
    check_accession(args.accession)

    logger.info("Fetching ENA FTP URLs")
    ena_metadata_fields = ["run_accession", "fastq_ftp"]
    ena_fetcher = ENAMetadataFetcher(ena_metadata_fields)

    ena_table = ena_fetcher.open_experiment_table(args.accession)
    logger.info(f"Found {len(ena_table)} ENA entries")
    logger.info(ena_table)

    logger.info("Writing FTP URLs to file")
    ftp_urls = [url for row in ena_table for url in row["fastq_ftp"].split(";")]
    with open(OUTFILE, "w") as fout:
        fout.writelines([f"ftp://{url}\n" for url in ftp_urls])

    logger.info("Done")


if __name__ == "__main__":
    main()
