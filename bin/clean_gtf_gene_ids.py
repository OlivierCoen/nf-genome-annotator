#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import csv
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

GTF_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

ID_WITH_QUOTES_PATTERN = r'("[a-zA-Z0-9._+-|,]+")'
ID_PATTERN = r"([a-zA-Z0-9._+-|,]+)"
TRANSCRIPT_ID_PATTERN = r"(g[0-9]+\.t[0-9]+)"
GENE_ID_PATTERN = r"(g[0-9]+)"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Clean GTF gene IDs")
    parser.add_argument(
        "--gtf", dest="gtf_file", type=Path, required=True, help="Input GTF file"
    )
    parser.add_argument(
        "--out", dest="outfile", type=Path, required=True, help="Outfile name"
    )
    return parser.parse_args()


def clean_gtf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Cleans dataframe corresponding to gtf file
    After BRAKER / TSEBRA pipeline, lots of entries have things like 'g12' or 'g236.t6' as attributes
    The purpose of this function is to retrieve gene ids and transcript ids for all entries
    So that all entries have a corresponding gene_id and transcript_id
    Write new cleaned gtf dataframe to outfile
    """
    # parsing gene_id from attribute field
    df["gene_id"] = make_intermediate_column(df, "gene_id")

    # getting entries where attributes is not correctly shaped
    nan_df = df[df["gene_id"].isna()]

    # guessing transcript_id from gene_id for those entries
    nan_df["transcript_id"] = nan_df["attribute"].str.extract(TRANSCRIPT_ID_PATTERN)

    # separating gene_id from transcript_id for these entries
    nan_df["gene_id_from_transcript"] = nan_df["transcript_id"].str.extract(
        GENE_ID_PATTERN
    )
    nan_df["gene_id"] = nan_df["gene_id_from_transcript"]
    nan_df["gene_id"] = np.where(
        nan_df["gene_id"].isna(), nan_df["attribute"], nan_df["gene_id"]
    )

    # making new attribute field
    nan_df["attribute"] = nan_df.apply(
        lambda row: make_new_attribute(row["gene_id"], row["transcript_id"]), axis=1
    )
    nan_df.drop(columns=["gene_id_from_transcript", "transcript_id"], inplace=True)

    # putting back the good entries and these fixed entries in a new dataframe
    new_df = pd.concat([df[~df["gene_id"].isna()], nan_df])
    new_df = new_df.sort_values(
        by=["seqname", "gene_id", "start", "end", "strand"], ascending=True
    )
    new_df = new_df.drop(columns=["gene_id"])

    return new_df


def make_intermediate_column(df: pd.DataFrame, attribute: str) -> str:
    whole_attribute_pattern = rf'({attribute} "[a-zA-Z0-9._+-|,]+")'
    # we ensure the transcript_id will be unique
    intermediate_col = df["attribute"].str.extract(whole_attribute_pattern)[0]
    intermediate_col = intermediate_col.str.extract(ID_WITH_QUOTES_PATTERN)[0]
    attribute_col = intermediate_col.str.extract(ID_PATTERN)[0]
    return attribute_col


def make_new_attribute(gene_id: str, transcript_id: str) -> str:
    new_attribute = ""
    if not pd.isnull(transcript_id):
        new_attribute += f'transcript_id "{transcript_id}";'
    new_attribute += f' gene_id "{gene_id}";'
    return new_attribute


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    df = pd.read_csv(str(args.gtf_file), sep="\t", names=GTF_COLUMNS, comment="#")

    df = clean_gtf(df)

    df.to_csv(
        str(args.outfile), sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )
