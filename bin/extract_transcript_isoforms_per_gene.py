#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from pathlib import Path

import logging
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GFF_COLUMNS = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gff", 
        dest="annot_file",
        required=True, 
        help="Annotation file in GFF3"
    )
    parser.add_argument(
        "--out", 
        dest="outfile",
        required=True, 
        help="Output GFF3"
    )
    return parser.parse_args()


def parse_gff3(file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file, 
        separator="\t", 
        comment_prefix="#", 
        has_header=False, 
        new_columns=GFF_COLUMNS
    )



def main():
    
    args = parse_args()

    annot_lf = parse_gff3(args.annot_file)

    (
    annot_lf
        .filter(pl.col("feature") == "transcript")
        .with_columns(
            pl.col("attribute").str.extract(r"ID=(.+?)(?:;|$)", 1).alias("transcript_id"),
            pl.col("attribute").str.extract(r"Parent=(.+?)(?:;|$)", 1).alias("gene_id")
        )
        .group_by("gene_id")
        .agg(pl.col("transcript_id").alias("isoforms"))
        .select(
            pl.col("isoforms").list.join(separator=";")
        )
        .sink_csv(args.outfile, include_header=False)
    )

    logger.info("Done")


if __name__ == "__main__":
    main()