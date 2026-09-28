#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""detect_interleaved.py - check whether a FASTQ file is interleaved (R1,R2,R1,R2...).

Usage: detect_interleaved.py reads.fastq[.gz] [n_pairs=1000]
"""
import gzip
import re
import sys
from itertools import islice
import logging
import argparse
from pathlib import Path


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

DEFAULT_NPAIRS = 1000
OUTFILE = "status.txt"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check whether a FASTQ file is interleaved (R1,R2,R1,R2...)"
    )
    parser.add_argument("--file", type=Path, required=True)
    parser.add_argument("--npairs", type=int, default=DEFAULT_NPAIRS, help="Number of pairs to check")
    return parser.parse_args()
    

def open_any(path):
    with open(path, "rb") as fin:
        gz = fin.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def parse_header(line):
    """Return (base_name, mate) where mate is 1, 2 or 0 (unknown)."""
    parts = line.rstrip("\n").lstrip("@").split(None, 1)
    name = parts[0]
    comment = parts[1] if len(parts) > 1 else ""

    mate = 0
    m = re.search(r"/([12])$", name)
    if m:                                   # name/1, name/2
        mate = int(m.group(1))
        name = name[:-2]
    else:
        m = re.match(r"([12])(:|\s|$)", comment)
        if m:                               # CASAVA "1:N:0:..." or SRA "1 length=.."
            mate = int(m.group(1))
    return name, mate


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

def main():
    args = parse_args()

    with open_any(args.file) as fh:
        headers = [l for l in islice(fh, 0, args.npairs * 2 * 4, 4)]

    pairs = len(headers) // 2
    if pairs == 0:
        raise RuntimeError("Not enough records to decide")

    same = marked = proper = 0
    for i in range(pairs):
        n1, m1 = parse_header(headers[2 * i])
        n2, m2 = parse_header(headers[2 * i + 1])
        if n1 == n2:
            same += 1
        if m1:
            marked += 1
        if n1 == n2 and m1 == 1 and m2 == 2:
            proper += 1

    logger.info(f"Pairs checked: {pairs} | same name: {same} | proper 1/2 mates: {proper}")

    interleaved = None
    if same == pairs and proper == pairs:
        logger.info("INTERLEAVED")
        interleaved = "true"
    elif same == pairs and marked == 0:
        logger.info("LIKELY INTERLEAVED (identical consecutive names, no mate markers)")
        interleaved = "true"
    elif same == 0:
        logger.info("NOT interleaved (consecutive records have different names)")
        interleaved = "false"
    else:
        raise ValueError("MIXED/BROKEN: some consecutive records pair up, others do not")

    with open(OUTFILE, 'w') as fout:
        fout.write(interleaved)

    logger.info(f"Interleaved: {interleaved}")
        

if __name__ == "__main__":
    main()