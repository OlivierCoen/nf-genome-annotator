#!/usr/bin/env python3

# Adapted from https://github.com/MultiQC/MultiQC/blob/main/multiqc/modules/rseqc/infer_experiment.py

"""
parse_rseqc_infer_experiment_report.py - Parse output from RSEQC infer_experiment

Parse text report and outputs a file containing the strandedness of RSEQC input file
"""

import logging
import re
import argparse
from enum import Enum
from pathlib import Path

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

OUTFILE = "strandedness.txt"

THRESHOLD = 0.9

class Stranded(Enum):
    FORWARD = "forward"
    REVERSE = "reverse"
    UNSTRANDED = "unstranded"
    UNKNOWN = "unknown"


class Report(Enum):
    PAIRED_END_SENSE = "paired_end_sense"
    PAIRED_END_ANTISENSE = "paired_end_antisense"
    SINGLE_END_SENSE = "single_end_sense"
    SINGLE_END_ANTISENSE = "single_end_antisense"
    FAILED = "failed"
    UNKNOWN = "unknown"


PATTERNS = {
    Report.PAIRED_END_SENSE: r"\"1\+\+,1--,2\+-,2-\+\": (\d\.\d+)",
    Report.PAIRED_END_ANTISENSE: r"\"1\+-,1-\+,2\+\+,2--\": (\d\.\d+)",
    Report.SINGLE_END_SENSE: r"\"\+\+,--\": (\d\.\d+)",
    Report.SINGLE_END_ANTISENSE: r"\+-,-\+\": (\d\.\d+)",
    Report.FAILED: r"Fraction of reads failed to determine: (\d\.\d+)",
    Report.UNKNOWN: r"Unknown Data type"
}

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", 
        type=Path, 
        required=True,
        help="Species name (or Taxon ID)"
    )
    return parser.parse_args()


def parse_report(file: Path) -> Stranded:
    results = {}
    
    with open(file, 'r') as fin:
        for line in fin:
            
            if line.strip() == PATTERNS[Report.UNKNOWN]:
                return Stranded.UNKNOWN
                
            for k, r in PATTERNS.items():
                match = re.search(r, line)
                if match:
                    results[k] = float(match.group(1))
    
    return infer_strandedness_from_results(results)


def get_strandedness(sense_pct: float | None, antisense_pct: float | None) -> Stranded:
    if sense_pct is not None and sense_pct >= THRESHOLD:
        return Stranded.FORWARD
    elif antisense_pct is not None and antisense_pct >= THRESHOLD:
        return Stranded.REVERSE
    else:
        return Stranded.UNSTRANDED


def is_paired_end(results: dict) -> bool:
    return Report.PAIRED_END_SENSE in results or Report.PAIRED_END_ANTISENSE in results


def is_single_end(results: dict) -> bool:
    return Report.SINGLE_END_SENSE in results or Report.SINGLE_END_ANTISENSE in results
    

def infer_strandedness_from_results(results: dict) -> Stranded:
    if Report.UNKNOWN in results:
        return Stranded.UNKNOWN
    # Paired end
    if is_paired_end(results): # assuming that no single end part was found
        return get_strandedness(
            results.get(Report.PAIRED_END_SENSE), 
            results.get(Report.PAIRED_END_ANTISENSE)
        )
    # Single end
    if is_single_end(results):
        return get_strandedness(
            results.get(Report.SINGLE_END_SENSE), 
            results.get(Report.SINGLE_END_ANTISENSE)
        )
    # should not happen
    return Stranded.UNKNOWN
    

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

def main():

    args = parse_args()

    strandedness = parse_report(args.report)
    
    with open(OUTFILE, 'w') as fout:
        fout.write(strandedness.value)
        

if __name__ == "__main__":
    main()