#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
infer_strand_from_xs.py - infer RNA-seq library strandedness from HISAT2's XS tag.

For each spliced read 1, compares its mapped strand to the transcript strand
inferred by HISAT2 from the splice-site motif (the XS:A: tag), and reports how
often they agree (forward-stranded) or disagree (reverse-stranded).
"""
import argparse
import sys
import logging
from enum import Enum
from pathlib import Path

import pysam

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

OUTFILE = "strandedness.txt"

class Stranded(Enum):
    FORWARD = "forward"
    REVERSE = "reverse"
    UNSTRANDED = "unstranded"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", type=Path, help="BAM file aligned with HISAT2 (needs XS:A: tag)")
    parser.add_argument("--min-mapq", type=int, default=0,
                     help="minimum MAPQ to consider (default: 0)")
    parser.add_argument("--unique-only", action="store_true",
                     help="only count reads with NH:i:1 (unique alignments: skip multi mapping)")
    parser.add_argument("--single-end", action="store_true",
                     help="data is single-end (don't filter on read1)")
    return parser.parse_args()


def parse_alignment(bam_file: Path, single_end: bool, min_mapq: int, unique_only: bool) -> tuple[int, int]:
    nb_same_strand = nb_opposite_strand = 0
    nb_total = nb_skipped_unmapped = nb_skipped_secondary = nb_skipped_suppl = nb_skipped_low_qual = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            nb_total += 1
            if read.is_unmapped:
                nb_skipped_unmapped += 1
                continue
            if read.is_secondary:
                nb_skipped_secondary += 1
                continue
            if read.is_supplementary:
                nb_skipped_suppl += 1
                continue
            if not single_end and not read.is_read1:
                continue
            if read.mapping_quality < min_mapq:
                nb_skipped_low_qual += 1
                continue
            if unique_only:
                nh = read.get_tag("NH") if read.has_tag("NH") else 1
                if nh != 1:
                    continue

            # only spliced reads (CIGAR contains an N/skip operation) carry
            # a splice-motif-derived strand
            if read.cigartuples is None or not any(op == 3 for op, _ in read.cigartuples):
                continue

            if not read.has_tag("XS"):
                continue
            xs = read.get_tag("XS")  # '+' or '-'

            r1_strand = "-" if read.is_reverse else "+"

            if r1_strand == xs:
                nb_same_strand += 1
            else:
                nb_opposite_strand += 1

    logger.info(f"# total: {nb_total}")
    logger.info(f"# unmapped (skipped): {nb_skipped_unmapped}")
    logger.info(f"# secondary (skipped): {nb_skipped_secondary}")
    logger.info(f"# supplementary (skipped): {nb_skipped_suppl}")
    logger.info(f"# low quality (skipped): {nb_skipped_low_qual}")
                
    return nb_same_strand, nb_opposite_strand


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

def main():

    args = parse_args()
    
    nb_same_strand, nb_opposite_strand = parse_alignment(args.bam, args.single_end, args.min_mapq, args.unique_only)

    total = nb_same_strand + nb_opposite_strand
    logger.info(f"# total kept for strandedness inference: {total}")
    
    if total == 0:
        raise ValueError("No spliced reads with an XS tag found.")

    logger.info(f"# of read1 same strand as transcript:     {nb_same_strand} ({100 * nb_same_strand / total:.1f}%)")
    logger.info(f"# of read1 opposite strand to transcript: {nb_opposite_strand} ({100 * nb_opposite_strand / total:.1f}%)")

    if nb_same_strand / total >= 0.9:
        verdict = "forward-stranded (ISF/SF) -> hisat2 --rna-strandness FR (paired) / F (single-end)"
        strandedness = Stranded.FORWARD
    elif nb_opposite_strand / total >= 0.9:
        verdict = "reverse-stranded (ISR/SR) -> hisat2 --rna-strandness RF (paired) / R (single-end)"
        strandedness = Stranded.REVERSE
    else:
        verdict = "unstranded (IU/U) -> omit --rna-strandness"
        strandedness = Stranded.UNSTRANDED
        
    logger.info(f"Verdict: {verdict}")

    with open(OUTFILE, 'w') as fout:
        fout.write(str(strandedness))
        

if __name__ == "__main__":
    main()