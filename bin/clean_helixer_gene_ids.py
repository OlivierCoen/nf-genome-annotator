#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, help="Annotation file in GFF3")
    parser.add_argument("--out", required=True, help="Output GFF3")
    return parser.parse_args()


def main():
    args = parse_args()

    logger.info(f"Cleaning gene IDs in {args.gff}")

    with open(args.gff, 'r') as fin, open(args.out, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                cleaned_line = line
            else:
                fields = line.split('\t')
                attributes = [re.sub(r'=_X', '=X', attr) for attr in fields[8].split(";")]
                cleaned_fields = fields[:8] + [";".join(attributes)]
                cleaned_line = '\t'.join(cleaned_fields)
            # writing to output file
            fout.write(cleaned_line)

    logger.info("Done")


if __name__ == "__main__":
    main()
            