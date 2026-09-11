#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

from pathlib import Path
import argparse
import random
import logging

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', type=Path, dest='input_file', required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--nb', dest='nb_to_sample', type=int, required=True)
    parser.add_argument('--out', type=Path, dest='outfile', required=True)
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    random.seed(args.seed)

    logger.info(f"Reading input file {args.input_file}")
    with open(args.input_file, 'r') as fin:
        items = fin.readlines()

    sampled_items = random.sample(items, k=args.nb_to_sample)

    logger.info(f"Writing output file {args.outfile}")
    with open(args.outfile, 'w+') as fout:
        fout.writelines(sampled_items)
        