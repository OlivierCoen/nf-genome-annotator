#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import yaml

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

MRNA_OUTFILE_SUFFIX = "mrna_gff_stats.csv"
RNA_OUTFILE_SUFFIX = "rna_gff_stats.csv"
TRANSCRIPT_OUTFILE_SUFFIX = "transcript_gff_stats.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Clean GTF gene IDs")
    parser.add_argument(
        "--gff", dest="gff_file", type=Path, required=True, help="Input GTF file"
    )
    parser.add_argument("--prefix", type=str, required=True, help="Outfile name prefix")
    return parser.parse_args()


def main():
    args = parse_args()
    logger.info(f"Parsing GTF file: {args.gff_file}")

    with open(args.gff_file, "r") as fin:
        gff_data = yaml.safe_load(fin)

        for feature, suffix in zip(
            ["mrna", "rna", "transcript"],
            [MRNA_OUTFILE_SUFFIX, RNA_OUTFILE_SUFFIX, TRANSCRIPT_OUTFILE_SUFFIX],
        ):
            if feature in gff_data:
                for isoform_type in ["without_isoforms", "with_isoforms"]:
                    if isoform_type in gff_data[feature]:
                        data = gff_data[feature][isoform_type].get("value")
                        if data is not None:
                            # prepending prefix
                            data = {"file": args.prefix} | data
                            # getting header and data
                            header = (
                                ",".join([f'"{key}"' for key in list(data.keys())])
                                + "\n"
                            )
                            data = ",".join(
                                [str(value) for value in list(data.values())]
                            )

                            outfile = f"{args.prefix}_{isoform_type}_{suffix}"
                            with open(outfile, "a") as fout:
                                fout.write(header)
                                fout.write(data)

    logger.info("Done")


if __name__ == "__main__":
    main()
