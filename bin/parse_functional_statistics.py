#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path
from io import TextIOWrapper

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

NB_HOLDED_BY_FEATURE_OUTFILE_SUFFIX = "nb_holded_by_features.csv"
ATTRIBUTE_STAT_OUTFILE_SUFFIX = "stats.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

def parse_args():
    parser = argparse.ArgumentParser(description="Parse functional annotation statistics")
    parser.add_argument(
        "--stat", dest="stat_file", type=Path, required=True, help="Statistics file"
    )
    parser.add_argument("--prefix", type=str, required=True, help="Outfile name prefix")
    return parser.parse_args()


def parse_line_with_numbers(line: str, key: str):
    """
    Parses a line like:
    Nb five_prime_utr without <name> attribute = 3
    """
    components = line.strip().split(" ")
    return {
        'feature': components[1], 
        'attribute': components[3].lstrip('<').rstrip('>'), 
        key: components[6]
    }

def parse_table(fin: TextIOWrapper) -> dict:
    # going to the first real line of the table
    line = next(fin)
    # parsing the feature mentionned in first line of the table header
    # such a line is like like:
    # |                         |      Nb holded by       |    Nb three_prime_utr   |
    components = [s for s in line.strip().split(' ') if s not in ['', '|']] 
    feature = components[-1]
    parsed_dict = {'feature': feature}
    # going over the second line of the header
    line = next(fin)
    # going over the second dash line
    line = next(fin)
    while True:
        line = next(fin)
        components = [s for s in line.strip().split(' ') if s not in ['', '|']] 
        attribute = components[0]
        parsed_dict[f"nb_{attribute}_holded"] = components[1]
        # the last attribute in such tables is 'dbxref'
        if attribute == "dbxref":
            # going over the last dash line
            line = next(fin)
            break
        else:
            # going over the next dash line
            line = next(fin)
    return parsed_dict


def main():
    args = parse_args()
    logger.info(f"Parsing stat file: {args.stat_file}")

    #####################################################
    # PARSING
    #####################################################

    with open(args.stat_file, "r") as fin:

        attribute_dicts = []
        nb_holded_by_feature_dicts = []
        
        while True:
            
            try:
                line = next(fin)
            except StopIteration:
                break

            if line.startswith(("\n", "Functional")):
                # head of the file
                continue

            line = line.strip()
            
            if line.startswith("Nb"):
                # lines like:
                # Nb cds = 1433
                # Nb cds with <name> attribute = 0
                # Nb cds without <name> attribute = 1433
                nb_components = len(line.split(" "))
                if nb_components == 4:
                    # we do not want to parse the number of features
                    # since it is performed in agat spstatistics
                    continue
                elif nb_components == 7:
                    # first parsing the 'with' line
                    parsed_with_dict = parse_line_with_numbers(line, key="with")
                    line = next(fin)
                    parsed_without_dict = parse_line_with_numbers(line, key="without")
                    # marging both dicts
                    parsed_dict = parsed_with_dict | parsed_without_dict
                    attribute_dicts.append(parsed_dict)
                else:
                    logger.warning(f"Unrecognised line: {line}")
            elif line.startswith("_"):
                # head of new table reached
                parsed_dict = parse_table(fin)
                nb_holded_by_feature_dicts.append(parsed_dict)
                
            else:
                logger.warning(f"Unrecognised line: {line}")

    #####################################################
    # WRITING PARSED STATS
    #####################################################
    
    found_attributes = list(set([d['attribute'] for d in attribute_dicts]))
    for attribute in found_attributes:
        
        attr_dicts = [d for d in attribute_dicts if d['attribute'] == attribute]
        attribute_stat_outfile = f"{args.prefix}.{attribute}.{ATTRIBUTE_STAT_OUTFILE_SUFFIX}"
        
        with open(attribute_stat_outfile, 'w') as fout:
            for i, d in enumerate(attr_dicts):
                # writing header
                if i == 0:
                    header = ','.join([k for k in d.keys() if k != 'attribute'])
                    fout.write(header + "\n")
                line = ','.join([v for k, v in d.items() if k != 'attribute'])
                fout.write(line + "\n")
                
    nb_holded_by_feature_outfile = f"{args.prefix}.{NB_HOLDED_BY_FEATURE_OUTFILE_SUFFIX}"
    with open(nb_holded_by_feature_outfile, 'w') as fout:
        for i, d in enumerate(nb_holded_by_feature_dicts):
            # writing header
            if i == 0:
                header = ','.join(list(d.keys()))
                fout.write(header + "\n")
            line = ','.join(list(d.values()))
            fout.write(line + "\n")

    logger.info("Done")


if __name__ == "__main__":
    main()
