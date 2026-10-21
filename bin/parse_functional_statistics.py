#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
Parse functional annotation statistics from AGAT's sp_functional_statistics.pl output.
"""

import argparse
import logging
from pathlib import Path
from io import TextIOWrapper

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

AGAT_STAT_OUTFILENAME = "table_per_feature_type.txt"

ATTRIBUTE_STAT_OUTFILE_SUFFIX = "agat_func_stats.n_holding.csv"
NB_HOLDED_BY_FEATURE_OUTFILE_SUFFIX = "agat_func_stats.n_holded.csv"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", dest="result_folder", type=Path, required=True, help="Output folder from agat_sp_functional_statistics.pl"
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
    # components is like:
    # ['Nb', 'holded', 'by', 'Nb', 'gene']
    # ['Nb', 'holded', 'by', 'Nb', 'transcript', 'Nb', 'gene', 'with']
    feature = components[4]
    parsed_dict = {'feature': feature}
    # going over the second line of the header
    line = next(fin)
    # going over the second dash line
    line = next(fin)
    while True:
        line = next(fin)
        components = [s for s in line.strip().split(' ') if s not in ['', '|']] 
        attribute = components[0]
        parsed_dict[attribute] = components[1]
        # the last attribute in such tables is 'dbxref'
        if attribute == "dbxref":
            # going over the last dash line
            line = next(fin)
            break
        else:
            # going over the next dash line
            line = next(fin)
    return parsed_dict


def parse_stat_file(stat_file: Path):
    attribute_dicts = []
    nb_holded_by_feature_dicts = []
    
    with open(stat_file, "r") as fin:
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

    return attribute_dicts, nb_holded_by_feature_dicts


def write_attribute_dicts(attribute_dicts: list[dict], prefix: str, suffix: str):
    found_attributes = list(set([d['attribute'] for d in attribute_dicts]))
    for attribute in found_attributes:
        
        attr_dicts = [d for d in attribute_dicts if d['attribute'] == attribute]

        outfile = f"{prefix}.{attribute}.{suffix}"
        with open(outfile, 'w') as fout:
            for i, d in enumerate(attr_dicts):
                # writing header
                if i == 0:
                    header = ','.join([k for k in d.keys() if k != 'attribute'])
                    fout.write(header + "\n")
                line = ','.join([v for k, v in d.items() if k != 'attribute'])
                fout.write(line + "\n")


def write_nb_holded_attributes_dicts(nb_holded_by_feature_dicts: list[dict], prefix: str, suffix: str):
    found_features = list(set([d['feature'] for d in nb_holded_by_feature_dicts]))
    for feature in found_features:
        outfile = f"{prefix}.{feature}.{suffix}"
        with open(outfile, 'w') as fout:
            for i, d in enumerate(nb_holded_by_feature_dicts):
                # writing header
                if i == 0:
                    header = ','.join(list(d.keys()))
                    fout.write(header + "\n")
                line = ','.join(list(d.values()))
                fout.write(line + "\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLASS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class StatParser:

    folder: Path
    stat_files: dict
    parsed_results: dict

    def __init__(self, folder: Path):
        self.folder = folder
        self.stat_files = self.get_stat_files()
        self.parsed_results = {}


    def parse(self):
        for level1, folder_dict in self.stat_files.items():
            self.parsed_results[level1] = {}
            for level2, file in folder_dict.items():
                self.parsed_results[level1][level2] = {}
                attribute_dicts, nb_holded_by_feature_dicts = parse_stat_file(file)
                self.parsed_results[level1][level2]['attributes'] = attribute_dicts
                self.parsed_results[level1][level2]['nb_holded_by_feature'] = nb_holded_by_feature_dicts


    def get_stat_files(self) -> dict:
        # getting all subfolders like 'gene@mrna' or 'gene@transcript'
        stat_subfolders = {}
        for item in self.folder.iterdir():
            if item.is_dir() and '@' in item.name:
                level1, level2 = item.name.split('@', 2)
                if level1 not in stat_subfolders:
                    stat_subfolders[level1] = {}
                for file in item.iterdir():
                    if file.is_file() and file.name == AGAT_STAT_OUTFILENAME:
                        stat_subfolders[level1][level2] = file
        return stat_subfolders


    def export_stats(self, prefix: str):
        # looping through first level features (typically, 'gene')
        for level1 in self.parsed_results:
            # looping through second level features (typically, 'transcript' or 'mrna')
            for level2 in self.parsed_results[level1]:
                
                attribute_dicts = self.parsed_results[level1][level2]['attributes']
                write_attribute_dicts(attribute_dicts, prefix, ATTRIBUTE_STAT_OUTFILE_SUFFIX)
                
                nb_holded_by_feature_dicts = self.parsed_results[level1][level2]['nb_holded_by_feature']
                write_nb_holded_attributes_dicts(nb_holded_by_feature_dicts, prefix, NB_HOLDED_BY_FEATURE_OUTFILE_SUFFIX)      


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    args = parse_args()
    logger.info(f"Parsing stat file: {args.result_folder}")

    #####################################################
    # PARSING
    #####################################################

    stat_parser = StatParser(args.result_folder)
    stat_parser.parse()

    #####################################################
    # WRITING PARSED STATS
    #####################################################

    stat_parser.export_stats(args.prefix)

    logger.info("Done")


if __name__ == "__main__":
    main()
