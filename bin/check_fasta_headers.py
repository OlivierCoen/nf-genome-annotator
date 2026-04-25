#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

from pathlib import Path
import argparse
import logging
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

ALLOWED_CHARACTERS = "a-zA-Z0-9.-_,;"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', type=Path, dest='fasta_file')
    parser.add_argument('--out', type=Path, dest='outfile')
    parser.add_argument("--fix", action="store_true")
    return parser.parse_args()

def get_records_with_cleaned_headers(fasta_file: Path) -> list[SeqRecord]:
    cleaned_records = []
    with open(fasta_file, 'r') as handle_in:
        for record in SeqIO.parse(handle_in, 'fasta'):
            new_id = re.sub(fr"[^{ALLOWED_CHARACTERS}]", '_', record.id)
            cleaned_record = SeqRecord(record.seq, id=new_id, name='', description='')
            cleaned_records.append(cleaned_record)
    return cleaned_records


def check_headers(fasta_file: Path) -> None:
    with open(fasta_file, 'r') as handle_in:
        for record in SeqIO.parse(handle_in, 'fasta'):
            if found := re.findall(fr"[^{ALLOWED_CHARACTERS}]", record.description):
                formated_found_characters = ", ".join([f"'{c}'" for c in found])
                raise ValueError(
                    f"Found unauthorised characters {formated_found_characters} in header '{record.description}'.\nAllowed characters: {ALLOWED_CHARACTERS}"
                )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    if args.fix:
        cleaned_records = get_records_with_cleaned_headers(args.fasta_file)

        with open(args.outfile, 'w+') as fout:
            for record in cleaned_records:
                SeqIO.write(record, fout, 'fasta')

    else:
        check_headers(args.fasta_file)
