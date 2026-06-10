#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

from pathlib import Path
import argparse
import logging
import re
import string
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord, Seq

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

ALLOWED_CHARACTERS = "a-zA-Z0-9.-_,;"
LETTERS = set(string.ascii_lowercase + string.ascii_uppercase)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', type=Path, dest='fasta_file', required=True)
    parser.add_argument('--out', type=Path, dest='outfile', required=True)
    return parser.parse_args()

def parse_fasta_and_check_headers(fasta_file: Path) -> list[SeqRecord]:
    records = []
    with open(fasta_file, 'r') as handle_in:
        for record in SeqIO.parse(handle_in, 'fasta'):
            if found := re.findall(fr"[^{ALLOWED_CHARACTERS}]", record.id):
                formated_found_characters = ", ".join([f"'{c}'" for c in set(found)])
                raise ValueError(
                    f"Found unauthorised characters {formated_found_characters} in header '{record.id}'.\nAllowed characters: {ALLOWED_CHARACTERS}"
                )
            # removing description
            record.description = ""
            records.append(record)
    return records


def fix_sequences(records: list[SeqRecord]) -> list[SeqRecord]:
    new_records = []
    for record in records:
        sequence = str(record.seq)
        seq_set = set(sequence)
        if not seq_set.issubset(LETTERS): # has invalid chars
            invalid_chars = seq_set.difference(LETTERS)
            for invalid_char in invalid_chars:
                sequence = sequence.replace(invalid_char, "N")
            record.seq = Seq(sequence)
        new_records.append(record)
    return new_records


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    records = parse_fasta_and_check_headers(args.fasta_file)

    records = fix_sequences(records)

    with open(args.outfile, 'w+') as fout:
        for record in records:
            SeqIO.write(record, fout, 'fasta')
