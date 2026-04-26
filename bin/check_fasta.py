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
    parser.add_argument('--in', type=Path, dest='fasta_file')
    parser.add_argument('--out', type=Path, dest='outfile')
    parser.add_argument("--type", choices=["dna", "protein"])
    parser.add_argument("--fix-headers", dest="fix_headers", action="store_true")
    parser.add_argument("--fix-sequences", dest="fix_sequences", action="store_true")
    return parser.parse_args()

def check_headers_and_fix_if_needed(fasta_file: Path, fix: bool) -> list[SeqRecord]:
    records = []
    with open(fasta_file, 'r') as handle_in:
        for record in SeqIO.parse(handle_in, 'fasta'):
            if fix:
                new_id = re.sub(fr"[^{ALLOWED_CHARACTERS}]", '_', record.id)
                cleaned_record = SeqRecord(record.seq, id=new_id, name='', description='')
                records.append(cleaned_record)
            else:
                if found := re.findall(fr"[^{ALLOWED_CHARACTERS}]", record.description):
                    formated_found_characters = ", ".join([f"'{c}'" for c in found])
                    raise ValueError(
                        f"Found unauthorised characters {formated_found_characters} in header '{record.description}'.\nAllowed characters: {ALLOWED_CHARACTERS}"
                    )
                records.append(record)
    return records


def get_generic_symbol(type: str) -> str:
    return 'N' if type == 'dna' else 'X'


def check_sequences_and_fix_if_needed(records: list[SeqRecord], seq_type: str, fix_sequences: bool):
    symbol = get_generic_symbol(seq_type)
    new_records = []
    for record in records:
        sequence = str(record.seq)
        seq_set = set(sequence)
        if not seq_set.issubset(LETTERS): # has invalid chars
            invalid_chars = seq_set.difference(LETTERS)
            if fix_sequences:
                for invalid_char in invalid_chars:
                    sequence = sequence.replace(invalid_char, symbol)
                record.seq = Seq(sequence)
            else:
                raise ValueError(f"Found invalid characters {invalid_chars} in sequence.")
        new_records.append(record)
    return new_records


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    records = check_headers_and_fix_if_needed(args.fasta_file, args.fix_headers)

    records = check_sequences_and_fix_if_needed(records, args.type, args.fix_sequences)

    if args.fix_headers or args.fix_sequences:
        with open(args.outfile, 'w+') as fout:
            for record in records:
                SeqIO.write(record, fout, 'fasta')
