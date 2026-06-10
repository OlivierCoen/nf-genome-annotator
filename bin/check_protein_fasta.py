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
    parser.add_argument("--minlen", type=int, dest="min_sequence_length", required=True)

    return parser.parse_args()

def parse_fasta_and_fix_headers_if_needed(fasta_file: Path) -> list[SeqRecord]:
    records = []
    with open(fasta_file, 'r') as handle_in:
        for record in SeqIO.parse(handle_in, 'fasta'):
            new_id = re.sub(fr"[^{ALLOWED_CHARACTERS}]", '_', record.id)
            cleaned_record = SeqRecord(record.seq, id=new_id, name='', description='')
            records.append(cleaned_record)
    return records


def check_sequences_and_fix_if_needed(records: list[SeqRecord]):
    new_records = []
    for record in records:
        sequence = str(record.seq)
        seq_set = set(sequence)
        if not seq_set.issubset(LETTERS): # has invalid chars
            invalid_chars = seq_set.difference(LETTERS)
            for invalid_char in invalid_chars:
                sequence = sequence.replace(invalid_char, "X")
            record.seq = Seq(sequence)
        new_records.append(record)
    return new_records


def filter_sequences_by_length(records: list[SeqRecord], min_length: int) -> list[SeqRecord]:
    return [record for record in records if len(record.seq) >= min_length]


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    records = parse_fasta_and_fix_headers_if_needed(args.fasta_file)

    records = check_sequences_and_fix_if_needed(records)

    if args.min_sequence_length >= 1:
        records = filter_sequences_by_length(records, args.min_sequence_length)

    with open(args.outfile, 'w+') as fout:
        for record in records:
            SeqIO.write(record, fout, 'fasta')
