#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from typing import List

from Bio import SeqIO


def find_identical_sequences(fasta_file: str) -> List[List[str]]:
    seq_dict = defaultdict(set)

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq).upper()
        seq_dict[seq_str].add(record.id)

    identical_groups = [ids for ids in seq_dict.values() if len(ids) > 1]
    identical_groups = [sorted(list(group)) for group in identical_groups]

    return identical_groups


def write_identical_sequences(identical_groups: List[List[str]], basename: str):
    output_file = basename + "_ident_seq.csv"
    with open(output_file, 'w') as f:
        for group in identical_groups:
            f.write(",".join(group) + "\n")
    print(f"Identical sequences groups written to {output_file}")


def write_unique_fasta(identical_groups: List[List[str]], input_fasta: str, basename: str):
    unique_file = basename + "_unique.fasta"
    written_ids = set()

    # Add the first identifier from each identical group to the set
    for group in identical_groups:
        if group:
            written_ids.add(group[0])

    with open(unique_file, 'w') as f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in written_ids or all(record.id not in group for group in identical_groups):
                SeqIO.write(record, f, "fasta")
                written_ids.add(record.id)

    print(f"Unique sequences fasta written to {unique_file}")


def main():
    parser = argparse.ArgumentParser(description="Znajdź identyczne sekwencje w pliku FASTA.")
    parser.add_argument("-i", "--input", required=True, help="Plik wejściowy w formacie FASTA")
    parser.add_argument("-o", "--output_dir", default="results", help="Katalog wyjściowy do zapisania wyników")
    args = parser.parse_args()

    basename = ''.join(os.path.splitext(args.input.split("/")[-1])[:-1])
    basename = os.path.join(args.output_dir, basename)

    identical = find_identical_sequences(args.input)
    write_identical_sequences(identical, basename)
    write_unique_fasta(identical, args.input, basename)


if __name__ == "__main__":
    main()
