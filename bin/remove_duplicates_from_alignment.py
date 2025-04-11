#!/usr/bin/env python3

import argparse
from typing import List

import Bio
from Bio.SeqIO import SeqRecord


def remove_duplicates(fasta_file: str, duplicated_ids_file: str) -> List[SeqRecord]:
    """
    Remove duplicate sequences from a FASTA file based on a CSV file with duplicated IDs.
    """
    sequences_to_remove = set()
    with open(duplicated_ids_file, 'r') as f:
        for line in f:
            # Extract the duplicated IDs from the line
            ids = line.strip().split(",")
            # Add all but the first ID to the set of reads to remove
            sequences_to_remove.update(ids[1:])

    # Read the FASTA file and filter out the duplicated sequences
    records = []
    for record in Bio.SeqIO.parse(fasta_file, "fasta"):
        if record.id not in sequences_to_remove:
            records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a FASTA file.")
    parser.add_argument("-a", "--alignment", required=True, help="Input FASTA file")
    parser.add_argument("-i", "--duplicated_ids", help="csv file with duplicated ids")
    parser.add_argument("-o", "--output", help="Output alignment file")
    args = parser.parse_args()

    records_to_save = remove_duplicates(args.alignment, args.duplicated_ids)
    # Write the unique sequences to the output FASTA file
    Bio.SeqIO.write(records_to_save, args.output, "fasta")
    print(f"{len(records_to_save)} records left in alignment. Alignment with unique sequences written to {args.output}")


if __name__ == '__main__':
    main()
