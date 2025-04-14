#!/usr/bin/env python3

import argparse
from typing import List
from collections import defaultdict

import Bio
from Bio.SeqIO import SeqRecord


def insert_duplicates(fasta_file: str, duplicated_ids_file: str):
    """
    Insert duplicated sequences into a FASTA file based on a CSV file with duplicated IDs.
    """
    sequences_to_insert = defaultdict(list)
    with open(duplicated_ids_file, 'r') as f:
        for line in f:
            # Extract the duplicated IDs from the line
            ids = line.strip().split(",")
            first_id = ids[0]
            for next_id in ids[1:]:
                sequences_to_insert[first_id].append(next_id)

    # Read the FASTA file and insert the duplicated sequences
    records = []
    counter = 0
    for record in Bio.SeqIO.parse(fasta_file, "fasta"):
        counter += 1
        records.append(record)
        if record.id in sequences_to_insert:
            for next_id in sequences_to_insert[record.id]:
                new_record = record[:]
                new_record.id = next_id
                new_record.description = ''
                records.append(new_record)
    print(f"{counter} records read from input alignment.")
    return records


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
    parser.add_argument("--action", choices=["insert", "remove"], help="Action to perform: insert or remove duplicates")
    args = parser.parse_args()

    if args.action == "insert":
        # Insert duplicated sequences into the alignment
        records_to_save = insert_duplicates(args.alignment, args.duplicated_ids)
    elif args.action == "remove":
        records_to_save = remove_duplicates(args.alignment, args.duplicated_ids)
    else:
        raise ValueError("Invalid action. Use 'insert' or 'remove'.")
    # Write the unique sequences to the output FASTA file
    Bio.SeqIO.write(records_to_save, args.output, "fasta")
    print(f"{len(records_to_save)} records in alignment after {args.action}ing duplicates. Alignment with unique sequences written to {args.output}")


if __name__ == '__main__':
    main()
