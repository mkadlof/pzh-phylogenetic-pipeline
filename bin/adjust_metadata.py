#!/usr/bin/env python3

import argparse
import csv


def adjust_metada(input_file: str) -> None:
    segments = ['chr1_PB2|',
                'chr2_PB1|',
                'chr3_PA|',
                'chr4_HA|',
                'chr5_NP|',
                'chr6_NA|',
                'chr7_MP|',
                'chr8_NS|', ]

    with open(input_file, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)
    adjusted_metadata = []
    for row in rows:
        for segment in segments:
            new_record = row.copy()
            new_record['strain'] = segment + row['strain']
            adjusted_metadata.append(new_record)
    # save the adjusted metadata to a new file
    output_file = input_file.replace('.tsv', '_adjusted.tsv')
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(adjusted_metadata)
    print(f"Adjusted metadata saved to {output_file}")


def main():
    parser = argparse.ArgumentParser("This scripts copies multiplies each metadata records by"
                                     " influenza segment to mach IDs in fasta files.")
    parser.add_argument('input_file', help='Path to the input file')
    args = parser.parse_args()
    adjust_metada(args.input_file)


if __name__ == '__main__':
    main()
