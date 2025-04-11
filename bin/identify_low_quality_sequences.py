#!/usr/bin/env python3

import argparse
import csv
import os
from typing import Dict, List


def validate_Ns_are_below_threshold(row, threshold=0.02):
    return int(row['N']) / int(row['length']) <= threshold


def validate_ambiguities_are_below_threshold(row, threshold=0):
    return int(row['other_IUPAC']) / int(row['length']) <= threshold


def is_valid(row):
    return validate_Ns_are_below_threshold(row) and validate_ambiguities_are_below_threshold(row)


def write_strains(strains: List[Dict[str, str]], output):
    with open(output, 'w') as f:
        f.write('\n'.join(row['strain'] for row in strains))


def main():
    parser = argparse.ArgumentParser(description='Identify low quality sequences')
    parser.add_argument('input', type=str, help='Input CSV augur index file')
    parser.add_argument('--output_dir', default="results", type=str, help='Output file with valid strains')
    parser.add_argument('--threshold_Ns', default=0.02, type=float, help='Threshold for Ns validation')
    parser.add_argument('--threshold_ambiguities', default=0, type=float, help='Threshold for ambiguities validation')
    args = parser.parse_args()

    valid_strains = []
    invalid_strains = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        next(reader)
        all_strains = 0
        valid = 0
        for row in reader:
            all_strains += 1
            if is_valid(row):
                valid += 1
                valid_strains.append(row)
            else:
                invalid_strains.append(row)
    write_strains(valid_strains, os.path.join(args.output_dir, 'valid_strains.txt'))
    write_strains(invalid_strains, os.path.join(args.output_dir, 'invalid_strains.txt'))
    print(f"Sequence quality validation:\n"
          f"Valid strains: {valid}, Dropped strains: {all_strains - valid}")


if __name__ == '__main__':
    main()
