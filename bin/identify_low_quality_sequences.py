#!/usr/bin/env python3

import argparse
import csv
import os
from typing import Dict, List

def restricted_float(x):
    """This function is used in argparse to validate the float value to be between 0.0 and 1.0"""
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def validate_Ns_are_below_threshold(row, threshold=0.02):
    n_proportion = int(row['N']) / int(row['length'])
    valid = n_proportion <= threshold
    print(f"{row['strain']} N_proportion: {n_proportion} <= {threshold} Result: {'valid' if valid else 'invalid'}")
    return valid


def validate_ambiguities_are_below_threshold(row, threshold=0):
    ambiguity_proportion = int(row['other_IUPAC']) / int(row['length'])
    valid = ambiguity_proportion <= threshold
    print(f"{row['strain']} Ambiguity_proportion: {ambiguity_proportion} <= {threshold} Result: {'valid' if valid else 'invalid'}")
    return valid

# def validate_ambiguities_are_below_threshold(row, threshold=0):
#     return int(row['other_IUPAC']) / int(row['length']) <= threshold


def is_valid(row, args):
    return validate_Ns_are_below_threshold(row, args.threshold_Ns) and validate_ambiguities_are_below_threshold(row, args.threshold_ambiguities)


def write_strains(strains: List[Dict[str, str]], output):
    with open(output, 'w') as f:
        f.write('\n'.join(row['strain'] for row in strains))


def main():
    parser = argparse.ArgumentParser(description='Identify low quality sequences')
    parser.add_argument('input', type=str, help='Input CSV augur index file')
    parser.add_argument('--output_dir', default="results", type=str, help='Output file with valid strains')
    parser.add_argument('--threshold_Ns', default=0.02, type=restricted_float, help='Threshold for Ns validation')
    parser.add_argument('--threshold_ambiguities', default=0, type=restricted_float, help='Threshold for ambiguities validation')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    valid_strains = []
    invalid_strains = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        all_strains = 0
        valid = 0
        for row in reader:
            all_strains += 1
            if is_valid(row, args):
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
