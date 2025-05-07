#!/usr/bin/env python3


import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Rozszerz drzewo Newick o brakujące identyfikatory.")
    parser.add_argument("--tree", required=True, help="Plik z drzewem Newick (np. tree.nwk)")
    parser.add_argument("--ids", required=True, help="Plik CSV z identyfikatorami (pierwszy znany, reszta do dodania)")
    parser.add_argument("--out", help="Plik wyjściowy (jeśli nie podano, wypisuje na stdout)")
    return parser.parse_args()

def load_id_map(path):
    id_map = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) >= 2:
                id_map[parts[0]] = parts[1:]
    return id_map

def extend_newick(newick_str, id_map):
    for main_id, extra_ids in id_map.items():
        # Dopasuj np. A:0.1 → A:0.1,E:0.1,...
        newick_str = re.sub(
            rf"(?<!\w){main_id}:(\d+(?:\.\d+)?)",
            lambda m: ",".join([f"{name}:{m.group(1)}" for name in [main_id] + extra_ids]),
            newick_str
        )
    return newick_str

def main():
    args = parse_args()
    with open(args.tree) as f:
        newick_str = f.read().strip()
    id_map = load_id_map(args.ids)
    extended = extend_newick(newick_str, id_map)

    if args.out:
        with open(args.out, "w") as f:
            f.write(extended + "\n")
    else:
        print(extended)

if __name__ == "__main__":
    main()
