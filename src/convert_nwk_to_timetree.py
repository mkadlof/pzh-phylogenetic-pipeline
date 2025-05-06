"""
Scirpit replaces original evolutionary distances with date-based branch lengths
proposed by treetime
"""

import click
from Bio import Phylo
import json
from typing import Dict


def prep_ndates_dict(branches_json: str) -> Dict[str, float]:
    node_dates = {}
    with open(branches_json) as f:
        branches = json.load(f)
    for name, data in branches['nodes'].items():
        node_dates[name] = data['numdate']
    return node_dates


def read_tree(tree_file: str):
    return Phylo.read(tree_file, "newick")


@click.command()
@click.option('--tree', help='[INPUT] a newick file with a tree')
@click.option('--branches', help='[INPUT] a json file with a dictionary of node names and their dates')
@click.option('--output', help='[OUTPUT] a newick file with a tree')
def main(tree, branches, output):
    tree = read_tree(tree)
    node_dates = prep_ndates_dict(branches)

    for clade in tree.find_clades():
        if clade.name in node_dates:
            clade.numdate = node_dates[clade.name]
        else:
            clade.numdate = None

    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            child.parent = parent

    # Replace evolutionary distances with date-based branch lengths
    for clade in tree.find_clades(order="postorder"):
        if hasattr(clade, "parent") and clade.parent:
            parent_date = getattr(clade.parent, "numdate", None)
            child_date = getattr(clade, "numdate", None)
            if parent_date is not None and child_date is not None:
                clade.branch_length = child_date - parent_date
            else:
                clade.branch_length = 0.0


    Phylo.write(tree, output, "newick")


if __name__ == "__main__":
    main()