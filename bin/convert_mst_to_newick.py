#!/usr/bin/env python3

"""Convert a cgMLST minimum spanning tree edge list to Newick.

The MST vertices are cgMLST sequence types, whereas Microreact associates tree
tips with metadata rows by sample identifier.  Each sequence-type vertex is
therefore represented by an internal clade with one zero-length leaf per sample
assigned to that sequence type.

The tree is rooted at its weighted graph centre: the sequence type whose
maximum shortest-path distance to every other sequence type is minimal.  This
root is chosen only to produce a balanced and deterministic visualization; it
does not imply ancestry.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import click
import networkx as nx
from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree


def load_sample_map(
    metadata_path: str | Path,
    sample_column: str = "strain",
    st_column: str = "cgMLST",
) -> dict[str, list[str]]:
    """Return a mapping from cgMLST sequence type to sample identifiers."""
    sample_map: dict[str, list[str]] = defaultdict(list)
    seen_samples: set[str] = set()
    unassigned: list[str] = []

    with open(metadata_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Metadata file does not contain a header.")

        missing_columns = {
            column for column in (sample_column, st_column)
            if column not in reader.fieldnames
        }
        if missing_columns:
            raise ValueError(
                "Metadata is missing required column(s): "
                + ", ".join(sorted(missing_columns))
            )

        for line_number, row in enumerate(reader, start=2):
            sample = (row.get(sample_column) or "").strip()
            st = (row.get(st_column) or "").strip()
            if not sample:
                raise ValueError(
                    f"Empty {sample_column!r} value in metadata line {line_number}."
                )
            if sample in seen_samples:
                raise ValueError(
                    f"Sample identifier {sample!r} occurs more than once in metadata."
                )
            seen_samples.add(sample)

            # The MST itself ignores samples without a sequence type, so they are
            # reported and skipped here instead of aborting the whole run.
            if not st:
                unassigned.append(sample)
                continue

            sample_map[st].append(sample)

    if unassigned:
        click.echo(
            f"Warning: {len(unassigned)} sample(s) without a {st_column} value were "
            "excluded from the MST Newick: " + ", ".join(sorted(unassigned))
        )

    if not sample_map:
        raise ValueError(
            f"Metadata does not contain any sample with a {st_column} value."
        )

    return dict(sample_map)


def load_mst_graph(
    edge_list_path: str | Path,
    sequence_types: Iterable[str],
) -> nx.Graph:
    """Load the MST edge list and validate it against metadata sequence types."""
    expected_nodes = set(sequence_types)
    graph = nx.Graph()

    with open(edge_list_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source", "target", "distance"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "MST edge list must contain the columns: source, target, distance."
            )

        for line_number, row in enumerate(reader, start=2):
            source = (row.get("source") or "").strip()
            target = (row.get("target") or "").strip()
            if not source or not target:
                raise ValueError(
                    f"Empty source or target in MST edge-list line {line_number}."
                )

            unknown = {source, target} - expected_nodes
            if unknown:
                raise ValueError(
                    "MST edge list contains sequence type(s) absent from metadata: "
                    + ", ".join(sorted(unknown))
                )

            try:
                distance = float(row["distance"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid distance in MST edge-list line {line_number}: "
                    f"{row.get('distance')!r}."
                ) from exc

            if distance < 0:
                raise ValueError(
                    f"Negative MST distance in edge-list line {line_number}."
                )
            if graph.has_edge(source, target):
                raise ValueError(
                    f"Duplicate MST edge between {source!r} and {target!r}."
                )

            graph.add_edge(source, target, weight=distance)

    if graph.number_of_nodes() == 0:
        # A single sequence type yields an MST without any edge.
        graph.add_nodes_from(expected_nodes)

    if not nx.is_forest(graph):
        raise ValueError(
            "MST edge list must be acyclic "
            f"(nodes={graph.number_of_nodes()}, edges={graph.number_of_edges()})."
        )

    unplaced = expected_nodes - set(graph)
    if unplaced:
        click.echo(
            f"Warning: {len(unplaced)} sequence type(s) present in metadata are "
            "absent from the MST and their samples are excluded from the Newick: "
            + ", ".join(sorted(unplaced))
        )

    return graph


def choose_weighted_center(graph: nx.Graph) -> tuple[str, float]:
    """Choose a deterministic vertex minimizing weighted eccentricity."""
    eccentricities: dict[str, float] = {}

    for node, distances in nx.all_pairs_dijkstra_path_length(
        graph, weight="weight"
    ):
        eccentricities[node] = max(distances.values(), default=0.0)

    minimum_eccentricity = min(eccentricities.values())
    candidates = [
        node
        for node, eccentricity in eccentricities.items()
        if eccentricity == minimum_eccentricity
    ]
    root = min(candidates, key=str)
    return root, minimum_eccentricity


def graph_to_newick_clade(
    graph: nx.Graph,
    sample_map: dict[str, list[str]],
    root: str,
) -> Clade:
    """Orient the MST from ``root`` and expand ST vertices into sample tips."""
    if root not in graph:
        raise ValueError(f"Root {root!r} is absent from the MST.")

    parents: dict[str, str | None] = {root: None}
    children: dict[str, list[str]] = defaultdict(list)
    traversal_order: list[str] = []
    stack = [root]

    while stack:
        node = stack.pop()
        traversal_order.append(node)
        node_children = sorted(
            (
                neighbour
                for neighbour in graph.neighbors(node)
                if neighbour != parents[node]
            ),
            key=str,
        )
        children[node] = node_children
        for child in reversed(node_children):
            parents[child] = node
            stack.append(child)

    clades: dict[str, Clade] = {}
    for node in reversed(traversal_order):
        sample_tips = [
            Clade(name=sample, branch_length=0.0)
            for sample in sorted(sample_map[node])
        ]
        child_clades = []
        for child_node in children[node]:
            child = clades[child_node]
            child.branch_length = float(graph.edges[node, child_node]["weight"])
            child_clades.append(child)

        clades[node] = Clade(
            name=f"cgMLST_{node}",
            clades=sample_tips + child_clades,
        )

    return clades[root]


def graph_to_newick_tree(
    graph: nx.Graph,
    sample_map: dict[str, list[str]],
) -> tuple[Tree, str, float]:
    """Build the Newick tree and report the root vertex and its eccentricity."""
    components = [graph.subgraph(nodes) for nodes in nx.connected_components(graph)]
    centers = [choose_weighted_center(component) for component in components]

    # Components are ordered by decreasing size so that the reported root belongs
    # to the largest one; the root name breaks size ties deterministically.
    ordered = sorted(
        zip(components, centers),
        key=lambda item: (-item[0].number_of_nodes(), str(item[1][0])),
    )

    if len(ordered) == 1:
        component, (root, eccentricity) = ordered[0]
        clade = graph_to_newick_clade(component, sample_map, root)
        return Tree(root=clade, rooted=True), root, eccentricity

    # Sequence types whose allelic distance rounds down to zero are dropped by the
    # sparse MST implementation, which splits the MST into a forest.  Grafting the
    # component roots onto a shared zero-length node keeps the file a valid single
    # tree without inventing distances.
    click.echo(
        f"Warning: the MST consists of {len(ordered)} disconnected component(s); "
        "their roots are joined by zero-length branches."
    )
    component_clades = []
    for component, (component_root, _) in ordered:
        clade = graph_to_newick_clade(component, sample_map, component_root)
        clade.branch_length = 0.0
        component_clades.append(clade)

    root, eccentricity = ordered[0][1]
    return Tree(root=Clade(clades=component_clades), rooted=True), root, eccentricity


def convert_mst_to_newick(
    mst_edges: str | Path,
    metadata: str | Path,
    output: str | Path,
    sample_column: str = "strain",
    st_column: str = "cgMLST",
) -> tuple[str, float]:
    """Convert input files to Newick and return the selected root information."""
    sample_map = load_sample_map(metadata, sample_column, st_column)
    graph = load_mst_graph(mst_edges, sample_map.keys())
    tree, root, eccentricity = graph_to_newick_tree(graph, sample_map)
    Phylo.write(tree, output, "newick")
    return root, eccentricity


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--mst-edges",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="MST edge list TSV with source, target and distance columns.",
)
@click.option(
    "--metadata",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Sample metadata TSV.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output Newick file.",
)
@click.option("--sample-column", default="strain", show_default=True)
@click.option("--st-column", default="cgMLST", show_default=True)
def main(
    mst_edges: Path,
    metadata: Path,
    output: Path,
    sample_column: str,
    st_column: str,
) -> None:
    """Convert a cgMLST MST edge list into a sample-level Newick tree."""
    root, eccentricity = convert_mst_to_newick(
        mst_edges=mst_edges,
        metadata=metadata,
        output=output,
        sample_column=sample_column,
        st_column=st_column,
    )
    click.echo(
        f"MST Newick saved to {output}; weighted-centre root={root}, "
        f"eccentricity={eccentricity:g} allelic differences."
    )


if __name__ == "__main__":
    main()
