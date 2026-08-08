import sys
from pathlib import Path

import pytest
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).parents[2] / "bin"))

import convert_mst_to_newick as MODULE


def test_converts_st_vertices_to_zero_length_sample_tips(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "strain\tcgMLST\n"
        "sample_A1\tA\n"
        "sample_A2\tA\n"
        "sample_B\tB\n"
        "sample_C\tC\n",
        encoding="utf-8",
    )
    edges = tmp_path / "mst.tsv"
    edges.write_text(
        "source\ttarget\tdistance\n"
        "A\tB\t2\n"
        "B\tC\t8\n",
        encoding="utf-8",
    )
    output = tmp_path / "mst.nwk"

    root, eccentricity = MODULE.convert_mst_to_newick(
        mst_edges=edges,
        metadata=metadata,
        output=output,
    )

    assert root == "B"
    assert eccentricity == 8

    tree = Phylo.read(output, "newick")
    terminals = {tip.name: tip for tip in tree.get_terminals()}
    assert set(terminals) == {"sample_A1", "sample_A2", "sample_B", "sample_C"}
    assert all(tip.branch_length == 0 for tip in terminals.values())
    assert tree.distance("sample_A1", "sample_A2") == 0
    assert tree.distance("sample_A1", "sample_B") == 2
    assert tree.distance("sample_A1", "sample_C") == 10

    internal_names = {
        clade.name for clade in tree.get_nonterminals() if clade.name is not None
    }
    assert internal_names == {"cgMLST_A", "cgMLST_B", "cgMLST_C"}


def test_supports_a_single_sequence_type(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "strain\tcgMLST\n"
        "sample_1\tlocal_1\n"
        "sample_2\tlocal_1\n",
        encoding="utf-8",
    )
    edges = tmp_path / "mst.tsv"
    edges.write_text("source\ttarget\tdistance\n", encoding="utf-8")
    output = tmp_path / "mst.nwk"

    root, eccentricity = MODULE.convert_mst_to_newick(edges, metadata, output)

    assert root == "local_1"
    assert eccentricity == 0
    tree = Phylo.read(output, "newick")
    assert {tip.name for tip in tree.get_terminals()} == {"sample_1", "sample_2"}
    assert tree.distance("sample_1", "sample_2") == 0


def test_weighted_center_tie_is_deterministic():
    graph = MODULE.nx.Graph()
    graph.add_edge("C", "D", weight=1)
    graph.add_edge("B", "C", weight=9)
    graph.add_edge("A", "B", weight=1)

    root, eccentricity = MODULE.choose_weighted_center(graph)

    assert root == "B"
    assert eccentricity == 10


def test_sequence_types_missing_from_the_mst_are_skipped(tmp_path):
    # calculate_allelic_distance_and_plot_MST.py only warns when a cgMLST profile
    # cannot be resolved, so such sequence types never reach the edge list.
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "strain\tcgMLST\n"
        "sample_A\tA\n"
        "sample_B\tB\n"
        "sample_unresolved\tZ\n",
        encoding="utf-8",
    )
    edges = tmp_path / "mst.tsv"
    edges.write_text("source\ttarget\tdistance\nA\tB\t3\n", encoding="utf-8")
    output = tmp_path / "mst.nwk"

    root, _ = MODULE.convert_mst_to_newick(edges, metadata, output)

    assert root in {"A", "B"}
    tree = Phylo.read(output, "newick")
    assert {tip.name for tip in tree.get_terminals()} == {"sample_A", "sample_B"}


def test_samples_without_a_sequence_type_are_skipped(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "strain\tcgMLST\n"
        "sample_A\tA\n"
        "sample_B\tB\n"
        "sample_missing\t\n",
        encoding="utf-8",
    )
    edges = tmp_path / "mst.tsv"
    edges.write_text("source\ttarget\tdistance\nA\tB\t3\n", encoding="utf-8")
    output = tmp_path / "mst.nwk"

    MODULE.convert_mst_to_newick(edges, metadata, output)

    tree = Phylo.read(output, "newick")
    assert {tip.name for tip in tree.get_terminals()} == {"sample_A", "sample_B"}


def test_disconnected_mst_is_joined_with_zero_length_branches():
    """A forest must still be written as a single valid Newick tree."""
    graph = MODULE.nx.Graph()
    graph.add_edge("A", "B", weight=4)
    graph.add_node("C")
    sample_map = {"A": ["sample_A"], "B": ["sample_B"], "C": ["sample_C"]}

    tree, root, eccentricity = MODULE.graph_to_newick_tree(graph, sample_map)

    assert root == "A"
    assert eccentricity == 4
    terminals = {tip.name for tip in tree.get_terminals()}
    assert terminals == {"sample_A", "sample_B", "sample_C"}
    assert tree.distance("sample_A", "sample_B") == 4
    assert tree.distance("sample_A", "sample_C") == 0


def test_edge_list_with_unknown_sequence_type_is_rejected(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("strain\tcgMLST\nsample_A\tA\n", encoding="utf-8")
    edges = tmp_path / "mst.tsv"
    edges.write_text("source\ttarget\tdistance\nA\tB\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="absent from metadata"):
        MODULE.convert_mst_to_newick(edges, metadata, tmp_path / "mst.nwk")
