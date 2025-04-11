import os

import pytest
from Bio import SeqIO

from find_identical_seqences import find_identical_sequences
from find_identical_seqences import write_identical_sequences, write_unique_fasta


@pytest.fixture
def expected_values():
    return {"test-1.fasta": [["seq1", "seq2"], ["seq4", "seq5"]],
            "test-2.fasta": [["a1", "a2"]],
            "test-3.fasta": [["x", "z"]],
            "test-4.fasta": [],
            "test-5.fasta": [["id1", "id2"]]}


test_cases = ["test-1.fasta", "test-2.fasta", "test-3.fasta", "test-4.fasta", "test-5.fasta"]


@pytest.mark.parametrize("filename", test_cases)
def test_find_identical_sequences_file(filename, expected_values):
    fasta_path = os.path.join(os.path.dirname(__file__), "test_fasta", filename)
    result = find_identical_sequences(fasta_path)
    sorted_result = [sorted(group) for group in result]
    expected = [sorted(group) for group in expected_values[filename]]
    assert sorted(sorted_result) == sorted(expected)


@pytest.mark.parametrize("test_case", test_cases)
def test_write_identical_sequences(test_case, expected_values, tmp_path):
    groups = expected_values[test_case]
    write_identical_sequences(groups, str(tmp_path))
    output_file = str(tmp_path) + "_ident_seq.csv"
    with open(output_file) as f:
        lines = f.read().splitlines()
    expected_lines = [','.join(sorted(group)) for group in groups]
    assert lines == expected_lines


@pytest.mark.parametrize("test_case", test_cases)
def test_write_unique_fasta(test_case, expected_values, tmp_path):
    """Test if we can write fasta file containing unique sequences."""
    groups = expected_values[test_case]
    fasta_path = os.path.join(os.path.dirname(__file__), "test_fasta", test_case)
    write_unique_fasta(groups, fasta_path, str(tmp_path))

    output_file = str(tmp_path) + f"_unique.fasta"
    records = list(SeqIO.parse(output_file, "fasta"))
    ids = {r.id for r in records}
    if test_case == "test-1.fasta":
        assert len(records) == 3
        assert "seq1" in ids
        assert "seq2" not in ids
        assert "seq3" in ids
        assert "seq4" in ids
        assert "seq5" not in ids
    elif test_case == "test-2.fasta":
        assert len(records) == 2
        assert "a1" in ids
        assert "a2" not in ids
        assert "b1" in ids
    elif test_case == "test-3.fasta":
        assert len(records) == 2
        assert "x" in ids
        assert "y" in ids
        assert "z" not in ids
    elif test_case == "test-4.fasta":
        assert len(records) == 3
        assert "s1" in ids
        assert "s2" in ids
        assert "s3" in ids
    elif test_case == "test-5.fasta":
        assert len(records) == 1
        assert "id1" in ids
        assert "id2" not in ids
