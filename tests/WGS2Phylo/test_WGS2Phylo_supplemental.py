"""
Integration tests for WGS2Phylo.py's --supplemental-file handling.

These run the real CLI (generate_metadata) against the checked-in example
WGS output directories (data/example_data/WGS2phylo/results_*), which is
the same data the README's WGS2Phylo examples use, to confirm:

- HierCC levels HC0/HC2/HC5/HC10/HC20 are emitted in normal mode, with a
  graceful "Unknown" fallback for levels a given organism/scheme does not
  report (Campylobacter has no HC0/HC2/HC20).
- Any supplemental column outside the date/region/country/division/city
  whitelist (e.g. age, gender, or an arbitrary column like sample_source)
  is passed through as-is, with missing values written as "N/A".
- Extra column values are sanitized against spreadsheet formula injection.
- --extra-fields mode keeps exactly one copy of every column (no
  duplicates/collisions) with the new columns ordered ahead of the
  WGS-derived AMR/QC block.
"""

import csv
import sys
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(REPO_ROOT / "bin"))

from WGS2Phylo import generate_metadata  # noqa: E402

CAMPYLOBACTER_RESULTS = REPO_ROOT / "data" / "example_data" / "WGS2phylo" / "results_Campylobacter"
RSV_RESULTS = REPO_ROOT / "data" / "example_data" / "WGS2phylo" / "results_rsv"


def write_supplemental(path, header, rows):
    """Write a tab-separated supplemental metadata file."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path):
    """Read a tab-separated file into a list of dicts keyed by header."""
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run_generate_metadata(tmp_path, output_dir, organism, header, rows, extra_fields=False):
    """Invoke the real WGS2Phylo CLI and return the aggregated output rows."""
    supplemental_file = tmp_path / "supplemental.tsv"
    write_supplemental(supplemental_file, header, rows)
    output_prefix = tmp_path / "out" / "metadata"

    arguments = [
        "--output_dir", str(output_dir),
        "--organism", organism,
        "--supplemental-file", str(supplemental_file),
        "--output-prefix", str(output_prefix),
        "--without-fasta",
    ]
    if extra_fields:
        arguments.append("--extra-fields")

    result = CliRunner().invoke(generate_metadata, arguments)
    assert result.exit_code == 0, result.output
    return read_tsv(f"{output_prefix}.tsv")


def test_campylobacter_hiercc_levels_and_extra_columns(tmp_path):
    header = ["id", "date", "region", "country", "division", "city", "age", "gender"]
    rows = [
        ["EQA-Camp-25-02_S28_L001", "2023-01-01", "Europe", "Poland", "-", "Warsaw", "34", "F"],
        ["EQA-Camp-25-05_S31_L001", "2023-01-02", "Europe", "Poland", "-", "Warsaw", "", "M"],
        ["EQA-Camp-25-04_S30_L001", "2024-01-01", "Europe", "Germany", "-", "Berlin", "58", "F"],
        ["EQA-Camp-25-01_S27_L001", "2000-01-02", "Europe", "Germany", "-", "Berlin", "27", "M"],
        ["EQA-Camp-25-03_S29_L001", "2000-12-12", "Europe", "Germany", "-", "Munchen", "45", "F"],
    ]
    rows_out = run_generate_metadata(tmp_path, CAMPYLOBACTER_RESULTS, "campylobacter", header, rows)
    by_strain = {r["strain"]: r for r in rows_out}

    # Campylobacter's cgMLST scheme does not report HC0/HC2/HC20 for this
    # dataset; missing levels must gracefully fall back to "Unknown", unlike
    # HC5/HC10 which ARE reported.
    row = by_strain["EQA-Camp-25-02_S28_L001"]
    assert row["HC0"] == "Unknown"
    assert row["HC2"] == "Unknown"
    assert row["HC20"] == "Unknown"
    assert row["HC5"] != "Unknown"
    assert row["HC10"] != "Unknown"

    # Extra supplemental columns pass through as-is; missing values become "N/A".
    assert row["age"] == "34"
    assert row["gender"] == "F"
    assert by_strain["EQA-Camp-25-05_S31_L001"]["age"] == "N/A"


def test_rsv_arbitrary_column_and_injection_sanitization(tmp_path):
    header = ["id", "date", "region", "country", "division", "city", "age", "gender", "sample_source"]
    rows = [
        ["SRR27383080_pass", "2000-01-02", "Europe", "Poland", "-", "Warsaw", "12", "M", "clinical"],
        ["SRR27383087_pass", "2023-01-02", "Europe", "Poland", "-", "Warsaw", "67", "F", "environmental"],
        ["SRR27383088_pass", "2024-01-01", "Europe", "Germany", "-", "Berlin", "", "M", "clinical"],
        ["SRR27383089_pass", "2000-01-02", "Europe", "Germany", "-", "Berlin", "39", "F", "=cmd|'/c calc'!A1"],
        ["SRR27383090_pass", "2000-12-12", "Europe", "Germany", "-", "Munchen", "51", "M", "clinical"],
    ]
    rows_out = run_generate_metadata(tmp_path, RSV_RESULTS, "rsv", header, rows)
    by_strain = {r["strain"]: r for r in rows_out}

    # "sample_source" is not part of any whitelist anywhere in WGS2Phylo -- it
    # must still be passed through untouched, proving the script doesn't
    # filter arbitrary supplemental columns.
    assert by_strain["SRR27383080_pass"]["sample_source"] == "clinical"
    assert by_strain["SRR27383080_pass"]["age"] == "12"

    # Missing value -> "N/A", and a spreadsheet-formula-injection payload is
    # defused with a leading apostrophe rather than passed through raw.
    assert by_strain["SRR27383088_pass"]["age"] == "N/A"
    assert by_strain["SRR27383089_pass"]["sample_source"] == "'=cmd|'/c calc'!A1"


def test_extra_fields_mode_keeps_new_columns_without_duplicates(tmp_path):
    header = ["id", "date", "region", "country", "division", "city", "age", "gender"]
    rows = [
        ["EQA-Camp-25-02_S28_L001", "2023-01-01", "Europe", "Poland", "-", "Warsaw", "34", "F"],
    ]
    rows_out = run_generate_metadata(
        tmp_path, CAMPYLOBACTER_RESULTS, "campylobacter", header, rows, extra_fields=True
    )
    fields = list(rows_out[0].keys())

    assert fields.count("age") == 1
    assert fields.count("gender") == 1
    for level in ("HC0", "HC2", "HC5", "HC10", "HC20"):
        assert fields.count(level) == 1
    # HierCC + supplemental columns must appear before the WGS --extra-fields
    # block (AMR/QC columns), matching the header assembly order.
    assert fields.index("gender") < fields.index("Azithromycin_opornos")
