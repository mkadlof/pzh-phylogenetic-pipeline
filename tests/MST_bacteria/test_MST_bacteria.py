import numpy as np
import pandas as pd
import pytest
import subprocess
from calculate_allelic_distance_and_plot_MST import (
    calculate_distance,
    calculte_distance_matrix,
    calculate_mst,
    load_profiles_to_dict,
    extract_profiles,
    visualize_mst
)

# -------------------------------------------------------------------
# 🧩 UNIT TESTS
# -------------------------------------------------------------------


def test_load_profiles_to_dict(tmp_path):
    public_path = "profiles_test.list"
    local_path = "profiles_local_test.list"

    pub_ST = load_profiles_to_dict(public_path)
    local_ST = load_profiles_to_dict(local_path)

    # Verify allele vectors
    np.testing.assert_array_equal(pub_ST['1'][:10], np.array(['1', 2., 3., 1., 3., 4., 3., 4., 3., 3.], dtype=object))
    np.testing.assert_array_equal( local_ST['local_3'][:10], np.array(['local_3', 2., 34., 4., 10., 9., 12., 1., 16., 53.], dtype=object) )



def test_extract_profiles(tmp_path):
    """
    Integration test: verify that extract_profiles() correctly merges
    public and local cgMLST profiles using real data files.
    """
    public_path = "profiles_test.list"
    local_path = "profiles_local_test.list"

    # --- Load both sets ---
    pub_ST_dict = load_profiles_to_dict(public_path)
    local_ST_dict = load_profiles_to_dict(local_path)

    # --- Extract single STs ---
    list_of_allels_public = extract_profiles(
        profiles_public=pub_ST_dict, profiles_local=local_ST_dict, selected_sts=['1']
    )
    list_of_allels_local = extract_profiles(
        profiles_public=pub_ST_dict, profiles_local=local_ST_dict, selected_sts=['local_3']
    )

    list_of_allels_both = extract_profiles(
        profiles_public=pub_ST_dict, profiles_local=local_ST_dict, selected_sts=['1', 'local_3']
    )


    # --- Verify content ---
    np.testing.assert_array_equal(
        list_of_allels_public[0, :10],
        np.array(['1', 2., 3., 1., 3., 4., 3., 4., 3., 3.], dtype=object)
    )

    np.testing.assert_array_equal(
        list_of_allels_local[0, :10],
        np.array(['local_3', 2., 34., 4., 10., 9., 12., 1., 16., 53.], dtype=object)
    )

    np.testing.assert_array_equal(
        list_of_allels_both[:2, :10],
        np.array([['1', 2., 3., 1., 3., 4., 3., 4., 3., 3.],
                  ['local_3', 2., 34., 4., 10., 9., 12., 1., 16., 53. ]], dtype=object)
    )


def test_calculate_distance_matrix_real_data(tmp_path):
    """
    Integration test: calculate allelic distance matrix for real profiles
    (public + local), using only the first 9 loci for simplicity.
    """
    public_path = "profiles_test.list"
    local_path = "profiles_local_test.list"

    # --- Load and merge ---
    pub_ST = load_profiles_to_dict(public_path)
    local_ST = load_profiles_to_dict(local_path)
    all_STs = list(pub_ST.keys()) + list(local_ST.keys())
    mat = extract_profiles(pub_ST, local_ST, all_STs)
    mat = mat[0:9, :10]


    # --- Compute distance matrix ---
    dist = calculte_distance_matrix(mat, cpus=1)

    assert dist[0,0] == 0
    assert dist[0, 1] == 9
    assert dist[0, 7] == 9
    assert dist[0, 4] == 8
    assert dist[6, 7] == 8


def test_calculate_distance_matrix_caps_workers_to_profile_count():
    """Ensure more requested workers than profiles do not create empty jobs."""
    mat = np.array([
        ["1", 1.0, 2.0],
        ["2", 2.0, 3.0],
    ], dtype=object)

    dist = calculte_distance_matrix(mat, cpus=3)

    assert dist.shape == (2, 2)
    assert np.array_equal(np.diag(dist), np.zeros(2, dtype=np.int16))


def test_visualize_mst_runs(tmp_path):
    """Ensure MST visualization runs and produces an HTML file."""
    edges = np.array([["A", "B", 1], ["B", "C", 2]], dtype=object)
    counts = {"A": 1, "B": 2, "C": 1}
    color_map = {"A": "1", "B": "1", "C": "2"}
    sample_map = {"A": ["s1"], "B": ["s2"], "C": ["s3"]}
    out_html = tmp_path / "mst.html"
    visualize_mst(edges, counts, color_map=color_map,
                  color_label="HC10", output_html=str(out_html),
                  sample_map=sample_map)
    assert out_html.exists()
    assert out_html.stat().st_size > 1000

def test_end_to_end_script(tmp_path):
    """
    End-to-end test: simulate full workflow of cgMLST distance + MST generation.
    Uses real input data, produces both TSV and HTML outputs, and validates results.
    """
    public_path = "profiles_test.list"
    local_path = "profiles_local_test.list"
    output_tsv = tmp_path / "dist_matrix.tsv"
    output_html = tmp_path / "mst_plot.html"

    # Step 1 — Load and merge
    pub_ST = load_profiles_to_dict(public_path)
    local_ST = load_profiles_to_dict(local_path)
    all_STs = list(pub_ST.keys()) + list(local_ST.keys())
    mat = extract_profiles(pub_ST, local_ST, all_STs)
    mat_reduced = np.column_stack([mat[:, 0], mat[:, 1:10].astype(float)])

    # Step 2 — Compute allelic distances
    dist = calculte_distance_matrix(mat_reduced, cpus=1)
    pd.DataFrame(dist, index=all_STs, columns=all_STs).to_csv(output_tsv, sep="\t")

    # Step 3 — Compute MST edges
    edges = calculate_mst(dist, all_STs)

    # Step 4 — Define metadata for visualization
    counts = {st: 1 for st in all_STs}
    color_map = {st: "1" if st.startswith("local_") else "0" for st in all_STs}
    sample_map = {st: [f"sample_{st}"] for st in all_STs}

    # Step 5 — Visualize MST
    visualize_mst(
        edges,
        counts,
        color_map=color_map,
        color_label="HC10",
        output_html='test.html',
        sample_map=sample_map
    )

    # Step 6 — Validate outputs
    assert output_tsv.exists(), "Distance matrix file was not created."
    assert output_tsv.stat().st_size > 100, "Distance matrix seems empty."
    assert output_html.exists(), "MST HTML plot was not created."
    assert output_html.stat().st_size > 2000, "HTML output too small — possibly empty graph."

    # Step 7 — Sanity check on MST distances
    assert np.allclose(dist, dist.T), "Distance matrix is not symmetric."
    assert np.all(np.diag(dist) == 0), "Diagonal of distance matrix is not zero."
    assert np.any(dist > 0), "All distances are zero — likely failed parsing."

    print(f"\n✅ End-to-end workflow executed successfully: {output_html}")
