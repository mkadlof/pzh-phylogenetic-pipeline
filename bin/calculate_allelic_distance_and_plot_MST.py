import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.sparse.csgraph import minimum_spanning_tree
import click
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px

# ---------- CORE DISTANCE FUNCTIONS ----------

def calculate_distance(mat: np.ndarray, start_row: int, end_row: int, allowed_missing: float = 0.05) -> np.ndarray:
    """Calculate a fragment of the allelic distance matrix for rows [start_row, end_row)."""
    n_loci = mat.shape[1] - 1
    distance_matrix = np.zeros(((end_row - start_row), mat.shape[0]), dtype=np.int16)

    for i in range(start_row, end_row):
        out_index_row = i - start_row
        wektor_x = mat[i, 1:]
        ql = np.sum(wektor_x > 0)

        for j in range(mat.shape[0]):
            wektor_y = mat[j, 1:]
            rl = np.sum(wektor_y > 0)
            shared_nonzero = (wektor_x > 0) & (wektor_y > 0)

            al = 1e-4 + np.sum(shared_nonzero)
            ad = 1e-4 + np.sum(wektor_x[shared_nonzero] != wektor_y[shared_nonzero])

            ll = max(ql, rl) - allowed_missing * n_loci
            if ll > al:
                ad += ll - al
                al = ll

            distance_matrix[out_index_row, j] = np.int16(ad / al * n_loci + 0.5)

    return distance_matrix


def calculte_distance_matrix(mat: np.ndarray, cpus: int = 1, allowed_missing: float = 0.05) -> np.ndarray:
    """Compute the full N×N allelic distance matrix for a given profile matrix."""
    n_samples = mat.shape[0]
    if n_samples == 0:
        raise ValueError("Cannot calculate distances without cgMLST profiles.")
    if cpus < 1:
        raise ValueError("The number of CPU threads must be at least 1.")

    cpus = min(cpus, n_samples)
    distance_matrix = np.zeros((n_samples, n_samples), dtype=np.int16)
    row_splits = np.array_split(np.arange(n_samples), cpus)

    with Pool(cpus) as pool:
        jobs = [
            pool.apply_async(calculate_distance, (mat, rows[0], rows[-1] + 1, allowed_missing))
            for rows in row_splits
        ]
        start = 0
        for job in jobs:
            partial_matrix = job.get()
            end = start + partial_matrix.shape[0]
            distance_matrix[start:end, :] = partial_matrix
            start = end

    return distance_matrix


def calculate_mst(distance_matrix: np.ndarray, labels: list[str]) -> np.ndarray:
    """Compute the Minimum Spanning Tree (MST) from the allelic distance matrix."""
    # scipy interprets a zero weight as a missing edge, which would disconnect
    # sequence types whose allelic distance rounds down to zero. Offsetting every
    # candidate edge by one keeps the graph spanning without changing which edges
    # are optimal, so the offset is removed again from the reported distances.
    mat = np.triu(distance_matrix.astype(np.int64) + 1, 1)
    mst = minimum_spanning_tree(mat).toarray().astype(np.int64)
    edges = [(labels[i], labels[j], int(mst[i, j]) - 1)
             for i in range(mst.shape[0]) for j in range(mst.shape[1]) if mst[i, j] > 0]
    return np.array(edges, dtype=object)


# ---------- UTILITIES ----------

def load_profiles_to_dict(profiles_path: str) -> dict[str, np.ndarray]:
    """
    Load all profiles from a profiles.list file into a dictionary.
    Returns: {ST_ID: np.array([ST_ID, allele_1, allele_2, ...])}
    """
    profiles = {}
    with open(profiles_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 2:
                continue
            st_id = parts[0]
            try:
                alleles = [float(x) if x not in ("", "-", "NA") else 0.0 for x in parts[1:]]
            except ValueError:
                raise ValueError(f"Invalid allele values in profile for {st_id}: {parts[1:]}")
            # prepend ST ID as first element (kept as string for labeling)
            row = np.array([st_id] + alleles, dtype=object)
            profiles[st_id] = row
    return profiles

def extract_profiles(profiles_public: dict[str, np.ndarray],
                     profiles_local: dict[str, np.ndarray] | None,
                     selected_sts: list[str]) -> np.ndarray:
    """
    Combine allele profiles for cgMLST IDs from public and local sources.
    IDs starting with 'local_' are fetched from the local profile dictionary.
    """
    rows = []
    for st in selected_sts:
        if st.startswith("local_") and profiles_local and st in profiles_local:
            rows.append(profiles_local[st])
        elif st in profiles_public:
            rows.append(profiles_public[st])
        else:
            click.echo(f"Warning: cgMLST ID {st} not found in either profile source.")
    if not rows:
        raise ValueError("No cgMLST profiles were found for provided metadata.")
    return np.vstack(rows)


def visualize_mst(edges: np.ndarray, counts: dict[str, int],
                  color_map: dict[str, str] | None = None,
                  color_label: str | None = None,
                  output_html: str | None = None,
                  sample_map:  dict | None = None) -> go.Figure:
    """
    Create an interactive MST visualization (Plotly preferred, static matplotlib fallback).
    Node size is proportional to sample count, color depends on metadata attribute.
    """
    # --- Build MST graph ---
    G = nx.Graph()
    for s, t, d in edges:
        # print(f'{s}\t{t}\t{d}')
        # Layout weights are clamped from below so that sequence types separated by
        # zero allelic differences do not collapse onto each other.
        G.add_edge(s, t, weight=min(max(d, 1), 50), true_weight=d)

    # Optimize position of nodes
    pos = nx.kamada_kawai_layout(G, weight="weight")

    # --- Build edge coordinates and distance hover text ---
    edge_x, edge_y, edge_text = [], [], []
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        dist = data.get("true_weight", 0)
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        edge_text.append(f"{u}–{v}: {dist} allelic differences")

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=1.5, color="#888"),  # constant color
        hoverinfo="text",
        text=edge_text,
    )

    # --- Determine if color values are numeric or categorical ---

    # --- Determine coloring scheme ---
    if color_label and color_label.upper().startswith("HC"):
        # HierCC is categorical
        is_numeric_color = False
        qualitative_palette = px.colors.qualitative.Plotly

        unique_clusters = sorted(set(color_map.values()))
        cluster_color_map = {
            clu: qualitative_palette[i % len(qualitative_palette)]
            for i, clu in enumerate(unique_clusters)
        }

        get_node_color = lambda n: cluster_color_map.get(color_map.get(n, "NA"), "#A0A0A0")

    else:
        # determine dynamically if numeric or categorical
        unique_colors = set(color_map.values()) if color_map else set()
        is_numeric_color = all(
            str(v).replace('.', '', 1).isdigit() for v in unique_colors if v != "NA"
        )

        if is_numeric_color:
            get_node_color = lambda n: float(color_map.get(n)) if color_map.get(n) not in ("NA", "", None) else 0.0
        else:
            qualitative_palette = px.colors.qualitative.Safe
            unique_vals = sorted(set(color_map.values()))
            cat_color_map = {
                val: qualitative_palette[i % len(qualitative_palette)]
                for i, val in enumerate(unique_vals)
            }
            get_node_color = lambda n: cat_color_map.get(color_map.get(n, "NA"), "#A0A0A0")

    # --- Now build nodes with precomputed color mapping ---
    # Define max
    NODE_MIN_SIZE, NODE_MAX_SIZE = 20, 60

    min_c, max_c = min(counts.values()), max(counts.values())

    # NODE_SIZE_SCALE = 2
    node_x, node_y, node_text, node_size, node_color = [], [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        c = counts.get(node, 1)
        if max_c == min_c:
            size = (NODE_MIN_SIZE + NODE_MAX_SIZE) / 2
        else:
            size = NODE_MIN_SIZE + (c - min_c) / (max_c - min_c) * (NODE_MAX_SIZE - NODE_MIN_SIZE)
        node_size.append(size)

        # node_size.append(10 + np.sqrt(c) * 6 * NODE_SIZE_SCALE)

        node_color.append(get_node_color(node))
        color_value = color_map.get(node, "NA") if color_map else "NA"
        samples = ", ".join(sample_map.get(node, []))
        # Definiton of a text "hovering" over node -> cgMLST, numbe of sample + samples id
        node_text.append(
            f"cgMLST {node}<br>"
            f"Samples: {c}<br>"
            f"Sample IDs: {samples if samples else 'NA'}<br>"
            f"{color_label or 'Attribute'}: {color_value}"
        )


    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode="markers+text",
        text=[str(n) for n in G.nodes()],
        textposition="top center",
        hoverinfo="text",
        marker=dict(
            size=node_size,
            color=node_color,
            colorscale="Viridis" if is_numeric_color else None,
            showscale=is_numeric_color,
            line=dict(width=1, color="darkblue"),
            colorbar=dict(title=color_label if is_numeric_color else "")
        ),
        hovertext=node_text,
    )

    # Hide legend
    marker = dict(
        size=node_size,
        color=node_color,
        colorscale="Viridis" if is_numeric_color else None,
        showscale=False,  # hide colorbar for HierCC/categorical
        line=dict(width=1, color="darkblue"),
    )

    # --- Add invisible midpoint markers for edge hover text ---
    mnode_x, mnode_y, mnode_text = [], [], []
    for u, v, data in G.edges(data=True):
        # Midpoint of u-v
        x_u, y_u = pos[u]
        x_v, y_v = pos[v]
        mx = (x_u + x_v) / 2
        my = (y_u + y_v) / 2
        mnode_x.append(mx)
        mnode_y.append(my)
        w = data.get("true_weight", 0)
        #mnode_text.append(f"{u}–{v}: {w} allelic differences (edges with allelic distance above 50 have identical length)")
        mnode_text.append(f"{w} allelic differences")


    mnode_trace = go.Scatter(
        x=mnode_x,
        y=mnode_y,
        mode="markers",
        showlegend=False,
        hoverinfo="text",
        hovertext=mnode_text,
        # marker=dict(
        #     size=12,  # larger circle
        #     color="rgba(255,255,255,0)",  # transparent fill
        #     line=dict(
        #         width=2.5,
        #         color="rgba(100,100,100,0.5)"  # soft gray outer rim
        #     ),
        #     symbol="circle"
        # ),
        marker=dict(
            size=10,
            symbol="line-ns",  # “north–south” line
            color="rgba(100,100,100,0.6)",
            line=dict(width=2, color="rgba(80,80,80,0.8)")
        )
    )

    fig = go.Figure(
        data=[edge_trace, node_trace, mnode_trace],
        layout=go.Layout(
            title=f"cgMLST Minimum Spanning Tree (wezly pokolorowane wedlug wartosci z kolumny {color_label})",
            showlegend=False,
            hovermode="closest",
            hoverlabel=dict(
                bgcolor="rgba(240,240,240,0.9)",  # light gray background
                font_color="black",  # dark text for readability
                bordercolor="rgba(200,200,200,0.8)"  # subtle border
            ),
            margin=dict(b=0, l=0, r=0, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
            paper_bgcolor="white"
        )
    )

    # --- Create toggleable edge label annotations ---
    offset_scale = 0.03  # how far to offset labels perpendicular to edge
    annotations = []

    for u, v, data in G.edges(data=True):
        # edge direction and midpoint
        x_u, y_u = pos[u]
        x_v, y_v = pos[v]
        dx, dy = x_v - x_u, y_v - y_u
        length = np.hypot(dx, dy) or 1.0
        mx, my = (x_u + x_v) / 2, (y_u + y_v) / 2

        # perpendicular unit vector
        pos_x, pos_y = -dy / length, dx / length

        # flip side deterministically (helps when multiple edges share similar orientation)
        if hash(u) % 2 == 0:
            pos_x, pos_y = -pos_x, -pos_y

        # offset midpoint slightly
        mx_shift = mx + pos_x * offset_scale
        my_shift = my + pos_y * offset_scale

        annotations.append(
            dict(
                x=mx_shift,
                y=my_shift,
                text=f"{data.get('true_weight', 0)} allelic differences",
                showarrow=False,
                font=dict(size=12, color="gray"),
                bgcolor="rgba(245,245,245,0.6)",
                bordercolor="rgba(180,180,180,1)",
                borderpad=2,
                opacity=0.9
            )
        )

    # --- Add interactive buttons ---
    fig.update_layout(
        updatemenus=[{
            "buttons": [
                {
                    "label": "Show Edge Labels",
                    "method": "relayout",
                    "args": [{"annotations": annotations}]
                },
                {
                    "label": "Hide Edge Labels",
                    "method": "relayout",
                    "args": [{"annotations": []}]
                }
            ],
            "direction": "down",
            "x": 1.05,
            "y": 1.0,
            "showactive": True,
            "xanchor": "left",
            "yanchor": "top"
        }]
    )

    if output_html:
        fig.write_html(output_html, auto_open=True)
        click.echo(f"🌐 Interactive MST visualization saved to: {output_html}")


# ---------- CLI ----------

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--metadata", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Metadata TSV containing cgMLST, HC10, Serovar, etc.")
@click.option("--profiles", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Public profiles.list with cgMLST allele data.")
@click.option("--local-profiles", type=click.Path(exists=True, dir_okay=False),
              help="Optional local profiles file containing local_x allele definitions.")
@click.option("--color-by", default="HC5", show_default=True,
              help="Metadata column used for node coloring.")
@click.option("--plot", "plot_html", type=click.Path(dir_okay=False),
              help="Output path for interactive MST visualization (HTML).")
@click.option("--output", "output_path", type=click.Path(dir_okay=False),
              help="Optional output path for the distance matrix (TSV).")
@click.option("--mst-output", "mst_output", type=click.Path(dir_okay=False),
              help="Optional output path for the MST edge list (TSV).")
@click.option("-t", "--threads", default=1, show_default=True,
              help="Number of CPU threads to use.")
@click.option("-m", "--missing", default=0.05, show_default=True,
              help="Allowed fraction of missing loci (0–1).")
def main(metadata, profiles, local_profiles, color_by,
         plot_html, output_path, mst_output, threads, missing):
    """
    Compute the cgMLST allelic distance matrix and MST using:
    - metadata.tsv (contains cgMLST, HC10, etc.)
    - profiles.list (public allele definitions)
    - optionally local profiles (for local_x IDs)
    """
    meta = pd.read_csv(metadata, sep="\t", dtype=str)


    if "cgMLST" not in meta.columns:
        raise ValueError("Missing required column 'cgMLST' in metadata file!")

    counts = meta["cgMLST"].value_counts().to_dict()

    # Group samples by their CgMLST (to pass to MST visualization)
    sample_map = (
        meta.groupby("cgMLST")["strain"]
        .apply(list)
        .to_dict()
    )


    selected_sts = list(counts.keys())

    # Load public and local profile sets
    click.echo(f"📚 Loading public profiles from: {profiles}")
    profiles_public = load_profiles_to_dict(profiles)

    profiles_local = None
    if local_profiles:
        click.echo(f"📂 Loading local profiles from: {local_profiles}")
        profiles_local = load_profiles_to_dict(local_profiles)

    # Prepare color mapping
    if color_by in meta.columns:
        color_map = meta.groupby("cgMLST")[color_by].first().to_dict()
        click.echo(f"🎨 Coloring nodes by metadata column: {color_by}")
    else:
        click.echo(f"⚠️ Column '{color_by}' not found in metadata file; using single color.")
        color_map = {k: "NA" for k in selected_sts}

    # Extract allele data
    click.echo(f"🔎 Extracting profiles for {len(selected_sts)} cgMLST IDs...")
    mat = extract_profiles(profiles_public, profiles_local, selected_sts)

    # Compute allelic distance matrix
    click.echo("⚙️  Calculating allelic distance matrix...")
    dist = calculte_distance_matrix(mat, cpus=threads, allowed_missing=missing)
    labels = mat[:, 0].astype(str).tolist()

    if output_path:
        dist_df = pd.DataFrame(dist, index=labels, columns=labels)
        dist_df.to_csv(output_path, sep="\t", index=True)
        click.echo(f"💾 Distance matrix with labels saved to: {output_path}")

    # Compute MST
    click.echo("🌳 Building MST...")

    edges = calculate_mst(dist, labels)

    if mst_output:
        np.savetxt(mst_output, edges, fmt="%s", delimiter="\t",
                   header="source\ttarget\tdistance", comments="")
        click.echo(f"💾 MST edge list saved to: {mst_output}")

    # Visualize MST
    visualize_mst(edges, counts, color_map=color_map, color_label=color_by, output_html=plot_html, sample_map = sample_map)
    click.echo("✅ Completed successfully.")


if __name__ == "__main__":
    main()

    # root@5666e7777f6f:/dane# python calculate_allelic_distance.py --metadata metadata_expanded_escherichia.tsv --profiles /db/cgmlst/Escherichia/profiles.list --local-profiles /db/cgmlst/Escherichia/local/profiles_local.list --plot mlst_test_plot.html --output mlst_test_distance.tsv -t 10
