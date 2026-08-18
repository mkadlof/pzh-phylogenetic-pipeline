import json
import sys
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).parents[2]
TEMPLATE_PATH = REPO_ROOT / "data" / "microreact_config_bacteria.microreact"
sys.path.insert(0, str(REPO_ROOT / "bin"))

import prepare_json_for_microreact as MODULE


def tab_ids(node):
    found = []
    if node.get("type") == "tab":
        found.append(node.get("id"))
    for child in node.get("children", []):
        found.extend(tab_ids(child))
    return found


def tabsets(node):
    found = [node] if node.get("type") == "tabset" else []
    for child in node.get("children", []):
        found.extend(tabsets(child))
    return found


def tree_tabset(layout):
    matching = [
        tabset
        for tabset in tabsets(layout)
        if any(tab["id"].startswith("tree-") for tab in tabset["children"])
    ]
    assert len(matching) == 1, "trees are expected to share a single tabbed pane"
    return matching[0]


def run_generator(tmp_path, include_mst):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("strain\tcgMLST\nsample_1\t1\n", encoding="utf-8")
    classical = tmp_path / "classical.nwk"
    classical.write_text("(sample_1:0);\n", encoding="utf-8")
    rescaled = tmp_path / "rescaled.nwk"
    rescaled.write_text("(sample_1:0);\n", encoding="utf-8")
    mst = tmp_path / "mst.nwk"
    mst.write_text("(sample_1:0)cgMLST_1;\n", encoding="utf-8")
    output = tmp_path / "project.microreact"

    arguments = [
        "--input_json", str(TEMPLATE_PATH),
        "--classical_tree", str(classical),
        "--rescaled_tree", str(rescaled),
        "--metadata", str(metadata),
        "--project_name", "test",
        "--output", str(output),
    ]
    if include_mst:
        arguments.extend(["--mst_tree", str(mst)])

    result = CliRunner().invoke(MODULE.main, arguments)
    assert result.exit_code == 0, result.output
    return json.loads(output.read_text(encoding="utf-8"))


def test_adds_mst_file_tree_and_panel(tmp_path):
    project = run_generator(tmp_path, include_mst=True)

    assert "mst0" in project["files"]
    assert project["files"]["mst0"]["blob"].startswith(
        "data:application/octet-stream;base64, "
    )
    assert project["trees"]["tree-3"]["file"] == "mst0"
    assert "tree-3" in tab_ids(project["panes"]["model"]["layout"])


def test_all_trees_share_one_tabbed_pane(tmp_path):
    project = run_generator(tmp_path, include_mst=True)

    tabset = tree_tabset(project["panes"]["model"]["layout"])
    assert [tab["id"] for tab in tabset["children"]] == ["tree-1", "tree-2", "tree-3"]
    assert tabset["selected"] < len(tabset["children"])


def test_removes_mst_panel_when_input_is_absent(tmp_path):
    project = run_generator(tmp_path, include_mst=False)

    assert "mst0" not in project["files"]
    assert "tree-3" not in project["trees"]
    assert "tree-3" not in tab_ids(project["panes"]["model"]["layout"])


def test_viral_layout_keeps_two_tabs_in_one_pane(tmp_path):
    project = run_generator(tmp_path, include_mst=False)

    tabset = tree_tabset(project["panes"]["model"]["layout"])
    assert [tab["id"] for tab in tabset["children"]] == ["tree-1", "tree-2"]
    assert tabset["selected"] < len(tabset["children"])
    assert all(tabset["children"] for tabset in tabsets(project["panes"]["model"]["layout"]))
