#!/usr/bin/env python3

"""
Skrypt do wygenerowania projektu mikroreat na podstawie wzorcowego pliku
w ktorym wstawaimy blob-y z danymi w zalozonych miejscach i podmieniamy timerstamp i naze
projektu
"""

import base64
import json
from datetime import datetime, timezone

import click


def text_file_to_base64(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_bytes = f.read()
    encoded = base64.b64encode(file_bytes).decode("utf-8")
    return encoded


def read_metadata_columns(metadata_path: str) -> list:
    """
    Read the column names of the (tab-separated) metadata file, in file order.

    The template project's table only shows the columns listed in
    tables.table-1.columns, so this list must be rebuilt from the metadata
    file actually produced for this run rather than hardcoded: metadata
    columns vary by organism (viral vs bacterial), by which HierCC levels a
    given cgMLST scheme reports, and by whatever extra columns the user put
    in WGS2Phylo's --supplemental-file (e.g. age, gender). Any column not in
    this list is otherwise hidden by default in the Microreact table pane.

    :param metadata_path: Path to the tab-separated metadata file.
    :return: Column names in file order.
    """
    with open(metadata_path, "r", newline="") as f:
        header_line = f.readline().rstrip("\r\n")
    if not header_line:
        return []
    return header_line.split("\t")


def remove_tree_panel(project: dict, tree_id: str) -> None:
    """Remove a tree and its pane (or tab) from a template project."""
    project.get('trees', {}).pop(tree_id, None)

    def prune_layout(node):
        children = node.get('children')
        if not isinstance(children, list):
            return

        retained_children = []
        for child in children:
            prune_layout(child)
            child_tabs = child.get('children')
            is_empty_tabset = child.get('type') == 'tabset' and child_tabs == []
            is_target_tab = child.get('type') == 'tab' and child.get('id') == tree_id
            if not is_empty_tabset and not is_target_tab:
                retained_children.append(child)
        node['children'] = retained_children

        if node.get('type') == 'tabset' and 'selected' in node:
            # A stale index would leave the tabset rendering nothing.
            node['selected'] = min(node['selected'], max(len(retained_children) - 1, 0))

    panes = project.get('panes', {}).get('model', {}).get('layout', {})
    prune_layout(panes)


@click.command()
@click.option('--input_json', help='[INPUT] Default project')
@click.option('--classical_tree', help='[OUTPUT] Classical tree in nwk format')
@click.option('--rescaled_tree', help='[OUTPUT] Rescaled tree in nwk format')
@click.option('--mst_tree', help='[OUTPUT] Sample-level MST in nwk format')
@click.option('--metadata', help='[OUTPUT] Metadata file')
@click.option('--project_name', help='[OUTPUT] Name of the project')
@click.option('--output', help='[OUTPUT] Output file')
def main(input_json, classical_tree, rescaled_tree, mst_tree, metadata, project_name, output):
    default_project = json.load(open(input_json))

    classical_tree_blob = text_file_to_base64(file_path=classical_tree)
    rescaled_tree_blob = text_file_to_base64(file_path=rescaled_tree)
    metadata_blob = text_file_to_base64(file_path=metadata)

    timestamp = datetime.now(timezone.utc).isoformat(timespec='milliseconds').replace('+00:00', 'Z')

    #  nazwy "sib4", "a8f0", "ecd0" oraz "mst0" sa na sztywno ustawione w projekcie
    #  nie wiem jak mikroreact zareaguje gdy bedzie N plikow z tymi samymi wartosciami
    default_project['files']['sib4']['blob'] = f"data:application/octet-stream;base64, {metadata_blob}"
    default_project['files']['a8f0']['blob'] = f"data:application/octet-stream;base64, {classical_tree_blob}"
    default_project['files']['ecd0']['blob'] = f"data:application/octet-stream;base64, {rescaled_tree_blob}"
    if mst_tree:
        mst_tree_blob = text_file_to_base64(file_path=mst_tree)
        default_project['files']['mst0']['blob'] = f"data:application/octet-stream;base64, {mst_tree_blob}"
    else:
        # Viral projects use this script and template without an MST input.
        default_project['files'].pop('mst0', None)
        remove_tree_panel(default_project, 'tree-3')

    default_project['meta']['name'] = project_name
    default_project['meta']['timestamp'] = timestamp

    # Rebuild the metadata table's column list from the actual metadata file
    # instead of relying on the template's hardcoded list, so every column
    # that is actually in this run's metadata is visible by default (and
    # nothing else stale is displayed as a dead/empty column).
    for table in default_project.get('tables', {}).values():
        table['columns'] = [{"field": field, "fixed": False} for field in read_metadata_columns(metadata)]

    #  dump json into a file
    with open(output, 'w') as f:
        json.dump(default_project, f, indent=4)


if __name__ == "__main__":
    main()
