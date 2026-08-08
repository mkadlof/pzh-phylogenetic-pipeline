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


def remove_tree_panel(project: dict, tree_id: str) -> None:
    """Remove a tree and its pane from a template project."""
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

    #  dump json into a file
    with open(output, 'w') as f:
        json.dump(default_project, f, indent=4)


if __name__ == "__main__":
    main()
