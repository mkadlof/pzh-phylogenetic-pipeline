"""
Skrypt do wygenerowania projektu mikroreat na podstawie wzorcowego pliku
w ktorym wstawaimy blob-y z danymi w zalozonych miejscach i podmieniamy timerstamp i naze
projektu
"""


import click
import json
import base64
from datetime import datetime, timezone


def text_file_to_base64(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_bytes = f.read()
    encoded = base64.b64encode(file_bytes).decode("utf-8")
    return encoded

@click.command()
@click.option('--input_json', help='[INPUT] Default project')
@click.option('--classical_tree', help='[OUTPUT] Classical tree in nwk format')
@click.option('--rescaled_tree', help='[OUTPUT] Rescaled tree in nwk format')
@click.option('--metadata', help='[OUTPUT] Metadata file')
@click.option('--project_name', help='[OUTPUT] Name of the project')
@click.option('--output', help='[OUTPUT] Output file')
def main(input_json, classical_tree, rescaled_tree, metadata, project_name, output):
    default_project = json.load(open(input_json))

    classical_tree_blob = text_file_to_base64(file_path=classical_tree)
    rescaled_tree_blob = text_file_to_base64(file_path=rescaled_tree)
    metadata_blob = text_file_to_base64(file_path=metadata)

    timestamp = datetime.now(timezone.utc).isoformat(timespec='milliseconds').replace('+00:00', 'Z')

    #  nazwy "sib4" "a8f0", oraz "ecd0" sa na sztywno ustawione w projekcie
    #  nie wiem jak mikroreact zareaguje gdy bedzie N plikow z tymi samymi wartosciami
    default_project['files']['sib4']['blob'] = f"data:application/octet-stream;base64, {metadata_blob}"
    default_project['files']['a8f0']['blob'] = f"data:application/octet-stream;base64, {classical_tree_blob}"
    default_project['files']['ecd0']['blob'] = f"data:application/octet-stream;base64, {rescaled_tree_blob}"

    default_project['meta']['name'] = project_name
    default_project['meta']['timestamp'] = timestamp

    #  dump json into a file
    with open(output, 'w') as f:
        json.dump(default_project, f, indent=4)

if __name__ == "__main__":
    main()