import click
from pathlib import Path
from Bio import SeqIO

@click.command()
@click.option('--input_dir', type=click.Path(exists=True, file_okay=False))
@click.option('--segment_name')
def extract_segment_sequences(input_dir, segment_name):
    """
    Extract sequences containing SEGMENT_NAME from all .fasta files in INPUT_DIR.
    Saves matching sequences to SEGMENT_NAME.fasta.
    """
    input_dir = Path(input_dir)
    output_file = f"{segment_name}.fasta"
    matching_records = []

    for fasta_file in input_dir.glob("*.fasta"):
        for record in SeqIO.parse(fasta_file, "fasta"):
            if segment_name in record.description:
                matching_records.append(record)

    if matching_records:
        SeqIO.write(matching_records, output_file, "fasta")
        print(f"Saved {len(matching_records)} matching sequences to '{output_file}'")
    else:
        print("No matching sequences found.")

if __name__ == '__main__':
    extract_segment_sequences()

