import hashlib
import click

def generate_color_code(feature: str, value: str) -> str:
    """Generate a hex color code based on a hash of feature and value."""
    hash_obj = hashlib.md5(f"{feature}_{value}".encode())
    return f"#{hash_obj.hexdigest()[:6]}"

def generate_colors_from_file(input_file: str, output_file: str) -> None:
    """
    Read two-column file and generate unique hex color for each (feature, value) pair,
    and write to output file.
    """
    seen = set()

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue  # Skip malformed lines

            feature, value = parts
            key = (feature, value)

            if key not in seen:
                color = generate_color_code(feature, value)
                outfile.write(f"{feature}\t{value}\t{color}\n")
                seen.add(key)

@click.command()
@click.option('--input_file', '-i', type=click.Path(exists=True), required=True,
              help='Input TSV file with two columns. First should be feature, second should be value.')
@click.option('--output_file', '-o', type=click.Path(), required=True,
              help='Output file to write colours for each feature/value combination.')
def main(input_file, output_file):
    """Generate a file with unique hex colors for each feature-value pair."""
    generate_colors_from_file(input_file, output_file)

if __name__ == '__main__':
    main()
