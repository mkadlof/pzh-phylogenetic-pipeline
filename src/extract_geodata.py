import pandas as pd
import requests
import time
import click

def query_nominatim_API(query_class: str, query_value: str) -> str:
    time.sleep(2)
    params = {
        query_class: query_value,
        'format': 'jsonv2'
    }

    headers = {
        'User-Agent': 'MyGeoApp/1.0 (mlazniewski@gmail.com)'  # Replace with your real info
    }

    url = "https://nominatim.openstreetmap.org/search"
    response = requests.get(url, params=params, headers=headers)

    if response.status_code == 200 and response.json():
        result = response.json()[0]
        return f'{query_class}\t{query_value}\t{result["lat"]}\t{result["lon"]}\n'
    else:
        return f'{query_class}\t{query_value}\t0\t0\n'

@click.command()
@click.option('--input_metadata', help='[INPUT] Tab-separated file with metadata including country/state/city columns.',
              type=click.Path(exists=True), required=True)
@click.option('--output_metadata', help='[OUTPUT] Output file name',
              type=str, required=True)
@click.option('--features', help='[INPUT] Geographic features to query (e.g. country, region, city). Use multiple times.',
              type=click.Choice(['country', 'region']), multiple=True, required=True)
def main(input_metadata: str, output_metadata: str, features: tuple) -> None:
    """
    Main function for extracting geographical data
    """
    visited_feature = {feature: set() for feature in features}

    metadata = pd.read_csv(input_metadata, sep='\t')

    with open(output_metadata, 'w') as output_f:
        for _, row in metadata.iterrows():
            for feature in features:
                if feature in row and pd.notnull(row[feature]):
                    value = row[feature]
                    if value not in visited_feature[feature]:
                        data = query_nominatim_API(query_class=feature, query_value=value)
                        output_f.write(data)
                        visited_feature[feature].add(value)

if __name__ == '__main__':
    main()
