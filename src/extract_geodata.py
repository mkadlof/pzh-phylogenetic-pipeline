import pandas as pd
import requests
import time
import click
from typing import Dict

def query_nominatim_API(query: Dict[str, str], query_key:str) -> str:
    """ query is a dict with key values, where key is understto as a valid class
    by the API  search engin https://nominatim.org/release-docs/latest/api/Search/
    query_key is the name of the feature we want actually extract e.g.
    query might be {'country': 'France', "city": "Paris"} and
    query_key is "city" as we want geodata for Paris, France
    """
    time.sleep(2)

    params = {
        **query,
        'format': 'jsonv2'
    }

    headers = {
        'User-Agent': 'MyGeoApp/1.0 (mlazniewski@gmail.com)'  # Replace with your real info
    }

    url = "https://nominatim.openstreetmap.org/search"
    response = requests.get(url, params=params, headers=headers)

    if response.status_code == 200 and response.json():
        result = response.json()[0]
        return f'{query_key}\t{query[query_key]}\t{result["lat"]}\t{result["lon"]}\n'
    else:
        return f'{query_key}\t{query[query_key]}\t0\t0\n'

@click.command()
@click.option('--input_metadata', help='[INPUT] Tab-separated file with metadata including country/state/city columns.',
              type=click.Path(exists=True), required=True)
@click.option('--output_metadata', help='[OUTPUT] Output file name',
              type=str, required=True)
@click.option('--features', help='[INPUT] Geographic features to query (e.g. country, region, city). Use multiple times.',
              type=click.Choice(['country', 'state', 'city']), multiple=True, required=True)
def main(input_metadata: str, output_metadata: str, features: tuple) -> None:
    """
    Main function for extracting geographical data
    """
    visited_feature = {feature: set() for feature in features}

    metadata = pd.read_csv(input_metadata, sep='\t')

    with open(output_metadata, 'w') as output_f:
        for _, row in metadata.iterrows():
            if 'country' in row and 'country' in features:
                feature = 'country'
                value = row[feature]
                if value not in visited_feature[feature]:
                    data = query_nominatim_API(query={'country': row['country']},
                                               query_key=feature)
                    output_f.write(data)
                    visited_feature[feature].add(value)
            if 'state' in row and 'state' in features:
                feature = 'state'
                value = row[feature]
                if value not in visited_feature[feature]:
                    data = query_nominatim_API(query={'state': row['state']},
                                               query_key=feature)
                    output_f.write(data)
                    visited_feature[feature].add(value)
            if 'city' in row and 'city' in features:
                feature = 'city'
                value = row[feature]
                if value not in visited_feature[feature]:
                    data = query_nominatim_API(query={'country': row['country'], 'city': row['city']},
                                               query_key=feature)
                    output_f.write(data)
                    visited_feature[feature].add(value)


if __name__ == '__main__':
    main()
