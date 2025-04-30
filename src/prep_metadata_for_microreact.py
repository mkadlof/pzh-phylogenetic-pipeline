import click
import pandas as pd

def read_metadata(plik:str) -> pd.DataFrame:
    df = pd.read_csv(plik, sep ="\t")
    return df

def read_coordintes(plik:str, feature_name:str) -> pd.DataFrame:
    """
    Function to read file with coordinates
    :param plik:
    :return:
    """
    dane = pd.read_csv(plik, sep='\t', engine='python', header=None)
    dane.columns = ['type', f'{feature_name}', 'latitude', 'longitude']
    dane = dane[dane.type == f'{feature_name}']
    dane = dane.drop(columns='type')
    return dane

def substitute_date(df:pd.DataFrame) -> pd.DataFrame:
    df['date'] = pd.to_datetime(df['date'], errors='coerce')
    df['year'] = df['date'].dt.year
    df['month'] = df['date'].dt.month
    df['day'] = df['date'].dt.day
    return df

@click.command()
@click.option('--metadata', type=click.Path(exists=True), help='[INPUT] File with metadata file', required=True)
@click.option('--coordinates', type=click.Path(exists=True), help='[INPUT] File with coordinates', required=True)
@click.option('--output', type=str, help='[OUTPUT] Output file')
@click.option('--level', help='[INPUT] Feature level can be country or city',  required=True,
              type=click.Choice(['country', 'city']))
def main(metadata:str, coordinates:str, output:str, level:str):
    metadata = read_metadata(metadata)
    if level not in metadata.columns:
        raise Exception(f'Invalid level: No {level} in metadata file')

    coordinates = read_coordintes(plik = coordinates,
                                  feature_name = level)
    metadata_new = metadata.merge(coordinates, on = level, how ='left')
    metadata_new[['latitude', 'longitude']] = metadata_new[['latitude', 'longitude']].fillna(0)
    metadata_new = substitute_date(metadata_new)
    metadata_new.to_csv(output, index=False, sep = "\t")
    # split date column into
    return True

if __name__ == "__main__":
    main()