from Bio import SeqIO
import click
import sys
from collections import Counter
from typing import Dict, Tuple, Any, List
import logging
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def extract_region_from_alignment(fasta_file:str,
                                  start:int,
                                  end:int,
                                  gene_name,
                                  max_gap:int,  gene_dict:Dict)-> List[Tuple[Dict[str, str], Dict[str, List]]]:
    """
    Extracts the region between positions start and end from all sequences in a FASTA alignment.

    :param fasta_file: Path to the input FASTA file.
    :param start: Start position (1-based index).
    :param end: End position (1-based index, inclusive).
    :param gene_dict: Dictionary holding updated alignment without constant sites
    """
    gene_temporal_dict = {}
    gene_constant_sites = {}
    with open(fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            extracted_seq = record.seq[start:end]
            gene_temporal_dict[record.id] = str(extracted_seq).lower()

    # If any sequence for a given gene contains too many gaps,
    # exclude that gene from the analysis by returning an empty dict with the same keys.
    for key, value in gene_temporal_dict.items():
        if value.count("-") / len(value) > max_gap:
            return {key: "" for key in gene_temporal_dict}, {gene_name:[0, 0, 0, 0]}

    #  remove constant sites
    # gene_temporal_dict = remove_constant_sites(gene_temporal_dict)
    # to calculate values for ASC_STAM correction for a gene
    gene_constant_sites[gene_name] = count_constant_for_dict(gene_temporal_dict)
    # remove non-SNPs from the alignment
    gene_temporal_dict = remove_constant_sites_amd_indels(gene_temporal_dict)
    gene_dict[gene_name] = gene_temporal_dict
    return  gene_dict, gene_constant_sites


def remove_constant_sites(slownik_alignmentu:Dict) -> Dict:
    """
    Funkcja usuwa z alignmentu te pozycje w ktorym sa albo te same nuklotydy + ewenutalne N-ki
    :param slownik_alignmentu: slownik z alignmentem
    :return:
    """
    position_counts = [Counter(bases) for bases in zip(*list(slownik_alignmentu.values()))]
    positions_to_keep = [i for i, counter in enumerate(position_counts) if ((len(counter.keys()) > 2) or(len(counter.keys()) == 2 and  "n" not in counter.keys()))]
    slownik_alignmentu_out = {key: ''.join(val[i] for i in positions_to_keep) for key, val in slownik_alignmentu.items()}

    return slownik_alignmentu_out

def count_constant_sites(plik_alignment:str) -> List[Any]:
    """
    Calculates number of constant sites in an alignment
    Returns a list of 4 ints number of constant sites with A, T, G and C
    :param slownik_alignmentu: slownik z alignmentem
    :return:
    """
    A, T, G, C = 0, 0, 0, 0
    temporal_dict = {}
    with open(plik_alignment, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            temporal_dict[record.id] = str(record.seq).upper()

    position_counts = [Counter(bases) for bases in zip(*list(temporal_dict.values()))]

    for counter in position_counts:
        if len(counter.keys()) == 1 and list(counter.keys())[0] == "A":
            A += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "T":
            T += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "G":
            G += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "C":
            C += 1


    return [A, C, G, T]


def count_constant_for_dict(slownik_alignmentu: Dict) -> List[Any]:
    """
    Calculates number of constant sites in an alignment
    Returns a list of 4 ints number of constant sites with A, T, G and C
    :param slownik_alignmentu: slownik z alignmentem
    :return:
    """
    A, T, G, C = 0, 0, 0, 0

    slownik_alignmentu_tmp = {}
    for key, value in slownik_alignmentu.items():
        slownik_alignmentu_tmp[key] = value.upper()

    position_counts = [Counter(bases) for bases in zip(*list(slownik_alignmentu_tmp.values()))]

    for counter in position_counts:
        if len(counter.keys()) == 1 and list(counter.keys())[0] == "A":
            A += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "T":
            T += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "G":
            G += 1
        elif len(counter.keys()) == 1 and list(counter.keys())[0] == "C":
            C += 1

    return [A, C, G, T]

def remove_constant_sites_amd_indels(slownik_alignmentu:Dict) -> Dict:
    """
    Funkcja usuwa z alignmentu te pozycje w ktorym sa albo te same nuklotydy + ewenutalne N-ki
    :param slownik_alignmentu: slownik z alignmentem
    :return:
    """
    position_counts = [Counter(bases) for bases in zip(*list(slownik_alignmentu.values()))]
    positions_to_keep = [i for i, counter in enumerate(position_counts) if ((len(counter.keys()) > 2 and "-" not in counter.keys())
                                                                            or (len(counter.keys()) == 2 and "n" not in counter.keys()) and "-" not in counter.keys())]
    slownik_alignmentu_out = {key: ''.join(val[i] for i in positions_to_keep) for key, val in slownik_alignmentu.items()}

    return slownik_alignmentu_out

@click.command()
@click.option('--input_fasta', help='[INPUT] a fasta of concatated core genes sequences from Roary',
              type=click.Path(), required=True)
@click.option('--input_fasta_annotation', help='[INPUT] a file with core genome annotation in embl format ',
              type=click.Path(), required=True)
@click.option('--model', help='[INPUT] a model that will be used for a phylogenetic analysis',
              type=str, required=True, default="GTR+G+FO")
@click.option('--max_gap', help='[INPUT] Max gap length observed in any sample with respect to the alignment '
                                'length to keep a gene for downstream phylogenetic analysis',
              type=int, required=True, default=50)
@click.option('--output_fasta', help='[OUTPUT] a fasta of concatated core genes sequences without constant sites',
              type=str, required=True)
@click.option('--output_partition', help='[OUTPUT] a partition file for raxml program',
              type=str, required=True)
@click.option('--cpus', help='[INPUT] a number of cpus to use',
              type=int, required=False, default=1)

def main(input_fasta, input_fasta_annotation, model, output_fasta, max_gap, output_partition, cpus) -> None:
    alignment_dict = {}
    ASAM_correction = {}
    pool = Pool(cpus)
    jobs = []
    gen_list = []

    # Count constant sites

    A_num, T_num, G_mum, C_num = count_constant_sites(plik_alignment=input_fasta)

    with open(input_fasta_annotation, "r") as input_annotation_handle:
        record = SeqIO.read(input_annotation_handle, "embl")  # Read the first record
        for feature in record.features:
            start = int(feature.location.start)
            end = int(feature.location.end)
            gene_name = feature.qualifiers["label"][0]
            gen_list.append(gene_name)
            #  print(f'Sending to pool for gen {gene_name}, {start}, {end}')
            jobs.append(pool.apply_async(extract_region_from_alignment, (input_fasta, start, end, gene_name, max_gap, alignment_dict)))

    # wykonujemy zadania asynchronicznie

    for job in jobs:
        alignment_dict = {**alignment_dict, **job.get()[0]}
        ASAM_correction = {**ASAM_correction, **job.get()[1]}
    pool.close()

    alignment_dict_final = {}
    old_end = 0
    partition_dict = {}

    # let the final alignment and partition file retain order from roary output
    for gen in gen_list:
        for sample_id, sekwencja in alignment_dict[gen].items():
            if sample_id not in alignment_dict_final.keys():
                alignment_dict_final[sample_id] = sekwencja
            else:
                alignment_dict_final[sample_id] += sekwencja
        #  add to partition data inf regarding genes that after
        #  removal of constant sites still have non-0 length
        if len(sekwencja) > 0:
            partition_dict[gen] = [old_end+1, old_end + len(sekwencja)]
            old_end = old_end + len(sekwencja)

        else:
            logging.info(f"{gen} has no variable positions among analyzed samples")


    #  dump new alignment to fasta file
    with open(output_fasta, "w") as output_handle:
        for klucz,wartosc in alignment_dict_final.items():
            output_handle.write(f">{klucz}\n{wartosc}\n")

    #  dump new partition into a file
    print(ASAM_correction)
    with open(output_partition, "w") as output_handle:
        for gen, partition in partition_dict.items():
            STAM_correction = f"ASC_STAM{{{'/'.join(map(str, ASAM_correction[gen]))}}}"
            output_handle.write(f"{model}+{STAM_correction}, {gen}={partition[0]}-{partition[1]}\n")

    with open("constant_sites.txt", "w") as output_handle:
        output_handle.write(f"{A_num}\t{T_num}\t{G_mum}\t{C_num}\n")

if __name__ == '__main__':
    main()


