import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "bin"))

from WGS2Phylo import get_fastqc_stats, get_contaminations_bacteria,  \
    get_sequencing_summary_bacteria, \
    get_amr_bacteria, \
    get_serovar_bacteria, \
    get_mlst_cgmlst, \
    get_viral_obligatory_data, \
    get_influenza_antiviral_data, \
    get_viral_kraken2_data, \
    get_viral_freyja_data, \
    get_viral_classification_data, \
    get_viral_genome_data

GOLDENS = {
    'campylo_illumina.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 375448,
            'number_of_bases_forward': 55200000,
            'reads_median_quality_forward': 37.0,
            'reads_median_length_forward': 150.0,
            'number_of_reads_reverse': 375448,
            'number_of_bases_reverse': 55200000,
            'reads_median_quality_reverse': 37.0,
            'reads_median_length_reverse': 150.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'Campylobacter coli',
            'kraken2_species_secondary' : 'Campylobacter jejuni',
            'kraken2_species_main_value' : 72.7,
            'kraken2_species_secondary_value' : 2.04,

            'metaphlan_species_main': 'Campylobacter coli',
            'metaphlan_species_secondary': 'None',
            'metaphlan_species_main_value': 100.0,
            'metaphlan_species_secondary_value': 0,

            'kmerfinder_species_main': 'Campylobacter coli',
            'kmerfinder_species_secondary': 'Campylobacter jejuni',
            'kmerfinder_species_main_value': 80.78,
            'kmerfinder_species_secondary_value': 16.59,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  1744267,
            'contigs_number': 33,
            'average_coverage': 57.53,
            'N50' : 232219,
            'L50' : 3,
            'Ns_value' : 257

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Ciprofloxacin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.T86I' } ,
                                   'Clindamycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Erythromycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Chloramphenicol': {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'gen',
                                                     'czynnik_nazwa': 'tet(O/32/O)',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'gen',
                                                     'czynnik_nazwa': "aph(2'')-Ib",
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.T86I' }
        },
        "test_get_serovar_bacteria" : 'Campylobacter_jejuni',
        "test_get_mlst_cgmlst" : {'mlst_id' : "828",
                                  'mlst_public' : "828",
                                  'cgmlst_id' : "local_41",
                                  'cgmlst_public' : "23885",
                                  'hc0': "Unknown",
                                  'hc2': "Unknown",
                                  'hc5': "local_41",
                                  'hc10' : "local_41",
                                  'hc20': "Unknown",
                                  "reasons" : []}
    },
    'campylo_nanopore.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': 57088,
            'number_of_bases': 159600000,
            'reads_median_quality': 38.0,
            'reads_median_length': 2000.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'Campylobacter coli',
            'kraken2_species_secondary' : 'Photobacterium leiognathi',
            'kraken2_species_main_value' : 81.29,
            'kraken2_species_secondary_value' : 11.05,


            'kmerfinder_species_main': 'Campylobacter coli',
            'kmerfinder_species_secondary': 'Campylobacter jejuni',
            'kmerfinder_species_main_value': 36.42,
            'kmerfinder_species_secondary_value': 10.7,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  1756767,
            'contigs_number': 3,
            'average_coverage': 36.91,
            'N50' : 1751062,
            'L50' : 1,
            'Ns_value' : 0

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Ciprofloxacin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' } ,
                                   'Clindamycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Erythromycin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': '23S',
                                                     'czynnik_mutacja': 'r.2075A>G' },
                                   'Chloramphenicol': {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'gen',
                                                     'czynnik_nazwa': 'tet(O)',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' }
        },
        "test_get_serovar_bacteria" : 'Campylobacter_jejuni',
        "test_get_mlst_cgmlst" : {'mlst_id' : "7125",
                                  'mlst_public' : "7125",
                                  'cgmlst_id' : "local_47",
                                  'cgmlst_public' : "44900",
                                  'hc0': "Unknown",
                                  'hc2': "Unknown",
                                  'hc5': "local_47",
                                  'hc10' : "local_47",
                                  'hc20': "Unknown",
                                  "reasons" : []}
    },
    'ecoli_illumina.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 832529,
            'number_of_bases_forward': 120100000,
            'reads_median_quality_forward': 37.0,
            'reads_median_length_forward': 150.0,
            'number_of_reads_reverse': 832529,
            'number_of_bases_reverse': 120200000,
            'reads_median_quality_reverse': 36.0,
            'reads_median_length_reverse': 150.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'Escherichia coli',
            'kraken2_species_secondary' : 'Escherichia albertii',
            'kraken2_species_main_value' : 11.73,
            'kraken2_species_secondary_value' : 0.14,

            'metaphlan_species_main': 'Escherichia coli',
            'metaphlan_species_secondary': 'None',
            'metaphlan_species_main_value': 100.0,
            'metaphlan_species_secondary_value': 0,

            'kmerfinder_species_main': 'Escherichia coli',
            'kmerfinder_species_secondary': 'None',
            'kmerfinder_species_main_value': 0.18,
            'kmerfinder_species_secondary_value': 0,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  5104809,
            'contigs_number': 238,
            'average_coverage': 44.98,
            'N50' : 96543,
            'L50' : 19,
            'Ns_value' : 241

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Ciprofloxacin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.S83A' } ,
                                   'Clindamycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Erythromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Chloramphenicol': {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.S83A' }
        },
        "test_get_serovar_bacteria" : 'O104:H4',
        "test_get_mlst_cgmlst" : {'mlst_id' : "678",
                                  'mlst_public' : "678",
                                  'cgmlst_id' : "231697",
                                  'cgmlst_public' : "231697",
                                  'hc0': "231697",
                                  'hc2': "231697",
                                  'hc5': "231697",
                                  'hc10' : "231697",
                                  'hc20': "1645",
                                  "reasons" : []}
    },
    'ecoli_nanopore.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': 208819,
            'number_of_bases': 498500000,
            'reads_median_quality': 41.0,
            'reads_median_length': 0.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'Escherichia coli',
            'kraken2_species_secondary' : 'Photobacterium leiognathi',
            'kraken2_species_main_value' : 55.3,
            'kraken2_species_secondary_value' : 8.76,


            'kmerfinder_species_main': 'Escherichia marmotae',
            'kmerfinder_species_secondary': 'Escherichia albertii',
            'kmerfinder_species_main_value': 18.49,
            'kmerfinder_species_secondary_value': 17.29,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  5267401,
            'contigs_number': 5,
            'average_coverage': 95.67,
            'N50' : 2017328,
            'L50' : 2,
            'Ns_value' : 0

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Ciprofloxacin' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.S83A' } ,
                                   'Clindamycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Erythromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Chloramphenicol': {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'mutacja_punktowa',
                                                     'czynnik_nazwa': 'gyrA',
                                                     'czynnik_mutacja': 'p.S83A' }
        },
        "test_get_serovar_bacteria" : 'O104:H4',
        "test_get_mlst_cgmlst" : {'mlst_id' : "678",
                                  'mlst_public' : "678",
                                  'cgmlst_id' : "local_22",
                                  'cgmlst_public' : "231697",
                                  'hc0': "local_22",
                                  'hc2': "231697",
                                  'hc5': "231697",
                                  'hc10' : "231697",
                                  'hc20': "1645",
                                  "reasons" : []}
    },
    'salmonella_illumina.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 1059384,
            'number_of_bases_forward': 155000000,
            'reads_median_quality_forward': 37.0,
            'reads_median_length_forward': 150.0,
            'number_of_reads_reverse': 1059384,
            'number_of_bases_reverse': 155000000,
            'reads_median_quality_reverse': 36.0,
            'reads_median_length_reverse': 150.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'Salmonella enterica',
            'kraken2_species_secondary' : 'Escherichia coli',
            'kraken2_species_main_value' : 5.4,
            'kraken2_species_secondary_value' : 0.22,

            'metaphlan_species_main': 'Salmonella enterica',
            'metaphlan_species_secondary': 'None',
            'metaphlan_species_main_value': 100.0,
            'metaphlan_species_secondary_value': 0,

            'kmerfinder_species_main': 'Salmonella enterica',
            'kmerfinder_species_secondary': 'Escherichia coli',
            'kmerfinder_species_main_value': 93.96,
            'kmerfinder_species_secondary_value': 5.4,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  5004068,
            'contigs_number': 78,
            'average_coverage': 55.41,
            'N50' : 225224,
            'L50' : 7,
            'Ns_value' : 840

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Ciprofloxacin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' } ,
                                   'Clindamycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Erythromycin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Chloramphenicol': {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'oporny',
                                                     'czynnik_typ': 'gen',
                                                     'czynnik_nazwa': 'tet(B)',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'wrazliw',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' }
        },
        "test_get_serovar_bacteria" : 'Typhimurium',
        "test_get_mlst_cgmlst" : {'mlst_id' : "34",
                                  'mlst_public' : "34",
                                  'cgmlst_id' : "local_241",
                                  'cgmlst_public' : "337784",
                                  'hc0': "local_241",
                                  'hc2': "337784",
                                  'hc5': "337784",
                                  'hc10' : "2",
                                  'hc20': "2",
                                  "reasons" : []}
    },
    'salmonella_nanopore.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': 39514,
            'number_of_bases': 88900000,
            'reads_median_quality': 40.0,
            'reads_median_length': 0.0,
        },
        "test_get_contaminations_bacteria" : {
            'kraken2_species_main' : 'brak',
            'kraken2_species_secondary' : 'brak',
            'kraken2_species_main_value' : -1,
            'kraken2_species_secondary_value' : -1,


            'kmerfinder_species_main': 'brak',
            'kmerfinder_species_secondary': 'brak',
            'kmerfinder_species_main_value': -1,
            'kmerfinder_species_secondary_value': -1,

        },
        "test_get_sequencing_summary_bacteria" : {
            'genome_length':  -1,
            'contigs_number': -1,
            'average_coverage': -1,
            'N50' : -1,
            'L50' : -1,
            'Ns_value' : -1

        },
        "test_get_amr_bacteria" : {'Azithromycin' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Ciprofloxacin' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' } ,
                                   'Clindamycin' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Erythromycin' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Chloramphenicol': {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Tetracycline' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Gentamicin' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' },
                                   'Nalidixic Acid' : {'opornos' : 'brak_danych',
                                                     'czynnik_typ': 'brak',
                                                     'czynnik_nazwa': 'brak',
                                                     'czynnik_mutacja': 'brak' }
        },
        "test_get_serovar_bacteria" : 'Unknown',
        "test_get_mlst_cgmlst" : {'mlst_id' : "Unknown",
                                  'mlst_public' : "Unknown",
                                  'cgmlst_id' : "Unknown",
                                  'cgmlst_public' : "Unknown",
                                  'hc0': "Unknown",
                                  'hc2': "Unknown",
                                  'hc5': "Unknown",
                                  'hc10' : "Unknown",
                                  'hc20': "Unknown",
                                  "reasons" : ['This module was eneterd with failed QC and poduced no valid output', 'This module was eneterd with failed QC and poduced no valid output']}

    },
    'influenza_illumina.json' :{
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 2866841,
            'number_of_bases_forward': 426800000,
            'reads_median_quality_forward': 33.0,
            'reads_median_length_forward': 150.0,
            'number_of_reads_reverse': 2866841,
            'number_of_bases_reverse': 426800000,
            'reads_median_quality_reverse': 33.0,
            'reads_median_length_reverse': 150.0,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025.INFL2.05",
            "virus": "influenza",
            "type" : "H1N1"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'S',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'S',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'S',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'S',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'S',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Alphainfluenzavirus influenzae',
            'kraken2_species_secondary': 'Homo sapiens',
            'kraken2_species_main_value': 96.45,
            'kraken2_species_secondary_value': 0.04,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': '5a.2a',
            'freyja_lineage_secondary': 'unk',
            'freyja_lineage_main_value': 1.0,
            'freyja_lineage_secondary_value': 0,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name_HA': '6B.1A.5a.2a',
            'nextclade_variant_name_NA': 'C.5.3.1'
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : 246.33,
            'total_length' : 13632,
            'number_of_Ns' : 200,
        }
    },
    'influenza_illumina_nodata.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': -1,
            'number_of_bases_forward': -1,
            'reads_median_quality_forward': -1,
            'reads_median_length_forward': -1,
            'number_of_reads_reverse': -1,
            'number_of_bases_reverse': -1,
            'reads_median_quality_reverse': -1,
            'reads_median_length_reverse': -1,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025.INFL2.13",
            "virus": "influenza",
            "type" : "Unknown"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Unknown',
            'kraken2_species_secondary': 'Unknown',
            'kraken2_species_main_value': -1,
            'kraken2_species_secondary_value': -1,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'Unknown',
            'freyja_lineage_secondary': 'Unknown',
            'freyja_lineage_main_value': -1,
            'freyja_lineage_secondary_value': -1,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name_HA' : 'Unknown',
            'nextclade_variant_name_NA' : 'Unknown'
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : -1,
            'total_length' : -1,
            'number_of_Ns' : -1,
        }


    },
    'influenza_nanopore.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': 149352,
            'number_of_bases': 154300000,
            'reads_median_quality': 28.0,
            'reads_median_length': 0.0,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025_INFL1_01",
            "virus": "influenza",
            "type" : "H3N2"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'S',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'R',
            'Zanamivir_mutation': 'A10G,A11T',
            'Peramivir_resistance_status': 'S',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'S',
            'Laninamivir_mutation': 'A10G',
            'Baloxavir_resistance_status': 'S',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Alphainfluenzavirus influenzae',
            'kraken2_species_secondary': 'Schaalia odontolytica',
            'kraken2_species_main_value': 95.94,
            'kraken2_species_secondary_value': 2.53,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': '2a.3a.1',
            'freyja_lineage_secondary': '3C.2a1b.2a',
            'freyja_lineage_main_value': 0.98,
            'freyja_lineage_secondary_value': 0.01,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name_HA' : '3C.2a1b.2a.2a.3a.1',
            'nextclade_variant_name_NA' : 'B.4.2'
        },
        'test_get_viral_genome_data' : {
            'average_coverage' :246.27,
            'total_length' : 13627,
            'number_of_Ns' : 198,
        }

    },
    'influenza_nanopore_nodata.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': -1,
            'number_of_bases': -1,
            'reads_median_quality': -1,
            'reads_median_length': -1,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025_INFL1_13",
            "virus": "influenza",
            "type" : "Unknown"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data": {
            'kraken2_species_main': 'Unknown',
            'kraken2_species_secondary': 'Unknown',
            'kraken2_species_main_value': -1,
            'kraken2_species_secondary_value': -1,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'Unknown',
            'freyja_lineage_secondary': 'Unknown',
            'freyja_lineage_main_value': -1,
            'freyja_lineage_secondary_value': -1,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name_HA' : 'Unknown',
            'nextclade_variant_name_NA' : 'Unknown'
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : -1,
            'total_length' : -1,
            'number_of_Ns' : -1,
        }

    },
    'rsv_illumina.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 978194,
            'number_of_bases_forward': 114500000,
            'reads_median_quality_forward': 34.0,
            'reads_median_length_forward': 2.0,
            'number_of_reads_reverse': 978194,
            'number_of_bases_reverse': 114700000,
            'reads_median_quality_reverse': 33.0,
            'reads_median_length_reverse': 2.0,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "SRR27383093_pass",
            "virus": "rsv",
            "type" : "A"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Orthopneumovirus hominis',
            'kraken2_species_secondary': 'Homo sapiens',
            'kraken2_species_main_value': 82.36,
            'kraken2_species_secondary_value': 1.85,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'A.D.5.2',
            'freyja_lineage_secondary': 'unk',
            'freyja_lineage_main_value': 1.0,
            'freyja_lineage_secondary_value': 0,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name' : 'A.D.5.2',
        },
'test_get_viral_genome_data' : {
            'average_coverage' : 245.26,
            'total_length' : 15233,
            'number_of_Ns' : 207,
        }
    },
    'rsv_illumina_nodata.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 1070780,
            'number_of_bases_forward': 120300000,
            'reads_median_quality_forward': 33.0,
            'reads_median_length_forward': 2.0,
            'number_of_reads_reverse': 1070780,
            'number_of_bases_reverse': 120300000,
            'reads_median_quality_reverse': 33.0,
            'reads_median_length_reverse': 2.0,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "SRR27383066_pass",
            "virus": "rsv",
            "type" : "unk"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Homo sapiens',
            'kraken2_species_secondary': 'Escherichia coli',
            'kraken2_species_main_value': 95.18,
            'kraken2_species_secondary_value': 0.03,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'Unknown',
            'freyja_lineage_secondary': 'Unknown',
            'freyja_lineage_main_value': -1,
            'freyja_lineage_secondary_value': -1,
        },
        'test_get_viral_classification_data' : {
            'nextclade_variant_name' : 'Unknown',
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : -1,
            'total_length' : -1,
            'number_of_Ns' : -1,
        }
    },
    'sars2_illumina.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': 80751,
            'number_of_bases_forward': 9900000,
            'reads_median_quality_forward': 39.0,
            'reads_median_length_forward': 150.0,
            'number_of_reads_reverse': 80751,
            'number_of_bases_reverse': 9900000,
            'reads_median_quality_reverse': 39.0,
            'reads_median_length_reverse': 150.0,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025_SARS2_11",
            "virus": "sars2",
            "type" : "sars2"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Severe acute respiratory syndrome-related coronavirus',
            'kraken2_species_secondary': 'Homo sapiens',
            'kraken2_species_main_value': 46.07,
            'kraken2_species_secondary_value': 40.32,
        },
        "test_get_viral_freyja_data": {
            'freyja_lineage_main': 'NH.3',
            'freyja_lineage_secondary': 'KP.3.1.1',
            'freyja_lineage_main_value': 0.01,
            'freyja_lineage_secondary_value': 0.01,
        },
        'test_get_viral_classification_data' : {
            'pangolin_variant_name' : 'Unassigned',
            'nextclade_variant_name' : '24A',
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : 22.09,
            'total_length' : 29871,
            'number_of_Ns' : 26313,
        }

    },
    'sars2_illumina_nodata.json': {
        'test_get_fastqc_stats': {
            'number_of_reads_forward': -1,
            'number_of_bases_forward': -1,
            'reads_median_quality_forward': -1,
            'reads_median_length_forward': -1,
            'number_of_reads_reverse': -1,
            'number_of_bases_reverse': -1,
            'reads_median_quality_reverse': -1,
            'reads_median_length_reverse': -1,
        },
        'test_get_viral_obligatory_data' : {
            "strain" : "ESIB_EQA_2025_SARS2_13",
            "virus": "sars2",
            "type" : "sars2"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Unknown',
            'kraken2_species_secondary': 'Unknown',
            'kraken2_species_main_value': -1,
            'kraken2_species_secondary_value': -1,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'Unknown',
            'freyja_lineage_secondary': 'Unknown',
            'freyja_lineage_main_value': -1,
            'freyja_lineage_secondary_value': -1,
        },
        'test_get_viral_classification_data': {
            'pangolin_variant_name': 'Unknown',
            'nextclade_variant_name': 'Unknown',
        },
        'test_get_viral_genome_data' : {
            'average_coverage' : -1,
            'total_length' : -1,
            'number_of_Ns' : -1,
        }
    },
    'sars2_nanopore.json': {
        'test_get_fastqc_stats': {
            'number_of_reads': 221267,
            'number_of_bases': 115700000,
            'reads_median_quality': 31,
            'reads_median_length': 400,
        },
        "test_get_viral_obligatory_data" : {
            "strain" : "ESIB_EQA_2025_SARS1_03",
            "virus": "sars2",
            "type" : "sars2"

        },
        "test_get_influenza_antiviral_data" : {
            'Oseltamivir_resistance_status': 'Unknown',
            'Oseltamivir_mutation': 'Unknown',
            'Zanamivir_resistance_status': 'Unknown',
            'Zanamivir_mutation': 'Unknown',
            'Peramivir_resistance_status': 'Unknown',
            'Peramivir_mutation': 'Unknown',
            'Laninamivir_resistance_status': 'Unknown',
            'Laninamivir_mutation': 'Unknown',
            'Baloxavir_resistance_status': 'Unknown',
            'Baloxavir_mutation': 'Unknown'
        },
        "test_get_viral_kraken2_data" : {
            'kraken2_species_main': 'Severe acute respiratory syndrome-related coronavirus',
            'kraken2_species_secondary': 'Serratia liquefaciens',
            'kraken2_species_main_value': 94.88,
            'kraken2_species_secondary_value': 0.21,
        },
        "test_get_viral_freyja_data" : {
            'freyja_lineage_main': 'KP.2.3.3',
            'freyja_lineage_secondary': 'unk',
            'freyja_lineage_main_value': 0.99,
            'freyja_lineage_secondary_value': 0,
        },
        'test_get_viral_classification_data': {
            'pangolin_variant_name': 'KP.2.3.3',
            'nextclade_variant_name': '24G',
        },
        'test_get_viral_genome_data' : {
            'average_coverage' :  225.18,
            'total_length' : 29845,
            'number_of_Ns' : 1482,
        }
    }

}

def load_json_content(json_path: Path):
    with open(json_path, 'r') as f:
        return json.load(f)

def test_get_fastqc_stats(json_file):
    content = load_json_content(json_file)

    result = {}
    if 'illumina' in str(json_file):
        result.update(get_fastqc_stats(content, "forward"))
        result.update(get_fastqc_stats(content, "reverse"))
    else:
        result.update(get_fastqc_stats(content, "forward"))

    expected = GOLDENS[json_file.name]['test_get_fastqc_stats']
    assert result == expected, f"{json_file.name}: FastQC stats output mismatch"

def test_get_contaminations_bacteria(bacterial_json):
    content = load_json_content(bacterial_json)
    result = get_contaminations_bacteria(content)

    expected = GOLDENS[bacterial_json.name]['test_get_contaminations_bacteria']

    assert result == expected, f"{bacterial_json.name}: Contamination analysis output mismatch"

def test_get_sequencing_summary_bacteria(bacterial_json):
    content = load_json_content(bacterial_json)
    result = get_sequencing_summary_bacteria(content)

    expected = GOLDENS[bacterial_json.name]['test_get_sequencing_summary_bacteria']

    assert result == expected, f"{bacterial_json.name}: Sequencing summary analysis output mismatch"

def test_get_amr_bacteria(bacterial_json):
    content = load_json_content(bacterial_json)
    result = get_amr_bacteria(content)

    expected = GOLDENS[bacterial_json.name]['test_get_amr_bacteria']

    assert result == expected, f"{bacterial_json.name}: Sequencing summary analysis output mismatch"


def test_get_serovar_bacteria(bacterial_json):
    content = load_json_content(bacterial_json)
    result = get_serovar_bacteria(content)['serovar']

    expected = GOLDENS[bacterial_json.name]['test_get_serovar_bacteria']

    assert result == expected, f"{bacterial_json.name}: Serovar extraction match"


def test_get_mlst_cgmlst(bacterial_json):
    content = load_json_content(bacterial_json)
    result = get_mlst_cgmlst(content)

    expected = GOLDENS[bacterial_json.name]['test_get_mlst_cgmlst']

    assert result == expected, f"{bacterial_json.name}: chMLST/MLST mismatch"

def test_get_viral_obligatory_data(viral_json):
    content = load_json_content(viral_json)
    result = get_viral_obligatory_data(content)

    expected = GOLDENS[viral_json.name]['test_get_viral_obligatory_data']

    assert result == expected, f"{viral_json.name}: viral obligatory columns mismatch"

def test_get_viral_kraken2_data(viral_json):
    content = load_json_content(viral_json)
    result = get_influenza_antiviral_data(content)

    expected = GOLDENS[viral_json.name]['test_get_influenza_antiviral_data']

    assert result == expected, f"{viral_json.name}: antiviral data columns mismatch"

def test_get_viral_kraken2_data(viral_json):
    content = load_json_content(viral_json)
    result = get_viral_kraken2_data(content)

    expected = GOLDENS[viral_json.name]['test_get_viral_kraken2_data']

    assert result == expected, f"{viral_json.name}: kraken2 data columns mismatch"

def test_get_viral_freyja_data(viral_json):
    content = load_json_content(viral_json)
    result = get_viral_freyja_data(content)

    expected = GOLDENS[viral_json.name]['test_get_viral_freyja_data']

    assert result == expected, f"{viral_json.name}: freyja data columns mismatch"


def test_get_viral_classification_data(viral_json):
    content = load_json_content(viral_json)
    result = get_viral_classification_data(content)

    expected = GOLDENS[viral_json.name]['test_get_viral_classification_data']

    assert result == expected, f"{viral_json.name}: viral classification data columns mismatch"

def test_get_viral_genome_data(viral_json):
    content = load_json_content(viral_json)
    result = get_viral_genome_data(content)

    expected = GOLDENS[viral_json.name]['test_get_viral_genome_data']

    assert result == expected, f"{viral_json.name}: viral genome data columns mismatch"
