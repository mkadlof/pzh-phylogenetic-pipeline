import os
import json
import csv
import click
from typing import Dict
import glob
from pathlib import Path
import gzip
import shutil

### Helper functions to extract specific info from JSON files ###

def get_viral_obligatory_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {})

    strain = output_data.get("sampleId", '') or "Unknown"
    virus =  output_data.get("pathogen", '') or "Unknown"

    type = "Unknwon"

    if virus == 'influenza':
        type = output_data.get("infl_data", {}).get('subtype_name', '') or "Unknown"
    elif virus == 'sars2':
        type = virus or "Unknown"
    elif virus == 'rsv':
        type = output_data.get("rsv_data", {}).get('type', '') or "Unknown"

    return {"strain": strain, "virus": virus, "type": type}


def get_influenza_antiviral_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {}).get('infl_data', {}).get('resistance_data', [])
    antivirals = ['Oseltamivir', 'Zanamivir', 'Peramivir', 'Laninamivir', 'Baloxavir']
    output_dict = {}
    for av in antivirals:
        output_dict[f"{av}_resistance_status"] = "Unknown"
        output_dict[f"{av}_mutation"] = "Unknown"

    for drug_data in output_data:
        drug_name = drug_data.get("drug_name") or 'Unknown'
        drug_resistance_status = drug_data.get('drug_resistance_status', '') or "Unknown"
        drug_resistance_mutations = ",".join(drug_data.get('mutation_list_data_reference_numbering')) or 'Unknown'
        if drug_name in antivirals:
            output_dict[f"{drug_name}_resistance_status"] =drug_resistance_status
            output_dict[f"{drug_name}_mutation"] = drug_resistance_mutations

    return output_dict

def get_viral_kraken2_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {}).get('contamination_data', [])
    output_dict = {
        'kraken2_species_main': 'Unknown',
        'kraken2_species_secondary': 'Unknown',
        'kraken2_species_main_value': -1,
        'kraken2_species_secondary_value': -1,

    }

    for program_data in output_data:
        program_name = program_data.get("program_name") or 'Unknown'
        status = program_data.get("status") or 'Unknown'
        if status == 'tak' and program_name == 'kraken2':
            output_dict = {
                'kraken2_species_main': program_data.get('main_species_name', '') or 'Unknown',
                'kraken2_species_secondary': program_data.get('secondary_species_name', '') or 'Unknown',
                'kraken2_species_main_value': program_data.get('main_species_value', ''),
                'kraken2_species_secondary_value': program_data.get('secondary_species_value', '')

            }

    return output_dict


def get_viral_freyja_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {}).get('freyja_data', {})
    output_dict = {
        'freyja_lineage_main': 'Unknown',
        'freyja_lineage_secondary': 'Unknown',
        'freyja_lineage_main_value': -1,
        'freyja_lineage_secondary_value': -1,

    }

    status = output_data.get("status") or 'Unknown'
    if status == 'tak' :
        output_dict = {
            'freyja_lineage_main': output_data.get('freyja_lineage1_name', '') or 'Unknown',
            'freyja_lineage_secondary': output_data.get('freyja_lineage2_name', '') or 'Unknown',
            'freyja_lineage_main_value': output_data.get('freyja_lineage1_abundance', ''),
            'freyja_lineage_secondary_value': output_data.get('freyja_lineage2_abundance', '') ,

        }

    return output_dict

def get_viral_classification_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {}).get('viral_classification_data', [])
    virus = json_data.get("output", {}).get("pathogen", '') or "Unknown"

    output_dict = {}

    # classification pango and/or nextclade species dependent
    if virus == 'sars2':
        output_dict['nextclade_variant_name'] = 'Unknown'
        output_dict['pangolin_variant_name'] = 'Unknown'
        for program_out in output_data:
            if program_out.get('status') == 'tak':
                if program_out.get('program_name') == 'scorpion':
                    output_dict['pangolin_variant_name'] = program_out.get('variant_name')
                elif program_out.get('program_name') == 'Nextclade':
                    output_dict['nextclade_variant_name'] = program_out.get('variant_name')
    elif virus == 'influenza':
        output_dict['nextclade_variant_name_HA'] = 'Unknown'
        output_dict['nextclade_variant_name_NA'] = 'Unknown'
        for program_out in output_data:
            if program_out.get('status')  == 'tak':
                if program_out.get('program_name')  == 'Nextclade' and program_out.get('sequence_source') == 'HA':
                    output_dict['nextclade_variant_name_HA'] = program_out.get('variant_name')
                elif program_out.get('program_name')  == 'Nextclade' and program_out.get('sequence_source') == 'NA':
                    output_dict['nextclade_variant_name_NA'] = program_out.get('variant_name')

    elif virus == 'rsv':
        output_dict['nextclade_variant_name'] = 'Unknown'
        for program_out in output_data:
            if program_out.get('status') == 'tak' and program_out.get('program_name') == 'Nextclade':
                    output_dict['nextclade_variant_name'] = program_out.get('variant_name')

    return output_dict

def get_viral_genome_data(json_data: Dict) -> Dict:
    output_data = json_data.get("output", {}).get('viral_genome_data', {})

    output_dict = {}
    output_dict['average_coverage'] = -1
    output_dict['total_length'] = -1
    output_dict['number_of_Ns'] = -1

    if output_data.get("status", '') == 'tak':
        output_dict['average_coverage'] = output_data.get('average_coverage_value', -1)
        output_dict['total_length'] = output_data.get('total_length_value', -1)
        output_dict['number_of_Ns'] = output_data.get('number_of_Ns_value', -1)

    return output_dict

# HierCC clustering levels reported in the output metadata, in normal mode
# (i.e. not gated behind --extra-fields). Not every organism/scheme emits
# every level (e.g. Campylobacter has no HC0/HC2/HC20), so missing levels
# fall back to 'Unknown' rather than raising an error.
HIERCC_LEVELS = ["0", "2", "5", "10", "20"]


def get_mlst_cgmlst(json_data: Dict) -> Dict:
    """
    Extract MLST and cgMLST profile IDs, plus HierCC clustering group IDs for
    each level in HIERCC_LEVELS.

    Returns a dict with mlst_id, cgmlst_id, mlst_public, cgmlst_public,
    hc<level> for each level in HIERCC_LEVELS, and reasons.
    """

    output_data = json_data.get("output", {})
    mlst_entries = output_data.get("mlst_data", [])

    reasons = []

    hiercc_by_level = {level: 'Unknown' for level in HIERCC_LEVELS}
    cgmlst_id = 'Unknown'
    cgmlst_public = "Unknown"
    mlst_id = 'Unknown'
    mlst_public = "Unknown"

    for entry in mlst_entries:
        scheme_name = str(entry.get("scheme_name", "")).lower()

        if scheme_name.startswith("cgmlst"):
            if entry.get('status').lower() != 'tak':
                reasons.append(entry.get('error_message'))
            else:
                cgmlst_id = entry.get('profile_id')
                cgmlst_public = entry.get('closest_external_profile_id')

                hiercc = entry.get('hiercc_clustering_internal_data', [])
                for level in hiercc:
                    level_name = str(level.get('level', ''))
                    if level_name in hiercc_by_level:
                        hiercc_by_level[level_name] = level.get('group_id')


        elif scheme_name.startswith("mlst"):
            if entry.get('status').lower() != 'tak':
                reasons.append(entry.get('error_message'))
            else:
                mlst_id = entry.get('profile_id')
                mlst_public = entry.get('closest_external_profile_id')


    result = {'mlst_id' : mlst_id,
              'cgmlst_id' : cgmlst_id,
              'mlst_public' : mlst_public,
              'cgmlst_public' : cgmlst_public,
              'reasons': reasons}
    for level in HIERCC_LEVELS:
        result[f'hc{level}'] = hiercc_by_level[level]
    return result


def get_serovar_bacteria(json_data:Dict) -> Dict:
    """
    Extract serovar (species/serotype) information based on genus.
    Supports Campylobacter, Escherichia, and Salmonella.
    Returns (serovar, reasons) where reasons is a list of warnings/errors.
    """
    serovar = 'Unknown'
    reason = []

    output_data = json_data.get("output", {})

    # Extract species from json to determine output of which program (ectype, seqser etc. to access)
    genus = output_data.get("pathogen_predicted_genus", "") or output_data.get("Genus", "")
    species = output_data.get("pathogen_predicted_species", "") or output_data.get("Species", "")

    if genus and genus.lower() == "campylobacter":
        if species:
            serovar = species
            if not serovar.lower().startswith(genus.lower()):
                serovar = f"{genus}_{species}"
        else:
            reason.append("Species not found for Campylobacter")

    elif genus and genus.lower() == "escherichia":
        ser = None
        for entry in output_data.get("antigenic_data", []):
            if entry.get("program_name") == "ectyper":
                if entry.get("status", "").lower() == "tak":
                    ser = entry.get("serotype_name")
                break
        if ser:
            serovar = ser
        else:
            reason.append("E. coli serotype (ECTyper) not found")

    elif genus and genus.lower() == "salmonella":
        ser = None
        for entry in output_data.get("antigenic_data", []):
            if entry.get("program_name") == "seqsero":
                if entry.get("status", "").lower() == "tak":
                    ser = entry.get("serotype_name")
                break
        if ser:
            serovar = ser
        else:
            reason.append("Salmonella serotype (SeqSero) not found")

    else:
        if genus:
            reason.append(f"Unsupported genus {genus}")
        else:
            reason.append("Genus not predicted")

    return {"serovar" : serovar,
            "reasons" : reason}

# Sections of bacterial related data
def get_sequencing_summary_bacteria(json_data):
    """Extract assembly summary stats (post-filtering) for bacterial genomes."""
    data_list = json_data.get("output", {}).get("bacterial_genome_data", [])
    for entry in data_list:
        if entry.get("step_name") == "post-filtering":
            # Return a dict of key assembly metrics (with safe defaults if missing)
            return {
                'genome_length': entry.get('total_length_value', -1),
                'contigs_number': entry.get('number_of_contigs_value', -1),
                'average_coverage': entry.get('average_coverage_value', -1),
                'N50': entry.get('N50_value', -1),
                'L50': entry.get('L50_value', -1),
                'Ns_value': entry.get('number_of_Ns_value', -1)
            }
    return {}

def get_amr_bacteria(json_data):
    """Extract antibiotic resistance results for key antibiotics (EQA panel)."""
    # List of antibiotics of interest (for EQA 2025 Campylobacter/Salmonella)
    KNOWN_ANTIBIOTICS = [
        'Azithromycin', 'Ciprofloxacin', 'Clindamycin', 'Erythromycin',
        'Chloramphenicol', 'Tetracycline', 'Gentamicin', 'Nalidixic Acid'
    ]
    output = {}
    # Initialize defaults for all known antibiotics: susceptible (wrażliwy) with no resistance factor
    # These defaults are kept if Resfinder was actually called
    for ab in KNOWN_ANTIBIOTICS:
        output[ab] = {
            'opornos': 'wrazliw',       # susceptible
            'czynnik_typ': 'brak',      # factor type
            'czynnik_nazwa': 'brak',    # factor name
            'czynnik_mutacja': 'brak'   # mutation detail
        }
    # Locate drug resistance data in JSON
    drug_data = json_data.get("output", {}).get("drug_resistance_data", [])
    for program in drug_data:
        if program.get("program_name") == "ResFinder/PointFinder" and program.get("status", "").lower() == "tak":
            for antibiotic in program.get("program_data", []):
                ab_name = antibiotic.get("antibiotic_name")
                if ab_name in output:  # only consider the known antibiotics
                    # Check if resistant
                    if antibiotic.get("antibiotic_status", "").lower() == "oporny":
                        # Get the first resistance factor if available (gene or mutation)
                        factor = {}
                        factors = antibiotic.get("antibiotic_resistance_data", [])
                        if factors:
                            factor = factors[0]
                        output[ab_name] = {
                            'opornos': 'oporny',  # resistant
                            'czynnik_typ': factor.get('factor_type_name', 'brak'),
                            'czynnik_nazwa': factor.get('factor_name', 'brak'),
                            'czynnik_mutacja': factor.get('factor_mutation', 'brak')
                        }
                    else:
                        # If not resistant (wrażliwy), ensure it stays marked as susceptible with no factor
                        output[ab_name] = {
                            'opornos': 'wrazliw',
                            'czynnik_typ': 'brak',
                            'czynnik_nazwa': 'brak',
                            'czynnik_mutacja': 'brak'
                        }
        elif program.get("program_name") == "ResFinder/PointFinder" and ( program.get("status", "").lower() == "nie" or program.get("status", "").lower() == "blad" ):
            # in case pipeline did not produce valid output for resfinde change defaults regarding resistance to "undetermined"
            for ab in KNOWN_ANTIBIOTICS:
                output[ab] = {
                    'opornos': 'brak_danych',  # poor sequencing no data
                    'czynnik_typ': 'brak',  # factor type
                    'czynnik_nazwa': 'brak',  # factor name
                    'czynnik_mutacja': 'brak'  # mutation detail
                }

    return output



def get_contaminations_bacteria(json_data):
    """Extract main and secondary species identification from contamination checks."""
    output = {}
    contamination_list = json_data.get("output", {}).get("contamination_data", [])
    for entry in contamination_list:
        program = entry.get('program_name')
        if not program:
            continue
        # Primary species and value (percentage or coverage)
        output[f"{program}_species_main"] = entry.get('main_species_name', 'brak')
        output[f"{program}_species_secondary"] = entry.get('secondary_species_name', 'brak')
        output[f"{program}_species_main_value"] = entry.get('main_species_value', entry.get('main_species_coverage', -1))
        output[f"{program}_species_secondary_value"] = entry.get('secondary_species_value', entry.get('secondary_species_coverage', -1))
    return output

def get_fastqc_stats(json_data, direction):
    """Extract FastQC read count and quality stats for forward (R1) or reverse (R2) reads."""
    output = {}
    seq_summary = json_data.get("output", {}).get("sequencing_summary_data", [])
    for entry in seq_summary:
        fname = entry.get('file_name', '')
        if direction == "forward" and ("R1" in fname or '_1.fastq' in fname):
            if entry.get("step_name", '') != 'post-filtering':
                # for viruses we have both pre- and post filtering results we keep pre-filtering for now
                # this key is not present in bacterial results so get should return ''
                output[f"number_of_reads_{direction}"] = entry.get('number_of_reads_value', -1)
                output[f"number_of_bases_{direction}"] = entry.get('number_of_bases_value', -1)
                output[f"reads_median_quality_{direction}"] = entry.get('reads_median_quality_value', -1)
                output[f"reads_median_length_{direction}"] = entry.get('reads_median_length_value', -1)
        elif direction == "reverse" and ("R2" in fname or '_2.fastq' in fname):
            if entry.get("step_name", '') != 'post-filtering':
                output[f"number_of_reads_{direction}"] = entry.get('number_of_reads_value', -1)
                output[f"number_of_bases_{direction}"] = entry.get('number_of_bases_value', -1)
                output[f"reads_median_quality_{direction}"] = entry.get('reads_median_quality_value', -1)
                output[f"reads_median_length_{direction}"] = entry.get('reads_median_length_value', -1)
        elif direction == "forward" and "R1" not in fname and "R2" not in fname and '_1.fastq' not in fname and '_2.fastq' not in fname:
            if entry.get("step_name", '') != 'post-filtering':
                # Placeholder gor nanopo re data that my lack R1/R2 in their filename
                output[f"number_of_reads"] = entry.get('number_of_reads_value', -1)
                output[f"number_of_bases"] = entry.get('number_of_bases_value', -1)
                output[f"reads_median_quality"] = entry.get('reads_median_quality_value', -1)
                output[f"reads_median_length"] = entry.get('reads_median_length_value', -1)
    return output


# Base metadata columns that are always populated by the pipeline itself.
# Supplemental-file columns using these names are ignored to avoid silently
# overwriting values the pipeline computes from the WGS JSON output.
RESERVED_METADATA_FIELDS = {"strain", "virus", "type", "serovar", "mlst", "cgmlst"} | {f"hc{level}" for level in HIERCC_LEVELS}


def sanitize_supplemental_value(value):
    """
    Sanitize a free-text value read from a user-supplied supplemental metadata file.

    The phylogenetic pipeline never filters extra metadata columns (they end up
    base64-blobbed into the Microreact project as-is), so arbitrary free-text
    fields such as age or gender are never executed as code anywhere in this
    codebase. The two concrete risks are therefore not code/shell injection but:
    (1) an embedded tab/newline/carriage-return silently shifting every column
    to its right in the output TSV, and (2) "CSV/formula injection" -- a value
    starting with '=', '+', '-' or '@' being interpreted as a formula by
    Excel/LibreOffice/Google Sheets if a human later opens the file there.

    :param value: Raw string value read from the supplemental file, or None.
    :return: A value safe to write as a single TSV cell; "N/A" if empty/missing.
    """
    if value is None:
        return "N/A"
    value = str(value).strip()
    if value == "":
        return "N/A"
    # Collapse characters that would otherwise break TSV column alignment.
    value = value.replace("\t", " ").replace("\r", " ").replace("\n", " ")
    # Drop other control characters defensively (keep normal whitespace/text).
    value = "".join(ch for ch in value if ch == " " or ch.isprintable()).strip()
    if value == "":
        return "N/A"
    # Defuse spreadsheet formula injection: prefix with a single quote, which
    # Excel/LibreOffice/Sheets treat as "force text" and simply strip on
    # display, instead of letting the leading character trigger a formula.
    if value[0] in ("=", "+", "-", "@"):
        value = "'" + value
    return value


#### CLI and main processing function ####

@click.command()
@click.option('--output_dir', required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory containing sample subdirectories with JSON/FASTA outputs.")
@click.option('--supplemental-file', type=click.Path(exists=True, readable=True), default=None,
              help="Optional TSV/CSV file with extra metadata (columns: date, region, country, division, city).")
@click.option('--id-column', default='id', show_default=True,
              help="Column name in the supplemental file that identifies the sample ID.")
@click.option('--output-prefix', default='metadata/metadata', show_default=True,
              help="Prefix for output files (aggregated will be <prefix>.tsv, per-sample will be <Sample>.<prefix>.tsv).")
@click.option('--extra-fields', is_flag=True, help="If set, include full WGS output fields (AMR, contamination, FastQC stats).")
@click.option('--organism', type=click.Choice(['sars-cov-2', 'influenza', 'rsv', 'salmonella', 'escherichia', 'campylobacter']), required=True,
              help="Name of the organism for which WGS sequencing was carried out.")
@click.option("--with-fasta/--without-fasta", default=True,
              help="Create fasta file directory used as an input for phylogentic pipline. Use --without-fasta to skip fasta creation")
def generate_metadata(output_dir, supplemental_file, id_column, output_prefix, extra_fields, organism, with_fasta):
    """
    Generate a phylogenetic metadata TSV from pipeline JSON outputs.
    """
    # Define optional metadata fields and read supplemental file if provided
    optional_fields = ["date", "region", "country", "division", "city"]
    optional_fields_present = []
    # Any other column the user puts in the supplemental file (e.g. age,
    # gender) is passed through to the output as-is, in file order, under its
    # original header text: it is not used anywhere in this pipeline besides
    # being base64-blobbed into the Microreact project, so it is not this
    # script's business to filter it. See sanitize_supplemental_value() for
    # the (non-code-execution) risks that are guarded against.
    extra_fields_present = []
    sup_data = {}
    if supplemental_file:
        # Determine delimiter from extension or content
        sup_ext = os.path.splitext(supplemental_file)[1].lower()
        if sup_ext in ('.tsv', '.tab'):
            delim = '\t'
        elif sup_ext == '.csv':
            delim = ','
        else:
            with open(supplemental_file, 'r') as sf:
                first_line = sf.readline()
            delim = '\t' if first_line.count('\t') > first_line.count(',') else ','
        # Read the supplemental file into a dictionary
        with open(supplemental_file, 'r') as sf:
            reader = csv.DictReader(sf, delimiter=delim)
            if reader.fieldnames:
                lower_fields = [h.lower().strip() for h in reader.fieldnames]
                for field in optional_fields:
                    if field in lower_fields:
                        optional_fields_present.append(field)
                    else:
                        # check for case variants of field name
                        for lf in lower_fields:
                            if lf == field:
                                optional_fields_present.append(field)
                                break
                # Preserve the defined order of optional fields
                optional_fields_present = [f for f in optional_fields if f in optional_fields_present]

                # Any remaining column (besides the id column and the fields
                # the pipeline itself computes) is a free-form extra field:
                # keep its header text exactly as the user wrote it, in file
                # order, and skip duplicate/blank header names defensively.
                reserved_lower = set(optional_fields) | RESERVED_METADATA_FIELDS | {id_column.lower().strip()}
                seen_extra_lower = set()
                for header_name in reader.fieldnames:
                    header_lower = header_name.lower().strip()
                    if not header_lower or header_lower in reserved_lower or header_lower in seen_extra_lower:
                        continue
                    seen_extra_lower.add(header_lower)
                    extra_fields_present.append(header_name)
            # Map sample ID to its supplemental metadata
            for row in reader:
                sample_id = (row.get(id_column) or row.get(id_column.capitalize()) or row.get(id_column.lower()))
                if not sample_id:
                    continue
                sample_id = str(sample_id)
                sup_data[sample_id] = {}
                for field in optional_fields_present:
                    val = None
                    # Find value in row (case-insensitive match)
                    for col, value in row.items():
                        if col.lower().strip() == field:
                            val = value.strip() if isinstance(value, str) else value
                            break
                    sup_data[sample_id][field] = val if val not in ['', None] else None
                for field in extra_fields_present:
                    sup_data[sample_id][field] = sanitize_supplemental_value(row.get(field))

    # Prepare output directory if needed
    out_dir = os.path.dirname(output_prefix)
    base_name = os.path.basename(output_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Open drop log for recording samples that are skipped
    drop_log_path = os.path.join(out_dir, "metadata_drop.txt") if out_dir else "metadata_drop.txt"
    drop_log = open(drop_log_path, 'w')

    # Prepare aggregated output TSV writer
    agg_path = os.path.join(out_dir, f"{base_name}.tsv") if out_dir else f"{base_name}.tsv"
    agg_file = open(agg_path, 'w', newline='')
    writer = csv.writer(agg_file, delimiter='\t')



    # Determine base header fields (for bacterial organisms vs viruses)
    if organism in ['salmonella', 'campylobacter', 'escherichia']:
        # Required field
        header = ["strain", "Serovar", "MLST", "cgMLST"] + [f"HC{level}" for level in HIERCC_LEVELS]

        # Add any optional metadata fields present in supplemental file
        header += optional_fields_present
        # Pass through any other user-supplied supplemental columns as-is
        # (e.g. age, gender); see sanitize_supplemental_value().
        header += extra_fields_present

        # extra columns were selected from full WGS json
        extra_columns = []
        if extra_fields:
            # Antibiotic resistance columns (8 antibiotics * 4 attributes each)
            antibiotics = ['Azithromycin', 'Ciprofloxacin', 'Clindamycin', 'Erythromycin',
                           'Chloramphenicol', 'Tetracycline', 'Gentamicin', 'Nalidixic Acid']
            for ab in antibiotics:
                extra_columns += [
                    f"{ab}_opornos", f"{ab}_czynnik_typ",
                    f"{ab}_czynnik_nazwa", f"{ab}_czynnik_mutacja"
                ]
            # Contamination columns for common programs
            for prog in ["kraken2", "metaphlan", "kmerfinder"]:
                extra_columns += [
                    f"{prog}_species_main", f"{prog}_species_secondary",
                    f"{prog}_species_main_value", f"{prog}_species_secondary_value"
                ]
            # FastQC stats columns (forward and reverse)
            extra_columns += [
                "number_of_reads_forward", "number_of_bases_forward",
                "reads_median_quality_forward", "reads_median_length_forward",
                "number_of_reads_reverse", "number_of_bases_reverse",
                "reads_median_quality_reverse", "reads_median_length_reverse",
                "cgMLST_public_id", "MLST_public_id"
            ]
        header += extra_columns
        # Write header to aggregated TSV
        writer.writerow(header)


        # Process each sample directory in output_dir
        for sample in os.listdir(output_dir):
            sample_dir = os.path.join(output_dir, sample)
            if not os.path.isdir(sample_dir):
                continue
            sample_id = sample
            json_file_path = os.path.join(sample_dir, f"{sample}.json")
            if not os.path.isfile(json_file_path):
                drop_log.write(f"{sample_id} - JSON file not found\n")
                continue

            # Load the JSON data for this sample
            try:
                with open(json_file_path, 'r') as jf:
                    data = json.load(jf)
            except Exception as e:
                drop_log.write(f"{sample_id} - JSON parse error: {e}\n")
                continue


            reasons = []  # to collect any drop reasons for this sample

            # Determine Serovar (species or serotype) based on genus
            serovar_dict = get_serovar_bacteria(data)

            serovar = serovar_dict.get("serovar")
            reasons.extend(serovar_dict.get("reasons"))

            # Extract MLST, chMLST, HierCC clustering

            extract_mlst_cgmlst_out = get_mlst_cgmlst(data)

            mlst_id = extract_mlst_cgmlst_out.get("mlst_id")
            cgmlst_id = extract_mlst_cgmlst_out.get('cgmlst_id')
            cgmlst_public_id = extract_mlst_cgmlst_out.get('cgmlst_public')
            mlst_public_id = extract_mlst_cgmlst_out.get('mlst_public')
            hiercc_values = {level: extract_mlst_cgmlst_out.get(f'hc{level}') for level in HIERCC_LEVELS}



            reasons.extend(extract_mlst_cgmlst_out.get('reasons'))

            # If any required data was missing or invalid, log and skip this sample
            if reasons:
                drop_log.write(f"{sample_id} - {'; '.join(reasons)}\n")
                continue

            # Assemble the row dictionary for this sample
            row = {
                "strain": sample_id,
                "Serovar": serovar,
                "MLST": mlst_id,
                "cgMLST": cgmlst_id,
            }
            for level in HIERCC_LEVELS:
                row[f"HC{level}"] = hiercc_values[level]
            # Include any supplemental metadata fields for this sample
            for field in optional_fields_present:
                val = sup_data.get(sample_id, {}).get(field)
                row[field] = "" if val is None else str(val)
            # Extra free-form supplemental columns are already sanitized and
            # default to "N/A" when missing for this sample.
            for field in extra_fields_present:
                row[field] = sup_data.get(sample_id, {}).get(field, "N/A")

            # Include extra fields if the flag is set
            if extra_fields:
                # AMR data (resistance to known antibiotics)
                # cgmlst / mlst public
                row['cgMLST_public_id'] = cgmlst_public_id
                row['MLST_public_id'] = mlst_public_id

                amr_results = get_amr_bacteria(data)

                for antibiotic, result in amr_results.items():
                    # Flatten each antibiotic's result into separate columns
                    row[f"{antibiotic}_opornos"] = result.get('opornos', '')
                    row[f"{antibiotic}_czynnik_typ"] = result.get('czynnik_typ', '')
                    row[f"{antibiotic}_czynnik_nazwa"] = result.get('czynnik_nazwa', '')
                    row[f"{antibiotic}_czynnik_mutacja"] = result.get('czynnik_mutacja', '')


                # Contamination data
                contam_results = get_contaminations_bacteria(data)
                for key, val in contam_results.items():
                    row[key] = "" if val is None else str(val)


                # FastQC statistics (forward and reverse if applicable)
                fastqc_forward = get_fastqc_stats(data, "forward")
                fastqc_reverse = get_fastqc_stats(data, "reverse")

                fastqc_stats = {}
                # Merge forward and reverse stats
                fastqc_stats.update(fastqc_forward)
                fastqc_stats.update(fastqc_reverse)
                # If no R1/R2 entries (e.g., single-end/Nanopore), use first entry as forward
                if not fastqc_stats:
                    seq_data = data.get("output", {}).get("sequencing_summary_data", [])
                    if seq_data:
                        first_entry = seq_data[0]
                        fastqc_stats = {
                            "number_of_reads_forward": first_entry.get('number_of_reads_value', -1),
                            "number_of_bases_forward": first_entry.get('number_of_bases_value', -1),
                            "reads_median_quality_forward": first_entry.get('reads_median_quality_value', -1),
                            "reads_median_length_forward": first_entry.get('reads_median_length_value', -1),
                            # Leave reverse fields empty in single-end case
                            "number_of_reads_reverse": "",
                            "number_of_bases_reverse": "",
                            "reads_median_quality_reverse": "",
                            "reads_median_length_reverse": ""
                        }
                # Add FastQC stats to the row
                for key, val in fastqc_stats.items():
                    row[key] = "" if val == "" else str(val)

            # Write this sample's data to the aggregated TSV (ensure correct column order)
            writer.writerow([row.get(col, "") for col in header])
            # Also write an individual TSV for the sample (with the same columns)
            sample_out_path = os.path.join(out_dir, f"{sample_id}.{base_name}.tsv") if out_dir else f"{sample_id}.{base_name}.tsv"
            with open(sample_out_path, 'w', newline='') as sf:
                sw = csv.writer(sf, delimiter='\t')
                sw.writerow(header)
                sw.writerow([row.get(col, "") for col in header])

            if with_fasta:
                # collect all .fasta files for the sample
                fasta_files = glob.glob(os.path.join(sample_dir, "fastas", "*.fasta"))

                # ensure output dir exists
                out_dir_for_fastas = Path(out_dir) / "fastas"
                out_dir_for_fastas.mkdir(parents=True, exist_ok=True)

                # open compressed output file
                out_file = out_dir_for_fastas / f"{sample_id}.fasta.gz"
                with gzip.open(out_file, "wt") as sf:
                    for fasta_file in fasta_files:
                        with open(fasta_file, "r") as f:
                            sf.write(f.read())


    elif organism in ['sars-cov-2', 'influenza', 'rsv']:
        # For viral organisms, implementation would go here (not provided in this context)
        # Creating a header
        header = ["strain", "virus", 'type']


        # Add any optional metadata fields present in supplemental file
        header += optional_fields_present
        # Pass through any other user-supplied supplemental columns as-is
        # (e.g. age, gender); see sanitize_supplemental_value().
        header += extra_fields_present
        extra_columns = []
        if extra_fields:
            # Resistance data (only for influenza)
            if organism == 'influenza':
                antivirals = ['Oseltamivir', 'Zanamivir', 'Peramivir', 'Laninamivir', 'Baloxavir']

                for av in antivirals:
                    extra_columns += [f"{av}_resistance_status", f"{av}_mutation",]

            # Contamination viral species are analyzed only with kraken2
            for prog in ["kraken2"]:
                extra_columns += [
                    f"{prog}_species_main", f"{prog}_species_secondary",
                    f"{prog}_species_main_value", f"{prog}_species_secondary_value"
                ]
            # FastQC stats columns (forward and reverse)
            extra_columns += [
                "number_of_reads_forward", "number_of_bases_forward",
                "reads_median_quality_forward", "reads_median_length_forward",
                "number_of_reads_reverse", "number_of_bases_reverse",
                "reads_median_quality_reverse", "reads_median_length_reverse",
            ]

            # freyja output
            extra_columns += ['freyja_lineage_main', 'freyja_lineage_secondary', 'freyja_lineage_main_value', 'freyja_lineage_secondary_value']

            # classification pango and/or nextclade species dependent
            if organism == 'sars-cov-2':
                extra_columns += ['nextclade_variant_name', 'pangolin_variant_name']
            elif organism == 'influenza':
                extra_columns += ['nextclade_variant_name_HA', 'nextclade_variant_name_NA']
            elif organism == 'rsv':
                extra_columns += ['nextclade_variant_name']
            # genome assembly stats
            extra_columns += ['average_coverage', 'total_length', 'number_of_Ns']

        header += extra_columns
        # Write header to aggregated TSV
        writer.writerow(header)


        # Now we iterate over each json in a dir to preapre a valid output
        for sample in os.listdir(output_dir):

            sample_dir = os.path.join(output_dir, sample)
            if not os.path.isdir(sample_dir):
                continue
            sample_id = sample
            json_file_path = os.path.join(sample_dir, f"{sample}.json")

            if not os.path.isfile(json_file_path):
                drop_log.write(f"{sample_id} - JSON file not found\n")
                continue
                # Load the JSON data for this sample
            try:
                with open(json_file_path, 'r') as jf:
                    data = json.load(jf)
            except Exception as e:
                drop_log.write(f"{sample_id} - JSON parse error: {e}\n")
                continue

            # Obligatory data

            get_viral_obligatory_data_out = get_viral_obligatory_data(data)

            row = {
                "strain": get_viral_obligatory_data_out.get('strain'),
                "virus": get_viral_obligatory_data_out.get('virus'),
                "type": get_viral_obligatory_data_out.get('type'),
            }

            if extra_fields:

                # Extract FastQC output
                fastqc_forward = get_fastqc_stats(data, "forward")
                fastqc_reverse = get_fastqc_stats(data, "reverse")
                fastqc_stats = {}
                # Merge forward and reverse stats
                fastqc_stats.update(fastqc_forward)
                fastqc_stats.update(fastqc_reverse)
                # If no R1/R2 entries (e.g., single-end/Nanopore), use first entry as forward
                if not fastqc_stats:
                    seq_data = data.get("output", {}).get("sequencing_summary_data", [])
                    if seq_data:
                        first_entry = seq_data[0]
                        fastqc_stats = {
                            "number_of_reads_forward": first_entry.get('number_of_reads_value', -1),
                            "number_of_bases_forward": first_entry.get('number_of_bases_value', -1),
                            "reads_median_quality_forward": first_entry.get('reads_median_quality_value', -1),
                            "reads_median_length_forward": first_entry.get('reads_median_length_value', -1),
                            # Leave reverse fields empty in single-end case
                            "number_of_reads_reverse": "",
                            "number_of_bases_reverse": "",
                            "reads_median_quality_reverse": "",
                            "reads_median_length_reverse": ""
                        }
                # Add FastQC stats to the row
                for key, val in fastqc_stats.items():
                    row[key] = "" if val == "" else str(val)

                # Add antivirals only for influenza
                if organism == 'influenza':
                    get_influenza_antiviral_data_out = get_influenza_antiviral_data(data)

                    for key, val in get_influenza_antiviral_data_out.items():
                        row[key] = "" if val == "" else str(val)

                # Add kraken2 data and add them to a row

                get_viral_kraken2_data_out = get_viral_kraken2_data(data)
                for key, val in get_viral_kraken2_data_out.items():
                    row[key] = "" if val == "" else str(val)

                # Add freyja2 lineage abundance
                get_viral_freyja_data_out = get_viral_freyja_data(data)
                for key, val in get_viral_freyja_data_out.items():
                    row[key] = "" if val == "" else str(val)

                get_viral_classification_data_out = get_viral_classification_data(data)
                for key, val in get_viral_classification_data_out.items():
                    row[key] = "" if val == "" else str(val)

                get_viral_genome_data_out = get_viral_genome_data(data)
                for key, val in get_viral_genome_data_out.items():
                    row[key] = "" if val == "" else str(val)


            # Include any supplemental metadata fields for this sample
            for field in optional_fields_present:
                val = sup_data.get(sample_id, {}).get(field)
                row[field] = "" if val is None else str(val)
            # Extra free-form supplemental columns are already sanitized and
            # default to "N/A" when missing for this sample.
            for field in extra_fields_present:
                row[field] = sup_data.get(sample_id, {}).get(field, "N/A")

            # Determine if sample should be dropped
            # Any field is Unknown or -1 (defualts for error when parsing json)
            valid_sample = True
            error_msg = ''
            for klucz, wartosc in row.items():
                if wartosc == -1 or wartosc == 'Unknwon':
                    error_msg += f'For column {klucz} the value is incorrect {wartosc}'
                    valid_sample = False
            if not valid_sample:
                drop_log.write(f"{sample_id} -  {error_msg}\n")
                continue

            # Write this sample's data to the aggregated TSV (ensure correct column order)
            writer.writerow([row.get(col, "") for col in header])
            # Also write an individual TSV for the sample (with the same columns)
            sample_out_path = os.path.join(out_dir, f"{sample_id}.{base_name}.tsv") if out_dir else f"{sample_id}.{base_name}.tsv"
            with open(sample_out_path, 'w', newline='') as sf:
                sw = csv.writer(sf, delimiter='\t')
                sw.writerow(header)
                sw.writerow([row.get(col, "") for col in header])

            if with_fasta:
                # for viruses there might be one or many fasta files (depending on number of segments)
                # phylo pipeline now simply requires a gzip fasta file with one-or-many sequences

                fasta_files = sorted(glob.glob(os.path.join(sample_dir, "*.fasta")))

                if not fasta_files:
                    drop_log.write(f"{sample_id} - no FASTA files found in {sample_dir}\n")
                    continue

                # Put output of this script in output_dir/fastas to avoid mixing with metadata
                out_dir_for_fastas = Path(out_dir) / "fastas"
                out_dir_for_fastas.mkdir(parents=True, exist_ok=True)

                # Define gzipped output file
                out_fasta_file = out_dir_for_fastas / f"{sample_id}.fasta.gz"

                # Merge all FASTAs into one gzipped file
                with gzip.open(out_fasta_file, "wt") as f_out:
                    for fasta_file in fasta_files:
                        with open(fasta_file, "r") as f_in:
                            shutil.copyfileobj(f_in, f_out)

                # Optional: verify compression worked (not required but clean)
                if not out_fasta_file.exists():
                    drop_log.write(f"{sample_id} - failed to write {out_fasta_file}\n")

    # Close output files
    agg_file.close()
    drop_log.close()

if __name__ == '__main__':
    generate_metadata()
