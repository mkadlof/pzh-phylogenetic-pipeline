// input can be one of the following:
// 1. A single FASTA file containing sequences - for single segment viruses like SARS-CoV-2
// 2. A directory containing multiple FASTA files - for multi-segment viruses like Influenza

input_fasta = file(params.input_fasta) // This var is overloaded (dir or file)
metadata = file(params.metadata)
organism = params.organism

src_dir = "${baseDir}/src"

// Core modules

include { augur_index_sequences } from './modules/augur_index_sequences.nf'
include { identify_low_quality_sequences } from './modules/identify_low_quality_sequences.nf'
include { augur_filter_sequences } from './modules/augur_filter_sequences.nf'
include { find_identical_sequences } from './modules/find_identical_sequences.nf'
include { augur_align } from './modules/augur_align.nf'
include { remove_duplicates_from_alignment } from './modules/remove_duplicates_from_alignment.nf'
include { iqtree } from './modules/iqtree.nf'
include { insert_duplicates_into_tree } from './modules/insert_duplicates_into_tree.nf'
include { insert_duplicates_into_alignment } from './modules/insert_duplicates_into_alignment.nf'
include { treetime } from './modules/treetime.nf'
include { augur_export } from './modules/augur_export.nf'

// influenza specific modules
include { transform_input } from './modules/transform_input.nf'

workflow core {
    augur_index_sequences(input_fasta)
    identify_low_quality_sequences(augur_index_sequences.out)
    augur_filter_sequences(input_fasta, augur_index_sequences.out, metadata, identify_low_quality_sequences.out)
    find_identical_sequences(augur_filter_sequences.out)
    augur_align(find_identical_sequences.out.uniq_fasta)
    iqtree(augur_align.out)
    insert_duplicates_into_tree(iqtree.out, find_identical_sequences.out.duplicated_ids)
    insert_duplicates_into_alignment(augur_align.out, find_identical_sequences.out.duplicated_ids)
    treetime(insert_duplicates_into_alignment.out, metadata, insert_duplicates_into_tree.out)
    augur_export(treetime.out.timetree, metadata, treetime.out.node_data)
}

workflow {
    if (organism.toLowerCase() in ['sars', 'sars2', 'sars-cov-2']) {
        core()
    }
    else if (organism.toLowerCase() in ['flu', 'infl','influenza']) {
        transform_input(input_fasta)
    }
    else if (organism.toLowerCase() in ['rsv']) {
        error "RSV is not supported yet. Please use 'sars-cov-2' or 'influenza'."
    }
    else {
        error "Organism not supported. Please use 'sars-cov-2' or 'influenza'."
    }
}
