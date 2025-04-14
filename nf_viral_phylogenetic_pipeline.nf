input_fasta = file(params.input_fasta)
metadata = file(params.metadata)

src_dir = "${baseDir}/src"

include { augur_index_sequences } from './modules/augur_index_sequences.nf'
include { identify_low_quality_sequences } from './modules/identify_low_quality_sequences.nf'
include { augur_filter_sequences } from './modules/augur_filter_sequences.nf'
include { find_identical_sequences } from './modules/find_identical_sequences.nf'
include { augur_align } from './modules/augur_align.nf'
include { remove_duplicates_from_alignment } from './modules/remove_duplicates_from_alignment.nf'
include { iqtree } from './modules/iqtree.nf'

workflow {
    augur_index_sequences(input_fasta)
    identify_low_quality_sequences(augur_index_sequences.out)
    augur_filter_sequences(input_fasta, augur_index_sequences.out, metadata, identify_low_quality_sequences.out)
    find_identical_sequences(augur_filter_sequences.out)
    augur_align(find_identical_sequences.out.uniq_fasta)
    iqtree(augur_align.out)

//     remove_duplicates_from_alignment(augur_align.out, find_identical_sequences.out.duplicated_ids)
//     insert_duplicates(build_tree.out, find_identical_sequences.ids, params.output_dir)
}

//
// process insert_duplicates {
//     input:
//     path tree
//     path ids
//     path outdir
//
//     output:
//     path "${outdir}/consensus_tree.nwk"
//
//     script:
//     """
//     python src/insert_missing_dupliated_sequences.py \
//         --tree ${tree} \
//         --ids ${ids} > ${outdir}/consensus_tree.nwk
//     """
// }
