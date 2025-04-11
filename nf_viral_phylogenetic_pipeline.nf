input_fasta = file(params.input_fasta)
metadata = file(params.metadata)

src_dir = "${baseDir}/src"

include { augur_index_sequences } from './modules/augur_index_sequences.nf'
include { identify_low_quality_sequences } from './modules/identify_low_quality_sequences.nf'
include { augur_filter_sequences } from './modules/augur_filter_sequences.nf'
include { find_identical_sequences } from './modules/find_identical_sequences.nf'

workflow {
    augur_index_sequences(input_fasta)
    identify_low_quality_sequences(augur_index_sequences.out)
    augur_filter_sequences(input_fasta, augur_index_sequences.out, metadata, identify_low_quality_sequences.out)
    find_identical_sequences(augur_filter_sequences.out)
//     align_no_dups(find_identical_sequences.out, params.output_dir)
//     align_with_dups(filter_sequences.out, params.output_dir)
//     build_tree(align_no_dups.out, params.metadata, params.output_dir)
//     insert_duplicates(build_tree.out, find_identical_sequences.ids, params.output_dir)
}

// process align_no_dups {
//     input:
//     path fasta
//     path outdir
//
//     output:
//     path "${outdir}/tree.fasta", emit: out
//
//     script:
//     """
//     augur align --sequences ${fasta} --output ${outdir}/tree.fasta
//     """
// }
//
// process align_with_dups {
//     input:
//     path fasta
//     path outdir
//
//     output:
//     path "${outdir}/tree_dups.fasta", emit: out
//
//     script:
//     """
//     augur align --sequences ${fasta} --output ${outdir}/tree_dups.fasta
//     """
// }
//
// process build_tree {
//     input:
//     path aln
//     path metadata
//     path outdir
//
//     output:
//     path "${outdir}/tree.fasta.contree", emit: out
//
//     script:
//     """
//     iqtree2 -nt AUTO \
//         -s ${aln} \
//         -m GTR+G \
//         -B 1000 \
//         -con \
//         -minsup 0.75 \
//         -redo \
//         --dating LSD \
//         --date ${metadata}
//     """
// }
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
//
// nextflow.config
//
// params.i = null
// params.m = null
// params.o = "results"
//
// process {
//     withLabel: 'local' {
//         executor = 'local'
//     }
// }
