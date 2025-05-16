process find_identical_sequences {
    tag "${segmentId}"
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    tuple val(segmentId), path(fasta)

    output:
    tuple val(segmentId), path("valid_strains_unique.fasta"), emit: uniq_fasta
    tuple val(segmentId), path("valid_strains_ident_seq.csv"), emit: duplicated_ids

    script:
    """
    find_identical_sequences.py --input ${fasta} --output_dir .
    """
}