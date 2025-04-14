process iqtree {
    input:
    path aln

    output:
    path "${aln}.contree", emit: out

    script:
    """
    iqtree2 -nt AUTO \
        -s ${aln} \
        -m GTR+G \
        -B 1000 \
        -con \
        -minsup 0.75 \
        -redo
    """
}