process treetime {
    input:
    path alignment
    path metadata
    path tree

    output:
    path "timetree.nexus", emit: out

    script:
    """
    treetime --aln ${alignment} \
             --dates ${metadata} \
             --tree ${tree} \
             --keep-polytomies
    mv *_treetime/* .
    rmdir *_treetime
    """
}