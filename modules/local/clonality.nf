process CLONALITY {
    label 'process_single'
    publishDir "${params.outdir}/clonality", mode: 'copy'

    input:
    path gt

    output:
    path "pairwise.tsv",     emit: pairwise
    path "group_splits.tsv", emit: splits, optional: true
    path "summary.txt",      emit: summary

    script:
    def grp = params.group ? "--group '${params.group}'" : ""
    """
    clonality_metrics.py --genotypes ${gt} \\
        --min-dp ${params.min_dp} --max-dp-factor ${params.max_dp_factor} ${grp}
    """
}
