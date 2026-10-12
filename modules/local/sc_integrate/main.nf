process SC_INTEGRATE {
    tag "all_samples"
    label 'process_high'

    input:
    path rds_files

    output:
    path "integrated.rds",     emit: rds
    path "cluster_counts.tsv", emit: cluster_counts
    path "umap.png",           emit: umap
    path "sessionInfo.txt",    emit: session_info

    script:
    """
    sc_integrate.R \\
        --batch_var ${params.batch_var ?: 'none'} \\
        --n_pcs ${params.n_pcs} \\
        --resolution ${params.resolution} \\
        ${rds_files}
    """

    stub:
    """
    touch integrated.rds cluster_counts.tsv umap.png sessionInfo.txt
    """
}
