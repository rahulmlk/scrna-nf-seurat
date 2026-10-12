process SC_QC {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("${meta.id}.qc.rds"), emit: rds
    path "${meta.id}.qc_metrics.tsv",           emit: metrics
    path "${meta.id}.sessionInfo.txt",          emit: session_info

    script:
    """
    sc_qc.R \\
        --sample ${meta.id} \\
        --condition ${meta.condition} \\
        --donor ${meta.donor} \\
        --matrix ${matrix_dir} \\
        --min_features ${params.min_features} \\
        --max_mito ${params.max_mito} \\
        --mito_pattern '${params.mito_pattern}'
    """

    stub:
    """
    touch ${meta.id}.qc.rds ${meta.id}.qc_metrics.tsv ${meta.id}.sessionInfo.txt
    """
}
