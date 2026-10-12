#!/usr/bin/env nextflow

/*
 * scrna-nf-seurat
 * Per-sample single-cell RNA-seq QC (Seurat + scDblFinder), then
 * Harmony integration and clustering. Downstream annotation and
 * pseudobulk DE live in analysis/downstream.qmd.
 */

include { SC_QC        } from './modules/local/sc_qc/main'
include { SC_INTEGRATE } from './modules/local/sc_integrate/main'

workflow {

    if (!params.input) {
        error "Please provide a samplesheet with --input (see README)"
    }

    // matrix_dir in the samplesheet is resolved relative to the samplesheet itself
    def base = file(params.input).parent

    ch_samples = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample, condition: row.condition, donor: row.donor],
                base.resolve(row.matrix_dir)
            )
        }

    SC_QC(ch_samples)

    ch_rds = SC_QC.out.rds
        .map { meta, rds -> rds }
        .collect()

    SC_INTEGRATE(ch_rds)
}
