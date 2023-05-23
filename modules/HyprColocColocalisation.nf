#!/bin/bash nextflow


process SampleOverlapMatrix {

    publishDir "${params.output}/sample_overlap", mode: 'copy', overwrite: true

    input:
        path inclusionDir

    output:
        path "sample_overlap_matrix.csv"

    script:
        """
        sample_overlap_table.py ${inclusionOutput}/filter_logs_full.log ${inclusionOutput}/filter_logs_full.log
        """
}


process RunCisTransHyprColoc {

    input:
        tuple val(locus_string), path(empirical), path(permuted)
        path geneCorrelations
        path sampleOverlap
        val posteriorThreshold
        val csThreshold
        val outputCsPip

    output:
        path "cis_trans_coloc.${locus_string}.cs.txt", emit: cs
        path "cis_trans_coloc.${locus_string}.variant_pips.txt", emit: pips

    script:
    """
    Rscript --vanilla $baseDir/bin/RunHyprColoc.R \
    --empirical ${empirical} \
    --permuted ${permuted} \
    --assoc ${gwas} \
    --gene-correlations ${geneCorrelations} \
    --sample-overlap ${sampleOverlap} \
    --output-cs-pip ${outputCsPip} \
    --posterior-threshold ${posteriorThreshold} \
    --cs-threshold ${csThreshold} \
    --output-prefix "cis_trans_coloc.${locus_string}"
    """
    if (outputCsPip == "FALSE")
      """
      touch cis_trans_coloc.${locus_string}.variant_pips.txt
      """
}

process RunHyprColoc {

    input:
        tuple val(locus_string), path(empirical), path(permuted)
        path geneCorrelations
        path sampleOverlap
        val posteriorThreshold
        val csThreshold
        val outputCsPip

    output:
        path "cis_trans_coloc.${locus_string}.cs.txt", emit: cs
        path "cis_trans_coloc.${locus_string}.variant_pips.txt", emit: pips

    script:
    """
    Rscript --vanilla $baseDir/bin/RunHyprColoc.R \
    --empirical ${empirical} \
    --permuted ${permuted} \
    --assoc ${gwas} \
    --gene-correlations ${geneCorrelations} \
    --sample-overlap ${sampleOverlap} \
    --output-cs-pip ${outputCsPip} \
    --posterior-threshold ${posteriorThreshold} \
    --cs-threshold ${csThreshold} \
    --output-prefix "cis_trans_coloc.${locus_string}"
    """
    if (outputCsPip == "FALSE")
      """
      touch cis_trans_coloc.${locus_string}.variant_pips.txt
      """
}



workflow CIS_TRANS_COLOCALIZATION {
    take:
        results_grouped_ch
        gene_correlations_ch
        inclusion_step_output_ch
        posterior_threshold
        cs_threshold
        output_cs_pip

    main:
        // Calculate the sample overlap matrix
        sample_overlap_ch = SampleOverlapMatrix(inclusion_step_output_ch)

        // For each cluster of overlapping significant loci, determine if the genes colocalize
        hypr_coloc_results = RunCisTransHyprColoc(
            results_grouped_ch, gene_correlations_ch, sample_overlap_ch, posterior_threshold, cs_threshold, output_cs_pip)

    emit:
        hypr_coloc_results
}


workflow GWAS_COLOCALIZATION {
    take:
        results_grouped_ch
        posterior_threshold
        cs_threshold
        output_cs_pip

    main:
        // For each cluster of overlapping significant loci, determine if the genes colocalize
        hypr_coloc_results = RunGwasHyprColoc(
            results_grouped_ch, posterior_threshold, cs_threshold, output_cs_pip)

    emit:
        hypr_coloc_results
}