#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { SampleOverlapMatrix; RunCisTransHyprColoc; RunHyprColoc; CIS_TRANS_COLOCALIZATION; GWAS_COLOCALIZATION } from './modules/HyprColocColocalisation.nf'


def helpmessage() {

log.info"""

HyprColocForEqtlCatalog v${workflow.manifest.version}"
===========================================================
Pipeline for running HyprColoc colocalisation analyses (https://www.nature.com/articles/s41467-020-20885-8) between GWAS locus and eQTL Catalogue datasets (https://www.ebi.ac.uk/eqtl/).

Usage:

nextflow run RunHyprColocOnGWAS.nf --gwas --window --Pthresh --OutputDir

Mandatory arguments:
--eqtls                   eQTLGen parquet dataset

Optional arguments:
--gwas                GWAS summary statistics file (by default, GWAMA format).
--gwasType                Whether your GWAS trait is binary (case-control) or continuous. One of two: `cc` or `cont`. Defaults to `cc`
--outputDir               Directory for colocalisation results and SNP PIPs from each GWAS-eQTL colocalisation. Defaults to "results" folder in work directory.
--window                  Window around each lead SNP to test the colocalization. Defaults to 1000000 (+/- 1Mb from lead SNP).
--pThresh                 P-value threshold to identify lead SNPs. Defaults to 5e-8.
--mafThresh               MAF threshold to filter the input GWAS data. Defaults to 0.01.
--posteriorThreshold      Posterior probability to declare colocalization. Defaults to 0.8.
--csThreshold             Threshold for credible set calculation. Defaults to 0.95 (95% credible sets).
--outputCsPip             Whether to write out PIPs for each colocalising locus: TRUE/FALSE, defaults to FALSE.
--regionList              File with regions to include to the analysis. Has to be in hg38 and in format 1:100000-12342145.
--liftOver                Whether to lift GWAS data into hg38 build. Defaults to "no", functional options are "Hg19ToHg38" and "Hg18ToHg38".
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.gwas = 'none'
params.gwas_type = 'cc'
params.window = 1000000
params.p_thresh = 5e-8
params.output_dir = ''
params.maf_thresh = 0.01
params.posterior_threshold = 0.8
params.cs_threshold = 0.95
params.output_cs_pip = false
params.region_list = 'none'
params.liftover = 'no'

params.inclusion_step_output = 'NO_FILE'

output_cs_pips = (params.output_cs_pip == "TRUE")

//Show parameter values
log.info """=======================================================
Colocalisation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['GWAS sumstats']                            = params.gwas
summary['GWAS type']                                = params.gwas_type
summary['Window']                                   = params.window
summary['P-value threshold']                        = params.pt_thresh
summary['MAF threshold']                            = params.maf_thresh
summary['Posterior probability threshold']          = params.posterior_threshold
summary['Credible set threshold']                   = params.cs_threshold
summary['Write out PIPs']                           = params.output_cs_pip
summary['empirical results']                        = params.empirical_results
summary['permuted results']                         = params.permuted_results
summary['Region list']                              = params.region_list
summary['Liftover']                                 = params.liftover
summary['Output directory']                         = params.output_dir

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get empirical channel
empirical_ch = Channel.fromPath(params.empirical)
    .map { file ->
           def key = file.name.toString().tokenize('.').get(1)
           return tuple(key, file) }
    .groupTuple()

// Get permuted channel
permuted_ch = Channel.fromPath(params.permuted)
    .map { file ->
           def key = file.name.toString().tokenize('.').get(1)
           return tuple(key, file) }
    .groupTuple()

// Get gene correlations channel
gene_correlations_ch = Channel.fromPath(params.gene_correlations).collect()

// Get dir with sample sizes per cohort and per gene
inclusion_step_output_ch = file(params.inclusion_step_output)

// Define parameter channels
window = Channel.value(params.window)
p_thresh = Channel.value(params.p_thresh)
maf_thresh = Channel.value(params.maf_thresh)
posterior_threshold = Channel.value(params.posterior_threshold)
cs_threshold = Channel.value(params.cs_threshold)

// Here specify optionally file with the list of regions to include to the analysis
// if(params.bed != 'none'){
//   Channel.fromPath(params.bed)
//   .ifEmpty { error "Cannot find regions from file: ${params.bed}" }
//   .splitCsv(header: true, sep: '\t', strip: true)
//   .map{row -> [ row.regions ]}
//   .set { loci_ch }
// }

// Here specify the list of eQTL Catalogue ftp files to query
// if(params.gwas == 'none'){
//   Channel.fromPath("$baseDir/bin/eQTLCatalogue_datasets.txt")
//   .ifEmpty { error "Cannot find regions from file: $baseDir/data/eQTLCatalogue_datasets.txt" }
//   .splitCsv(header: false, sep: '\t', strip: true)
//   .map{row -> tuple(row)}
//   .set { eqtl_files }
// } else {
//   Channel.fromPath(params.gwas)
//   .ifEmpty { error "Cannot find regions from file: ${params.gwas}" }
//   .splitCsv(header: false, sep: '\t', strip: true)
//   .map{row -> tuple(row)}
//   .set { InpGwas }
// }

// Identify the eQTL blood genes for which we should do colocalisations
// Prepare the eQTL blood genes for which we should do colocalisations

// Identify external phenotypes for which we should do colocalisations
// Prepare the external phenotypes for which we should do colocalisations

// If there are external phenotypes for which we should do colocalisations:
//   Collect combinations of eQTL blood genes for which we should do colocalisations

workflow {
    results_grouped_ch = empirical_ch.join(permuted_ch).view()

    CIS_TRANS_COLOCALIZATION(
        results_grouped_ch,
        gene_correlations_ch, inclusion_step_output_ch,
        posterior_threshold, cs_threshold, output_cs_pips)

    CIS_TRANS_COLOCALIZATION.out.cs.flatten()
        .collectFile(name: 'HyprColocResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")

    if (output_cs_pips == true) {
        CIS_TRANS_COLOCALIZATION.out.pips.flatten()
            .collectFile(name: 'GwasColocSnpPipResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
    }
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
