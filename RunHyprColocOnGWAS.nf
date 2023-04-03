#!/usr/bin/env nextflow

def helpmessage() {

log.info"""

HyprColocForEqtlCatalog v${workflow.manifest.version}"
===========================================================
Pipeline for running HyprColoc colocalisation analyses (https://www.nature.com/articles/s41467-020-20885-8) between GWAS locus and eQTL Catalogue datasets (https://www.ebi.ac.uk/eqtl/).

Usage:

nextflow run RunHyprColocOnGWAS.nf --gwas --window --Pthresh --OutputDir

Mandatory arguments:
--gwas                    GWAS summary statistics file (by default, GWAMA format).

Optional arguments:
--GwasType                Whether your GWAS trait is binary (case-control) or continuous. One of two: `cc` or `cont`. Defaults to `cc` 
--OutputDir               Directory for colocalisation results and SNP PIPs from each GWAS-eQTL colocalisation. Defaults to "results" folder in work directory.
--window                  Window around each lead SNP to test the colocalization. Defaults to 1000000 (+/- 1Mb from lead SNP).
--Pthresh                 P-value threshold to identify lead SNPs. Defaults to 5e-8.
--MafThresh               MAF threshold to filter the input GWAS data. Defaults to 0.01.
--PosteriorThreshold      Posterior probability to declare colocalization. Defaults to 0.8.
--CsThreshold             Threshold for credible set calculation. Defaults to 0.95 (95% credible sets).
--OutputCsPip             Whether to write out PIPs for each colocalising locus: TRUE/FALSE, defaults to FALSE.
--eQtlFile                Local file listing eQTL Catalogue files. This file can be used to run analysis on the local data folder (not to download from ftp). In the latter case, the folder has to also contain .tbi index for each listed file.
--RegionList              File with regions to include to the analysis. Has to be in hg38 and in format 1:100000-12342145.
--LiftOver                Whether to lift GWAS data into hg38 build. Defaults to "no", functional options are "Hg19ToHg38" and "Hg18ToHg38".
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.gwas = ''
params.GwasType = 'cc'
params.window = 1000000
params.Pthresh = 5e-8
params.OutputDir = ''
params.MafThresh = 0.01
params.PosteriorThreshold = 0.8
params.CsThreshold = 0.95
params.OutputCsPip = false
params.eQtlFile = 'none'
params.RegionList = 'none'
params.LiftOver = 'no'

if (params.OutputCsPip == false){OutputCsPip = "FALSE"} else {OutputCsPip = "TRUE"}

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
summary['GWAS type']                                = params.GwasType
summary['Window']                                   = params.window
summary['P-value threshold']                        = params.Pthresh
summary['MAF threshold']                            = params.MafThresh
summary['Posterior probability threshold']          = params.PosteriorThreshold
summary['Credible set threshold']                   = params.CsThreshold
summary['Write out PIPs']                           = params.OutputCsPip
summary['eQTL file']                                = params.eQtlFile
summary['Region list']                              = params.RegionList
summary['Liftover']                                 = params.OutputDir
summary['Output directory']                         = params.LiftOver

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
InpGwas = Channel.fromPath(params.gwas).ifEmpty { exit 1, "No input GWAS summary statistics found" } 
OutpDir = Channel.fromPath(params.OutputDir)

Window = Channel.value(params.window)
Pthresh = Channel.value(params.Pthresh)
MafThresh = Channel.value(params.MafThresh)
PosteriorThreshold = Channel.value(params.PosteriorThreshold)
CsThreshold = Channel.value(params.CsThreshold)

// Here specify optionally file with the list of regions to include to the analysis
if(params.RegionList != 'none'){
  Channel.fromPath(params.RegionList)
  .ifEmpty { error "Cannot find regions from file: ${params.RegionList}" }
  .splitCsv(header: true, sep: '\t', strip: true)
  .map{row -> [ row.regions ]}
  .set { Regions }
}


// Here specify the list of eQTL Catalogue ftp files to query

if(params.eQtlFile == 'none'){
  Channel.fromPath("$baseDir/bin/eQTLCatalogue_datasets.txt")
  .ifEmpty { error "Cannot find regions from file: $baseDir/data/eQTLCatalogue_datasets.txt" }
  .splitCsv(header: false, sep: '\t', strip: true)
  .map{row -> tuple(row)}
  .set { eqtl_files }
} else {
  Channel.fromPath(params.eQtlFile)
  .ifEmpty { error "Cannot find regions from file: ${params.eQtlFile}" }
  .splitCsv(header: false, sep: '\t', strip: true)
  .map{row -> tuple(row)}
  .set { eqtl_files }
}

process AdjustGwasFile {

    tag {"$gwas"}

    input:
      path gwas from InpGwas
      val id from params.id
      val chr from params.chr
      val pos from params.pos
      val effect_allele from params.effect_allele
      val other_allele from params.other_allele
      val effect from params.effect
      val standard_error from params.standard_error
      val p_value from params.p_value
      val allele_frequency from params.allele_frequency
      val liftover from params.LiftOver

    output:
      file "adjusted_gwas.txt" into GwasAdjusted


    """
    Rscript --vanilla $baseDir/bin/AdjustGwasFile.R ${gwas} ${id} ${chr} ${pos} ${effect_allele} ${other_allele} ${effect} ${standard_error} ${p_value} ${allele_frequency} ${liftover}
    """
}

GwasAdjusted.into{InpGwasRegions; InpGwasSubset}

process FindRegions {

    tag {"$gwas"}

    when:
      params.RegionList == 'none'

    input:
      path gwas from InpGwasRegions
      val Win from Window
      val Pthresh from Pthresh
      val maf from params.MafThresh

    output:
      file "Regions.txt" into LeadSNPsFile

    """
    Rscript --vanilla $baseDir/bin/IdentifyLeadSnps.R ${gwas} ${maf} ${Win} ${Pthresh}
    """
}

if(params.RegionList == 'none'){
LeadSNPsFile
.splitCsv(header: true, sep: '\t', strip: true)
.map{row -> tuple(row.Region)}
.flatten()
.set { LeadSnpCh }
} else {
Channel.fromPath(params.RegionList)
.ifEmpty { error "Cannot find regions from file: ${params.RegionList}" }
.splitCsv(header: false, sep: '\t', strip: true)
.map{row -> tuple(row)}
.set { LeadSnpCh }
}


LeadSnpCh.combine(InpGwasSubset).into{LeadSnpCh1; LeadSnpCh2; LeadSnpCh3}

process SubsetGwas {

    tag {"$Region"}

    input:
      tuple val(Region), file(InpGwas) from LeadSnpCh1
      val Maf from params.MafThresh

    output:
      tuple val(Region), file("*.out.gz") into GwasSubset

    """
    Rscript --vanilla $baseDir/bin/ProcessGwasFile.R ${InpGwas} ${Region} ${Maf} ${Region}.out
    gzip -f ${Region}.out
    """
}

LeadSnpCh2.combine(eqtl_files).into{eqtl_files2; eqtl_files3}

process DownloadEqtl {
    
    tag {"$Region $eqtl_file"}

    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val(Region), file(gwas_file), val(eqtl_file) from eqtl_files3

    output:
      tuple val(Region), path("*tsv.gz.txt.gz") into eqtl_region

    shell:
    '''
    dat_name="$(echo !{eqtl_file})"
    dat_name=${dat_name}_eqtl_region.txt

    echo Processing: !{eqtl_file}

    FileName="$(echo !{eqtl_file} |sed -e "s/.*\\///")"
    echo ${FileName}

    
    len_file2=2

    tabix !{eqtl_file} !{Region} > ${FileName}.txt
    len_file1="$(wc -l ${FileName}.txt | awk '{print $1}')"

    # Make sure that downloaded file is not empty, if exit with error
    if [[ $len_file1 -gt $len_file2 ]]
      then
        gzip ${FileName}.txt
        echo "Download successful!"
      else
        echo "Download NOT successful!"
        exit 1
    fi
    '''
}

eqtl_region = eqtl_region.combine(GwasSubset, by: 0)

process RunHyprColoc {
    
    tag {"$Region $eqtl_region"}

    input:
      tuple val(Region), file(eqtl_region), file(GwasSubset) from eqtl_region
      val posterior from PosteriorThreshold
      val cs from CsThreshold
      val GwasType from params.GwasType
      val OutputCsPip from OutputCsPip

    output:
      path "output.txt" into eqtl_coloc
      path "SNP_PIPs.txt" into SnpPips

    script:
    if (OutputCsPip == "TRUE")
      """
      Rscript --vanilla $baseDir/bin/CombineRunHyprColoc.R \
      ${GwasSubset} \
      ${GwasType} \
      ${eqtl_region} \
      ${Region} \
      ${OutputCsPip} \
      ${posterior} \
      ${cs} \
      output.txt
      """
    else if (OutputCsPip == "FALSE")
      """
      Rscript --vanilla $baseDir/bin/CombineRunHyprColoc.R \
      ${GwasSubset} \
      ${GwasType} \
      ${eqtl_region} \
      ${Region} \
      ${OutputCsPip} \
      ${posterior} \
      ${cs} \
      output.txt

      touch SNP_PIPs.txt
      """
}

eqtl_coloc.flatten().collectFile(name: 'HyprColocResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
if (OutputCsPip == "TRUE"){SnpPips.flatten().collectFile(name: 'GwasColocSnpPipResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
