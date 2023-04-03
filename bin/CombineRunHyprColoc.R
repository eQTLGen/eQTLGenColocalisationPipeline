#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

gwas <- args[1]
GwasType <- args[2]
eqtl <- args[3]
region <- args[4]
cred <- args[5]
PosteriorThreshold <- args[6]
CsThreshold <- args[7]
output <- args[8]

# Load packages ----
library(data.table)
library(tidyverse)
library(dplyr)
library(Matrix)
library(hyprcoloc)

# Parse region
CHR <- str_replace(region, ":.*", "")
POS <- str_replace(region, ".*:", "")
POS1 <- as.numeric(str_replace(POS, "-.*", ""))
POS2 <- as.numeric(str_replace(POS, ".*-", ""))

# Define functions ----
#' Convert SNP information to more unique SNP ID.
#'
#' This function takes SNP chromosome, position, allele information and constructs the new SNP ID in the form chr:pos_FirsAllele_SecondAllele, where alleles are ordered alphabetically. This is useful when comparing data (e.g. effect directions) from SNPs which come from different studies, potentially using different dbSNP builds. Pay attention that function does not take care of differences in genome build (hg18, hg19, hg38). It also does not take care of potential strand issues, therefore care is needed when working with older studies and A/T or C/G SNPs.
#'
#' @param chr Vector with chromosome position for each SNP. Expects 1:22, X, Y, MT. Chr/chr at the beginning is allowed but is removed from output SNP ID.
#' @param pos Vector with genomic positions for each SNP. Function does not check any genomic build (hg18, hg19, hg38) or whether data is in 0-based/1-based format.
#' @param allele1 Vector with first allele for each SNP. Program does not make any checks, so indels are also supported. "Missing" alleles should be denoted as "-" (E.g. alleles: "-/A")
#' @param allele2 Vector with second allele for each SNP.
#'
#' @return Returns vector with new SNP IDs in the format chr:pos_OneAllele_SecondAllele, where alleles are ordered alphabetically. This function might be useful when one needs to compare SNPs from different studies.
#' @export
#'
#' @examples
#' sample_data <- data.frame(
#'   rs = c("rs100", "rs101", "rs102", "rs103", "rs104", "rs105", "rs100"),
#'   chr = c("1", "5", "X", "12", "chr3", "chr3", "1"),
#'   pos = c(10034, 1341351, 13515, 23413153, 1342525, 1342525, 10034),
#'   all1 = c("T", "A", "C", "G", "A", "A", "T"),
#'   all2 = c("C", "T", "G", "A", "-", "AT", "C")
#' )
#'
#' sample_data$SNPID <- ConvUniqSNPName(chr = sample_data$chr, pos = sample_data$pos, allele1 = sample_data$all1, allele2 = sample_data$all2)
ConvUniqSNPName <- function(chr = chr, pos = pos, allele1 = a1, allele2 = a2) {
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(data.table)

  .datatable.aware = TRUE

  # Clean "chr" and/or "Chr" from the start of each chromosome name
  message('Removing "chr/Chr" from chromosome name')
  chr <- str_replace(chr, "chr", "")
  chr <- str_replace(chr, "Chr", "")

  # Make chr to character

  chr <- as.character(chr)

  # Check if chromosomes are 1:22, X, Y, MT
  if (!all(chr %in% c(as.character(1:22), "X", "Y", "MT"))) {
    stop("Error: some of the chromosomes are not 1:22, X, Y, MT. Please check!")
  }

  # Combine chr and pos
  chr_pos <- paste(chr, pos, sep = ":")

  # Sort alleles

  allele1 <- as.character(allele1)
  allele2 <- as.character(allele2)

  # Check if deletions are not annotated as empty
  if (any(allele1 == "") | any(allele1 == " ")) {
    stop("Error: some alleles are empty or designated as NA. Those have to be specified as " - " in the data. Please check!")
  }

  # Check if chromosome position is integer
  if (!all(round(pos) == pos)) {
    stop("Error: chromosome position is not integer")
  }

  temp_table <- data.table(chr_pos = chr_pos, all1 = as.character(allele1), all2 = as.character(allele2))

  # If same position happens >1 in the data, add separate indicator
  temp_table$id <- rowidv(temp_table, cols = c("chr_pos"))

  temp_table <- melt(temp_table, measure.vars = c("all1", "all2"))

  temp_table$value <- as.character(temp_table$value)

  # temp_table <- temp_table %>%
  #   group_by(chr_pos, id) %>%
  #   mutate(sort_alleles = paste(sort(unique(value)), collapse = "_"))

  # Attempt to speed up with data.table
  temp_table <- as.data.table(temp_table)
  temp_table[, sort_alleles := paste(sort(unique(as.character(value))), collapse = "_"), by = .(chr_pos, id)]

  temp_table <- temp_table[temp_table$variable != 'all2', ]

  # new SNP ID
  SNPID <- paste(temp_table$chr_pos, temp_table$sort_alleles, sep = "_")
  return(SNPID)
}


#' Parser for eQTL datasets
#' 
#' Takes eQTL region as input and constructs two tables: in rows are SNPs and in the columns are different QTL features. 
#' In the cells are eQTL betas and se(betas), respectively. Additionally it gives out the table with variant names, ref and alt alleles: 
#' this is useful for aligning those tables with additional datasets.
#'
#' @param eQTL Link to tabixed eQTL dataset (eQTL Catalogue).
#' @param region_chr Chromosome of your region (2:22, X, Y).
#' @param reg_pos1 Start position of the region (hg38!).
#' @param reg_pos2 End position of the region (hg38!).
#' @param overlap_frac Keeps in only those features where the fraction of variants is equal or larger than specified value. Defaults to 0.8.
#'
#' @return List with three tables: betas, se(betas) and ref/alt alleles for each variant.
#' @export
#'
#' @examples
ParseEqtlCatalogueRegion <- function(eQTL, region_chr, reg_pos1, reg_pos2, overlap_frac = 0.8){

library(data.table)
library(tidyverse)
  
  
CHR <- region_chr
POS1 <- reg_pos1
POS2 <- reg_pos2

eqtl <- fread(cmd = paste("gunzip -c", eQTL))
colnames(eqtl) <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

length(unique(eqtl$variant))

# Filter in features which have 80% of SNPs in this window tested
nr_of_snps_per_feature <- as.data.table(table(eqtl$molecular_trait_id))
filtered_features <- nr_of_snps_per_feature[nr_of_snps_per_feature$N >= overlap_frac * length(unique(eqtl$variant))]$V1

eqtl <- eqtl %>% filter(eqtl$molecular_trait_id %in% filtered_features)

ref_alt <- unique(data.table(variant = eqtl$variant, chr = eqtl$chr, pos = eqtl$position, ref = eqtl$ref, alt = eqtl$alt))

# Filter relevant columns
eqtl <- eqtl %>% select(c("molecular_trait_id", "variant", "ref", "alt", "beta", "se"))

# There are duplicate lines with several rsids for same position: take only unique
eqtl <- unique(eqtl)

beta <- eqtl %>% select(c("molecular_trait_id", "variant", "beta"))
se <- eqtl %>% select(c("molecular_trait_id", "variant", "se"))

beta <- spread(beta, molecular_trait_id, beta) %>% drop_na()
se <- spread(se, molecular_trait_id, se) %>% drop_na()
ref_alt <- ref_alt %>% filter(variant %in% beta$variant)

# Order SNP names
beta <- beta[order(beta$variant), ]
se <- se[order(se$variant), ]
ref_alt <- ref_alt[order(ref_alt$variant), ]

return(list(beta = beta, se = se, ref_alt = ref_alt))

}

#' Function for calculating credible sets
#'
#' @param pips Data frame with PIPs. Four columns: dataset, traits, SNP, PIP.
#' @param thresh Threshold for credible sets. Defaults to 0.95 (95%).
#'
#' @return Combined credible sets in the combined tibble.
#' @export
#'
#' @examples
calc_cred_set <- function(input_pips, thresh = 0.95){
  input_pips <- data.table(input_pips)

  res_sort <- input_pips[order(input_pips$PIP, decreasing = TRUE), ]
  cred <- res_sort %>% mutate(cum_PP = cumsum(PIP))

  res <- data.table(locus = "abi", dataset = "abi", traits = "abi", SNP = "abi", PIP = "abi")
  res <- res[-1, ]
  if (nrow(cred) > 0){
  for (i in 1:nrow(cred)){
      if (cred$cum_PP[i] < thresh){
        res <- rbind(res, cred[i, -ncol(cred), with = FALSE])
      } else {
        res <- rbind(res, cred[i, -ncol(cred), with = FALSE])
        break
      }
    }
  }
  return(res)
}

############
# Analysis #
############

gwas <- fread(cmd = paste("gunzip -c", gwas))

res <- ParseEqtlCatalogueRegion(eQTL = eqtl, region_chr = CHR, reg_pos1 = POS1, reg_pos2 = POS2)

if (nrow(res$ref_alt) > 0){

  head(res$ref_alt)

  ### Check if the reference allele match ----
  gwas <- data.table(SNP = gwas$SNP, effect_allele = gwas$EA, non_effect_allele = gwas$NEA, beta = gwas$BETA, se = gwas$SE)

  res$ref_alt$variant <- str_replace(res$ref_alt$variant, "chr", "")
  res$ref_alt$variant <- ConvUniqSNPName(chr = res$ref_alt$chr, pos = res$ref_alt$pos, allele1 = res$ref_alt$ref, allele2 = res$ref_alt$alt)

  # Make help table for comparing effect directions
  help_table <- merge(gwas, res$ref_alt, by.x = "SNP", by.y = "variant")
  help_table[!help_table$effect_allele == help_table$alt, ]$beta <- - help_table[!help_table$effect_allele == help_table$alt, ]$beta
  help_table <- help_table[, c(1, 4, 5), with = FALSE]

  # Merge with original res ----
  res$beta$variant <- res$ref_alt$variant
  res$se$variant <- res$ref_alt$variant
  res$beta <- merge(data.table(SNP = help_table$SNP, GWAS = help_table$beta), res$beta, by.x = "SNP", by.y = "variant")
  res$se <- merge(data.table(SNP = help_table$SNP, GWAS = help_table$se), res$se, by.x = "SNP", by.y = "variant")
  res$ref_alt <- res$ref_alt %>% filter(variant %in% help_table$SNP)

  # This works as HyprColoc input
  beta <- as.matrix(res$beta[, -1])
  se <- as.matrix(res$se[, -1])

  # Specify if your GWAS is for case-control or continuous trait
  bin_outcomes <- rep(0, ncol(beta))
  if (GwasType == "cc"){
  bin_outcomes[1] <- 1
  } else if (GwasType == "cont") {bin_outcomes[1] <- 0} else {
    stop("GWAS has to be cc or cont!")
  }

  message(bin_outcomes)

  # Run Hyprcoloc
  # Sample overlap matrix- GWAS has no sample overlap, all others have all samples overlapping

  # sample_overlap <- matrix(1, nrow = ncol(beta), ncol = ncol(beta))
  # sample_overlap[1, ] <- 0
  # sample_overlap[, 1] <- 0
  # diag(sample_overlap) <- 1

  hyprcoloc_res <- hyprcoloc(
  effect.est = beta,
  effect.se = se,
  binary.outcomes = bin_outcomes,
  # sample.overlap = sample_overlap,
  snp.id = res$ref_alt$variant,
  trait.names = colnames(se),
  snpscores = TRUE,
  bb.alg = TRUE
  )

  h_res <- hyprcoloc_res$results
  eqtl <- str_replace_all(eqtl, ".tsv.gz.*", "")
  eqtl <- str_replace_all(eqtl, ".txt.gz.*", "")
  eqtl <- str_replace_all(eqtl, ".all", "")

  study <- str_replace(eqtl, "_ge.*", "")
  study <- str_replace(study, "_exon.*", "")
  study <- str_replace(study, "_tx.*", "")
  study <- str_replace(study, "_txrev.*", "")
  study <- str_replace(study, "_microarray.*", "")

  feature <- str_replace(eqtl, ".*_ge_.*", "ge")
  feature <- str_replace(feature, ".*_exon_.*", "exon")
  feature <- str_replace(feature, ".*_tx_.*", "tx")
  feature <- str_replace(feature, ".*_txrev_.*", "txrev")
  feature <- str_replace(feature, ".*_microarray_.*", "microarray")

  tissue <- str_replace(eqtl, ".*_ge_", "")
  tissue <- str_replace(tissue, ".*_exon_", "")
  tissue <- str_replace(tissue, ".*_tx_", "")
  tissue <- str_replace(tissue, ".*_txrev_", "")
  tissue <- str_replace(tissue, ".*_microarray_", "")

  h_res <- data.table(locus = region,
  study = study,
  feature = feature,
  tissue = tissue,
  h_res,
  n_snps_tested = nrow(beta))

  # Write out only rows where there is trait "GWAS"
  h_res_output <- h_res[str_detect(h_res$traits, "GWAS") & h_res$posterior_prob >= PosteriorThreshold, ]
  fwrite(h_res_output, output, sep = '\t', quote = FALSE)

  if (cred == TRUE){

    if (nrow(h_res_output[!is.na(h_res_output$posterior_prob), ]) > 0){

    message("Outputting PIPs!")

    traits <- h_res_output$traits

    snpscores <- hyprcoloc_res$snpscores
    snpscores <- snpscores[str_detect(h_res$traits, "GWAS")]
    
    names(snpscores) <- traits

    pips <- data.table(locus = 'temp', dataset = 'temp', traits = 'temp', SNP = 'temp', PIP = 'temp')
    pips <- pips[-1, ]
    for(i in 1:length(snpscores)){
      print(head(pips))
      print(region)
      print(eqtl)
      print(head(names(snpscores)[i]))
      print(head(snpscores[[i]]))
      temp <- data.table(locus = region, dataset = eqtl, traits = names(snpscores)[i], SNP = res$ref_alt$variant, PIP = snpscores[[i]])
      print(head(temp))
      print(ncol(pips))
      print(ncol(temp))
      pips <- rbind(pips, temp)
    }
    pips$PIP <- as.numeric(pips$PIP)

    fwrite(pips, "temp_PIPs.txt", sep = "\t", quote = FALSE)

    cs <- pips %>%
    group_by(locus, dataset, traits) %>%
    do(calc_cred_set(input_pips = ., thresh = as.numeric(CsThreshold))) %>%
    mutate(CS_size = length(unique(SNP)))

    print(head(cs))

    cs$comb <- paste(cs$locus, cs$dataset, cs$traits, cs$SNP, sep = "_")
    pips$comb <- paste(pips$locus, pips$dataset, pips$traits, pips$SNP, sep = "_")

    pips$CS <- "no"
    pips[pips$comb %in% cs$comb, ]$CS <- "yes"


    pips <- merge(pips, cs[, c(7, 6)], by = "comb", all.x = TRUE, allow.cartesian = TRUE)
    colnames(pips) <- c("comb", "locus", "dataset", "traits", "SNP", "PIP", "CS", "CS_size")
    fwrite(pips[, -1, with = FALSE], "SNP_PIPs.txt", sep = '\t', quote = FALSE, row.names = FALSE)
    } else {
    message("No colocalisations with GWAS!")
    pips <- data.table(locus = 'temp', dataset = 'temp', traits = 'temp', SNP = 'temp', PIP = 'temp', CS = 'temp', CS_size = 'temp')
    pips <- pips[-1, ]
    fwrite(pips, "SNP_PIPs.txt", sep = '\t', quote = FALSE, row.names = FALSE)
    }
  }

} else {
  # empty data table
  h_res_output <- data.table(
  locus= NA,
  study = NA,
  feature = NA,
  tissue = NA,
  iteration = NA,
  traits = NA,
  posterior_prob = NA,
  regional_prob = NA,
  candidate_snp = NA,
  posterior_explained_by_snp = NA,
  dropped_trait = NA,
  n_snps_tested = NA)

  h_res_output <- h_res_output[-1, ]

  fwrite(h_res_output, output, sep = '\t', quote = FALSE)
    if (cred == TRUE){
      pips <- data.table(locus = 'temp', dataset = 'temp', traits = 'temp', SNP = 'temp', PIP = 'temp', CS = 'temp', CS_size = 'temp')
      pips <- pips[-1, ]
      fwrite(pips, "SNP_PIPs.txt", sep = '\t', quote = FALSE, row.names = FALSE)
    }
}
