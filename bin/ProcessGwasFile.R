#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Load packages and functions ----
library(data.table)
library(dplyr)
library(tidyverse)

setDTthreads(4)

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
  temp_table$id <- rowidv(temp_table, cols=c("chr_pos"))

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

############
# Analysis #
############

# Reading in GWAS file ----
gwas <- fread(args[1])

# MAF filter ----
MafThreshold <- as.numeric(args[3])
gwas <- gwas[gwas$EAF > as.numeric(MafThreshold) & gwas$EAF < (1 - as.numeric(MafThreshold))]

gwas$POS <- as.numeric(gwas$POS)

# Define locus and settings ----
chr <- str_replace(args[2], ":.*", "")
pos <- str_replace(args[2], ".*:", "")
pos1 <- as.numeric(str_replace(pos, "-.*", ""))
pos2 <- as.numeric(str_replace(pos, ".*-", ""))

if(is.na(pos1)){stop("Region is not in the correct format!")}
if(is.na(pos2)){stop("Region is not in the correct format!")}
if(!as.numeric(chr) %in% c(1:22)){stop("Position chromosome is not 1:22!")}

# Data processing ----
gwas <- gwas[as.character(gwas$CHR) == as.character(chr)]
gwas <- gwas[!is.na(gwas$CHR)]
gwas <- gwas[gwas$CHR %in% c(as.character(c(1:22)), "X", "Y")]

gwas$MAF <- gwas$EAF
gwas[gwas$EAF > 0.5, ]$MAF <- 1 - gwas[gwas$EAF > 0.5, ]$MAF

gwas <- gwas[gwas$CHR == chr & gwas$POS > pos1 & gwas$POS < pos2]

gwas$SNP <- ConvUniqSNPName(chr = gwas$CHR, pos = gwas$POS, allele1 = gwas$EA, allele2 = gwas$NEA)

# Remove SNPs where beta or se(beta) is NA 
message(paste(nrow(gwas[is.na(gwas$BETA) | is.na(gwas$SE), ]), "removed"))
gwas <- gwas[!(is.na(gwas$BETA) | is.na(gwas$SE)), ]
gwas <- gwas[!(gwas$BETA == 0 | gwas$SE == 0), ]

fwrite(gwas, args[4], sep = '\t', quote = FALSE, row.names = FALSE)

message("GWAS file subsetted!")
message(paste(nrow(gwas), "lines in GWAS file!"))
