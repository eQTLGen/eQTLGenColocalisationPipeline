#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(IGUtilityPackage)
library(vroom)
library(AnnotationDbi)
library(rtracklayer)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  parser <- ArgumentParser(description = 'Compare eQTls and pQTLs')
  parser$add_argument('--id', metavar = 'str', type = 'character',
                      help = 'Table with pqtl data')
  parser$add_argument('--protein', metavar = 'file', type = 'character',
                      help = 'Table with pqtl data')
  parser$add_argument('--expression', metavar = 'file', type = 'character',
                      help = 'Table with eqtl data')

  args <- parser$parse_args(argv)
  # Process input
  protein_qtls <- fread(args$protein)
  expr_qtls <- fread(args$expression)

  # Perform method

  ManhPlot(protein_qtls,
           InputChrCol = "CHR_num",
           InputPosCol = "BP",
           InputPvalCol = "p-value",
           InputRsCol = "rs_number",
           build = "hg38",
           title = paste(args$id, "pQTL"),
           AnnotateClosestGenes = TRUE,
           gtf = params$gtf, alpha = 1e-200)

  ManhPlot(expr_qtls,
           InputChrCol = "CHR_num",
           InputPosCol = "BP",
           InputPvalCol = "p-value",
           InputRsCol = "rs_number",
           build = "hg38",
           title = paste(args$id, "pQTL"),
           AnnotateClosestGenes = TRUE,
           gtf = params$gtf, alpha = 1e-200)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}