#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2021
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(data.table)
library(tidyverse)
library(argparse)
library(dplyr)
library(Matrix)
library(hyprcoloc)

# Declare constants
## Argument parser
parser <- ArgumentParser(description='Run HyprColoc to detect cis-trans colocalizations')
parser$add_argument('--empirical', metavar='file', type = 'character',
                    help='A tab-delimited file with eqtl results for a specified region.')
parser$add_argument('--permuted', metavar='file', type = 'character',
                    help='A tab-delimited file with permuted eqtl results for a specified region.')
parser$add_argument('--assoc', metavar='file', type = 'character',
                    help="A file with GWAS associations")
parser$add_argument('--assoc-desc', metavar='file', type = 'character',
                    help="A tab-delimited file with traitnames for assoc, and gwas types")
parser$add_argument('--gene-correlations', metavar='file', type = 'character',
                    help='A tab-delimited file with gene-gene correlations.')
parser$add_argument('--sample-overlap', metavar='file', type = 'character',
                    help='A tab-delimited file with sample overlap per gene.')
parser$add_argument('--region', metavar='character', type = 'character',
                    help='The region of interest')
parser$add_argument('--output-cs-pip', type = 'character',
                    help='Output credible set probability')
parser$add_argument('--posterior-threshold', type = 'float',
                    help='Posterior probability threshold')
parser$add_argument('--cs-threshold', type = 'float',
                    help='Credible set threshold')
parser$add_argument('--output-prefix', metavar='path', type='character')


# Main

#' Execute main
#' Expected contents of argv:
#'
#' $baseDir/bin/RunHyprColoc.R \
#' --empirical ${empirical} \
#' --permuted ${permuted} \
#' --gene-correlations ${geneCorrelations} \
#' --sample-overlap ${sampleOverlap} \
#' --output-cs-pip ${outputCsPip} \
#' --posterior-threshold ${posteriorThreshold} \
#' --cs-threshold ${csThreshold} \
#' --output-prefix "cis_trans_coloc.${locus_string}"
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Parse command-line arguments
  args <- parser$parse_args(argv)

  # extract parameters
  output_cs_pip <- args$output_cs_pip
  posterior_threshold <- args$posterior_threshold
  cs_threshold <- args$cs_threshold
  output_prefix <- args$output_prefix
  region <- args$region

  # Read eqtl output to do colocalizations for
  empirical <- fread(args$empirical)

  # Read permuted output to do colocalizations for
  permuted <- fread(args$permuted)

  # Specify GWAS
  assoc <- fread(args$assoc)
  assoc_desc <- fread(args$assoc_desc)

  binary_outcomes <- assoc_desc$type

  # Specify sample overlap matrix
  sample_overlap_matrix <- fread(args$sample_overlap)

  # Get zscores as a wide pivotted table
  z_scores <- data.table::dcast(permuted, variant ~ phenotype, value.var = "z_score")

  # Get effect sizes as a wide pivotted table
  eqtl_effect_sizes <- data.table::dcast(empirical, variant ~ phenotype, value.var = "effect_size")

  # Get standard errors around the betas as a wide pivotted table
  eqtl_standard_error <- data.table::dcast(empirical, variant ~ phenotype, value.var = "standard_error")

  # Specify ld matrix
  ld_matrix <- cor(z_scores)^2

  # Run HyprColoc
  hyprcoloc_res <- hyprcoloc(
    effect.est = effect_sizes,
    effect.se = standard_error,
    binary.outcomes = binary_outcomes,
    sample.overlap = sample_overlap_matrix,
    ld.matrix = ld_matrix,
    trait.cor = trait_correlation,
    snp.id = rownames(effect_sizes),
    trait.names = colnames(effect_sizes),
    snpscores = TRUE,
    bb.alg = TRUE
  )

  # Write output for HyprColoc

  # Format output table
  h_res <- data.table(
    locus = region,
    feature = feature,
    tissue = tissue,
    hyprcoloc_res$results,
    n_snps_tested = nrow(beta))

  # Write out only rows where there is trait "GWAS"
  h_res_output <- h_res[h_res$posterior_prob >= posterior_threshold, ]
  fwrite(h_res_output, paste0(output_prefix, ".cs.txt"), sep = '\t', quote = FALSE)

  # If output_cs_pip is true, write credible set posterior probabilities
  if (output_cs_pip == TRUE) {
    cs_pips()
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

