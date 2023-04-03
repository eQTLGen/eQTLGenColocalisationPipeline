args = commandArgs(trailingOnly = TRUE)

#' Convert Z-score to P-value
#'
#' This function converts Z-score to two-sided P-value. It also enables to deal with very significant effect sizes for which P-values hit the numerical limit (P<2.22e-308). In those cases it enables to calculate logarithmed P-value instead.
#'
#' @param Z Vector with Z-scores.
#'
#' @param largeZ Boolean which indicates whether the input Z-score vector contains very large absolute Z-scores. These would yield very small P-values hitting the numerical limit (P<2.22e-308), and are rounded to 0. If this is the case, it would make sense to use logartihmed P-values instead. When TRUE and log10P is set to FALSE, it gives the exponent of P-value. Defaults to FALSE.
#' @param log10P Boolean which specifies whether to output -log10(P) instead. Those can then be used for standard visualizations, e.g. on Manhattan plot. Evaluated only if largeZ=TRUE and it defaults to TRUE.
#' @return Vector with P values, log P-values or minus log10 P-values. log P-values can be converted back to P-values with command 10^([-log10(P)]), however P-values smaller than 2.22e-308 will be rounded to 0.
#' @note This implementation is based on the thread here: https://stackoverflow.com/questions/46416027/how-to-compute-p-values-from-z-scores-in-r-when-the-z-score-is-large-pvalue-muc
#' @export
#'
#' @examples
#' z1 <- c(3, 5, 3, 2.3, -3, -0.03, 1.96)
#' z2 <- c(3, 5, 3, 2.3, -3, "-")
#' z3 <- c(3, 5, 3, 2.3, -3, 1.96, -1.2, 0.2, 48)
#'
#' ZtoP(z1)
#' ZtoP(z2)
#' ZtoP(z3)
#' ZtoP(z3, largeZ = TRUE)
ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    message("Converting the non-numeric vector to numeric vector.")
    Z <- as.numeric(Z)
  }
  
  if (largeZ == TRUE) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)
    
    if (largeZ == TRUE & log10P == TRUE) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)
    
    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
    }
  }
  
  return(P)
}

IdentifyLeadSNPs <- function(data, window = 1000000, Pthresh = 5e-8, RemoveHLA = TRUE) {

  data$P <- ZtoP(data$Z)
  data_f <- data[data$P < as.numeric(Pthresh), ]

  if (nrow(data_f) == 0){stop("No  SNPs in summary statistics which are below significance threshold!")}

  
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]
  res$min_reg <- 0
  res$max_reg <- 0
  
  while (nrow(data_f) > 0) {
    
    lead_snp <- data_f[abs(data_f$Z) == max(abs(data_f$Z)), ][1, ] # If multiple SNPs have exactly the same effect, take first

    # Adjust region loci which overlap chromosome boundaries
    min_pos <- min(data[data$chr == lead_snp$chr, ]$pos)
    max_pos <- max(data[data$chr == lead_snp$chr, ]$pos)

    # Region boundaries
    min_reg <- lead_snp$pos - window
    max_reg <- lead_snp$pos + window
    
    if (min_reg <= min_pos){min_reg <- min_pos}
    if (max_reg >= max_pos){max_reg <- max_pos}

    lead_snp$min_reg <- min_reg
    lead_snp$max_reg <- max_reg

    res <- rbind(res, lead_snp)

    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos >= min_reg & data_f$pos <= max_reg),]
    message(paste("Added:", lead_snp$chr, lead_snp$pos))
    
  }
  
  # Remove SNPs for which the region overlaps with HLA region (hg38)
  res <- res[!(res$chr == 6 & ((res$pos - window > 28510120 & res$pos - window < 33480577) | (res$pos + window > 28510120 & res$pos + window < 33480577))), ]
  message("SNPs overlapping hg38 MHC region removed!")

  res <- paste0(res$chr, ":", res$min_reg, "-", res$max_reg)
  
  return(res)
}

library(data.table)

setDTthreads(8)

gwas <- fread(args[1])
message("Data read in!")


gwas$MAF <- gwas$EAF
gwas[gwas$MAF > 0.5, ]$MAF <- 1 - gwas[gwas$MAF > 0.5, ]$EAF
gwas <- gwas[gwas$MAF > as.numeric(args[2])]

gwas$Z <- gwas$BETA / gwas$SE
gwas <- data.table(chr = gwas$CHR, pos = gwas$POS, allele1 = gwas$EA, allele2 = gwas$NEA, Z = gwas$Z)
message("Data prepared!")

LeadSNPs <- IdentifyLeadSNPs(gwas, window = as.numeric(args[3]), Pthresh = as.numeric(args[4]))

fwrite(data.table(Region = LeadSNPs), "Regions.txt", sep = "\t")
