#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(rtracklayer)
library(stringr)
library(argparse)

# Create an argument parser
parser <- ArgumentParser(description = "Script description")

# Add the arguments
parser <- parser$add_argument("-i", "--input", help = "ID argument")
parser <- parser$add_argument("-id", "--id", help = "ID argument")
parser <- parser$add_argument("-chr", "--chr", help = "Chr argument")
parser <- parser$add_argument("-bp", "--pos", help = "Pos argument")
parser <- parser$add_argument("-ea", "--effect_allele", help = "Effect allele argument")
parser <- parser$add_argument("-oa", "--other_allele", help = "Other allele argument")
parser <- parser$add_argument("-b", "--effect", help = "Effect argument")
parser <- parser$add_argument("-se", "--standard_error", help = "Standard error argument")
parser <- parser$add_argument("-p", "--p_value", help = "P-value argument")
parser <- parser$add_argument("-mlp", "--mlog10_p_value", help = "P-value argument")
parser <- parser$add_argument("-af", "--allele_frequency", help = "Allele frequency argument")
parser <- parser$add_argument("-l", "--liftover", help = "Liftover argument")
parser <- parser$add_argument("-out", "--out_file", default = "adjusted_gwas.txt", help = "Output file argument")

# Parse the command-line arguments
args <- parse_args(parser)

# Access the parsed arguments
message(args$args)
out_file <- args$out_file
print(args)
input <- args$input
id <- args$id
chr <- args$chr
pos <- args$pos
effect_allele <- args$effect_allele
other_allele <- args$other_allele
effect <- args$effect
standard_error <- args$standard_error
p_value <- args$p_value
mlog10_p_value <- args$mlog10_p_value
allele_frequency <- args$allele_frequency
liftover <- args$liftover

# Read in GWAS sumstats file
gwas <- fread(input)

message("Data read in!")

if (length(gwas[[id]]) < 1) { message(paste(id, "is not present in the file!")) }
if (length(gwas[[chr]]) < 1) { message(paste(chr, "is not present in the file!")) }
if (length(gwas[[pos]]) < 1) { message(paste(pos, "is not present in the file!")) }
if (length(gwas[[effect_allele]]) < 1) { message(paste(effect_allele, "is not present in the file!")) }
if (length(gwas[[other_allele]]) < 1) { message(paste(other_allele, "is not present in the file!")) }
if (length(gwas[[p_value]]) < 1) { message(paste(p_value, "is not present in the file!")) }
if (length(gwas[[mlog10_p_value]]) < 1) { message(paste(mlog10_p_value, "is not present in the file!")) }
if (length(gwas[[effect]]) < 1) { message(paste(effect, "is not present in the file!")) }
if (length(gwas[[standard_error]]) < 1) { message(paste(standard_error, "is not present in the file!")) }
if (length(gwas[[allele_frequency]]) < 1) { message(paste(allele_frequency, "is not present in the file!")) }

print(head(gwas))

# Change file into format usable by remainder of the pipeline
gwas_reformat <- data.table(MARKERNAME = gwas[[id]],
CHR = gwas[[chr]],
POS = gwas[[pos]],
EA = gwas[[effect_allele]],
NEA = gwas[[other_allele]],
P = ifelse(length(gwas[[p_value]]) > 0, gwas[[p_value]], 10^-gwas[[mlog10_p_value]]),
BETA = gwas[[effect]],
SE = gwas[[standard_error]],
EAF = gwas[[allele_frequency]]
)
print(head(gwas_reformat))
message("Data adjusted!")


if (liftover == "Hg19ToHg38"){

    message("LiftOver hg19 to hg38")

    gwas_reformat$CHR <- as.character(gwas_reformat$CHR)
    gwas_reformat[gwas_reformat$CHR == "23", ]$CHR <- "X"

    gr <- GRanges(seqnames = paste0("chr", gwas_reformat$CHR), strand = "*",
              ranges = IRanges(start = gwas_reformat$POS, width = 1))
    values(gr) <- DataFrame(MARKERNAME = gwas_reformat$MARKERNAME, 
    EA = gwas_reformat$EA, NEA = gwas_reformat$NEA, P = gwas_reformat$P, 
    BETA = gwas_reformat$BETA, SE = gwas_reformat$SE, EAF = gwas_reformat$EAF)

    download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
    system("gunzip hg19ToHg38.over.chain.gz")
    ch <- import.chain("hg19ToHg38.over.chain")

    lifted <- as.data.table(unlist(liftOver(gr, ch)))
    print(head(lifted))

    gwas_reformat <- data.table(
        MARKERNAME = lifted$MARKERNAME,
        CHR = str_replace(lifted$seqnames, "chr", ""),
        POS = lifted$start,
        EA = lifted$EA,
        NEA = lifted$NEA,
        P = lifted$P,
        BETA = lifted$BETA,
        SE = lifted$SE,
        EAF = lifted$EAF)


} else if (liftover == "Hg18ToHg38"){

    message("LiftOver hg18 to hg38")

    gwas_reformat$CHR <- as.character(gwas_reformat$CHR)
    gwas_reformat[gwas_reformat$CHR == "23", ]$CHR <- "X"

    gr <- GRanges(seqnames = paste0("chr", gwas_reformat$CHR), strand = "*",
              ranges = IRanges(start = gwas_reformat$POS, width = 1))

    values(gr) <- DataFrame(MARKERNAME = gwas_reformat$MARKERNAME, 
    EA = gwas_reformat$EA, NEA = gwas_reformat$NEA, P = gwas_reformat$P, 
    BETA = gwas_reformat$BETA, SE = gwas_reformat$SE, EAF = gwas_reformat$EAF)

    download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz", "hg18ToHg38.over.chain.gz")
    system("gunzip hg18ToHg38.over.chain.gz")
    ch <- import.chain("hg18ToHg38.over.chain")

    lifted <- as.data.table(unlist(liftOver(gr, ch)))

    gwas_reformat <- data.table(
        MARKERNAME = lifted$MARKERNAME,
        CHR = str_replace(lifted$seqnames, "chr", ""),
        POS = lifted$start,
        EA = lifted$EA,
        NEA = lifted$NEA,
        P = lifted$P,
        BETA = lifted$BETA,
        SE = lifted$SE,
        EAF = lifted$EAF)
    print(head(gwas_reformat))    
    

}


fwrite(gwas_reformat, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
