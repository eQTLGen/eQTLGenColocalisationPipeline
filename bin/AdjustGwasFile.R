#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

message(args[1])
out_file <- "adjusted_gwas.txt"
print(args)

id <- args[2]
chr <- args[3]
pos <- args[4]
effect_allele <- args[5]
other_allele <- args[6]
effect <- args[7]
standard_error <- args[8]
p_value <- args[9]
allele_frequency <- args[10]
liftover <- args[11]
if (length(args) == 12) {
    out_file <- args[12]
}

library(data.table)
library(rtracklayer)
library(stringr)

# Read in GWAS sumstats file
gwas <- fread(args[1])

message("Data read in!")

if (length(gwas[[id]]) < 1) { message(paste(id, "is not present in the file!")) }
if (length(gwas[[chr]]) < 1) { message(paste(chr, "is not present in the file!")) }
if (length(gwas[[pos]]) < 1) { message(paste(pos, "is not present in the file!")) }
if (length(gwas[[effect_allele]]) < 1) { message(paste(effect_allele, "is not present in the file!")) }
if (length(gwas[[other_allele]]) < 1) { message(paste(other_allele, "is not present in the file!")) }
if (length(gwas[[p_value]]) < 1) { message(paste(p_value, "is not present in the file!")) }
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
P = gwas[[p_value]],
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
