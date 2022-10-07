# file: runFrag.R
# author: Derek Wong, Ph.D
# date: Aug 3rd, 2022

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--fasta"), type = "character", help = "Path to fasta file. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=999, stringsAsFactors=F)

## Load required packages
library(tidyverse)
library(data.table)
library(Biostrings)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

## Get variables from input script
id <- opt$id
fasta_file <- opt$fasta
outdir <- opt$outdir

## Read in fasta
fasta <- read.delim(fasta_file, header = FALSE)

### Prune and format fasta
fasta$V4 <- toupper(fasta$V4)
fasta <- as.data.table(fasta)
fasta <- fasta[!(fasta$V4 %like% "N"), ]

### Find contexts
#fasta <- fasta[1:1000, ] ### Use this for quick tests
sequences <- DNAStringSet(fasta$V4)
contexts <- dinucleotideFrequency(sequences, step = 1, with.labels = TRUE)
contexts <- as.data.frame(t(contexts))

sums <- data.frame(context = row.names(contexts),
                   count = rowSums(contexts))

sums_freq <- data.frame(context = row.names(contexts),
                        freq = sums$count/sum(sums$count))

### Write file
write.table(sums_freq, file.path(outdir, paste0(id, "_contexts.txt")), sep = "\t", row.names = FALSE)
write.table(sums, file.path(outdir, paste0(id, "_raw.txt")), sep = "\t", row.names = FALSE)
