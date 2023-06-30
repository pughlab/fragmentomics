# file: runFrag.R
# author: Derek Wong, Ph.D
# date: Aug 3rd, 2022

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--fasta_5"), type = "character", help = "Path to fasta file. Required."),
  make_option(c("--fasta_3"), type = "character", help = "Path to fasta file. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=999, stringsAsFactors=F)

## Load required packages
library(tidyverse)
library(data.table)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

## Get variables from input script
id <- opt$id
fasta_file_5 <- opt$fasta_5
fasta_file_3 <- opt$fasta_3
outdir <- opt$outdir

## Read in fasta
fasta_5 <- read.delim(fasta_file_5, header = FALSE)
fasta_3 <- read.delim(fasta_file_3, header = FALSE)

### Prune and format fastas
fasta_5$V4 <- toupper(fasta_5$V4)
fasta_5 <- as.data.table(fasta_5)
fasta_5 <- fasta_5[!(fasta_5$V4 %like% "N"), ]
fasta_5$length <- nchar(fasta_5$V4)
#fasta_5 <- fasta_5[fasta_5$length == 4, ]

fasta_3$V4 <- toupper(fasta_3$V4)
fasta_3 <- as.data.table(fasta_3)
fasta_3 <- fasta_3[!(fasta_3$V4 %like% "N"), ]
fasta_3$length <- nchar(fasta_3$V4)
#fasta_3 <- fasta_3[fasta_3$length == 4, ]

### Reverse 3' end
#fasta_3 <- fasta_3[1:1000, ] ### Use this for quick tests
splits <- strsplit(fasta_3$V4, "")
reversed <- lapply(splits, rev)
concat <- lapply(reversed, function(x){paste(x, collapse = "")})
fasta_3$V4 <- unlist(concat)
fasta_3$V4 <- chartr("ATGC","TACG", fasta_3$V4)

### Count the nucleotides
fasta <- bind_rows(fasta_5, fasta_3)
motif <- strsplit(fasta$V4, "")
motif <- do.call(rbind, motif)

counts <- data.frame(nucleotide = c("A", "T", "G", "C"))

for (i in 1:ncol(motif)) {
  x <- as.data.frame(table(motif[, i]))
  
  counts <- merge(counts, x, by.x = "nucleotide", by.y = "Var1", all = TRUE)
}

counts[is.na(counts)] <- 0
colnames(counts) <- c("nucleotide", seq(-15, -15 + (ncol(counts) - 2)))

total <- sum(counts$`-15`)

frequency <- counts[, 2:ncol(counts)]/total
ratio <- frequency/c(0.2951856, 0.2039184, 0.2047607, 0.2961353)

frequency$nucleotide <- counts$nucleotide
ratio$nucleotide <- counts$nucleotide

frequency <- frequency %>%
  select(nucleotide, everything())
ratio <- ratio %>%
  select(nucleotide, everything())

### Write file
write.table(counts, file.path(outdir, paste0(id, "_count.txt")), sep = "\t", row.names = FALSE)
write.table(frequency, file.path(outdir, paste0(id, "_frequency.txt")), sep = "\t", row.names = FALSE)
write.table(ratio, file.path(outdir, paste0(id, "_ratio.txt")), sep = "\t", row.names = FALSE)
