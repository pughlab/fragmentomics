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
options(scipen=0, stringsAsFactors=F)

## Load required packages
library(tidyverse)
library(data.table)
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
fasta$length <- nchar(fasta$V4)
fasta <- fasta[fasta$length == 267, ]

#fasta <- fasta[1:1000, ] ### Use this for quick tests
sequences <- as.data.table(lapply(fasta$V4, function(x) sapply(seq(from=1, to=nchar(x), by=2), function(i) substr(x, i, i+1))))
sequences <- as.data.table(t(sequences))
sequences <- sequences[, 1:(ncol(sequences)-1)]

sequences2 <- sub('.', '', fasta$V4) 
sequences2 <- as.data.table(lapply(sequences2, function(x) sapply(seq(from=1, to=nchar(x), by=2), function(i) substr(x, i, i+1))))
sequences2 <- as.data.table(t(sequences2))

### Count the contexts
contexts <- sapply(sequences, table)
colnames(contexts) <- seq(1, ncol(contexts)*2, by = 2)

context_percent <- sweep(contexts, 2, colSums(contexts), "/")
colnames(context_percent) <- seq(1, 266, by = 2)

contexts2 <- sapply(sequences2, table)
colnames(contexts2) <- seq(2, ncol(contexts)*2, by = 2)

context_percent2 <- sweep(contexts2, 2, colSums(contexts2), "/")
colnames(context_percent2) <- seq(2, 266, by = 2)

### Combine raw contexts table
context_raw <- merge(contexts, contexts2, by = "row.names")
names(context_raw)[names(context_raw) == "Row.names"] <- "context"
order <- as.character(c("context", 1:266))
context_raw <- context_raw[, order]
colnames(context_raw) <- c("context", c(1:266))

### Combine percent contexts table
context_table <- merge(context_percent, context_percent2, by = "row.names")
names(context_table)[names(context_table) == "Row.names"] <- "context"
order <- as.character(c("context", 1:266))
context_table <- context_table[, order]
colnames(context_table) <- c("context", c(1:266))

### Write file
write.table(context_table, file.path(outdir, paste0(id, "_contexts.txt")), sep = "\t", row.names = FALSE)
write.table(context_raw, file.path(outdir, paste0(id, "_raw.txt")), sep = "\t", row.names = FALSE)
