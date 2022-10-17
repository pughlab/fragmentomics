# file: runFrag.R
# author: Derek Wong, Ph.D
# date: Aug 3rd, 2022

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--path"), type = "character", help = "Path to bedpe file. Required."),
  make_option(c("--peaks"), type = "character", help = "Path to nucleosome peaks. Required"),
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
path <- opt$path
peaks_path <- opt$peaks
outdir <- opt$outdir

## Find files
files <- list.files(path, id, full.names = TRUE)
peaks <- list.files(peaks_path, "chr", full.names = TRUE)
chrs <- files[grep("chr", files)]
chrs <- sapply(strsplit(chrs, "_\\s*"), tail, 1)
chrs <- gsub(".bedpe", "", chrs)

### Read in files for each chromosome
table_prox <- data.frame(length = c(1:600))
table_dist <- data.frame(length = c(1:600))

for (i in c(1:length(chrs))) {
  ### Set variables
  chr <- chrs[[i]]
  file <- files[grep(paste0(chr, ".bed"), files)]
  peak <- peaks[grep(paste0(chr, ".txt"), peaks)]
  
  ### Read in files
  data <- read.delim(file, header = FALSE)
  nucleosomes <- read.delim(peak)
  
  ### Format bedpe file
  data <- data[data$V1 == data$V4, ]
  data <- data[data$V1 == chr, ]
  data$V2 <- ifelse(data$V2 < data$V3, data$V2, data$V3)
  data$V5 <- ifelse(data$V5 < data$V6, data$V5, data$V6)
  data <- data[, c(1,2,6)]
  data$mid <- round((data$V6 + data$V2)/2)
  data$length <- data$V6 - data$V2
  
  nucleosomes <- nucleosomes[nucleosomes$Chrom == chr, ]
  
  ### Find nearest peak
  setDT(data)
  setDT(nucleosomes)
  
  data[,merge:=mid]
  nucleosomes[,merge:=peak]
  
  setkeyv(data, c("merge"))
  setkeyv(nucleosomes, c("merge"))
  
  merged <- nucleosomes[data, roll = "nearest"]
  merged$distance <- abs(merged$mid - merged$peak)
  
  ### Seperate out into reads within 500bp from peak
  proximal <- merged[merged$distance < 500 & merged$length < 600, ]
  distal <- merged[merged$distance >= 500 & merged$length < 600, ]

  ### Count length frequencies
  lengths <- data.frame(length = c(1:600))
  freq_prox <- as.data.frame(table(proximal$length))
  freq_dist <- as.data.frame(table(distal$length))
  
  lengths_prox <- merge(lengths, freq_prox, by.x = "length", by.y = "Var1", all = TRUE)
  lengths_dist <- merge(lengths, freq_dist, by.x = "length", by.y = "Var1", all = TRUE)
  lengths_prox[is.na(lengths_prox)] <- 0
  lengths_dist[is.na(lengths_dist)] <- 0
  
  ### Merge with table
  table_prox <- merge(table_prox, lengths_prox, by = "length")
  table_dist <- merge(table_dist, lengths_dist, by = "length")
}

table_prox$count <- rowSums(table_prox[,-1])
table_prox <- table_prox[, c("length", "count")]

table_dist$count <- rowSums(table_dist[,-1])
table_dist <- table_dist[, c("length", "count")]

write.table(table_prox, file.path(outdir, paste0(id, "_length_proximal.txt")), sep = "\t", row.names = FALSE)
write.table(table_dist, file.path(outdir, paste0(id, "_length_distal.txt")), sep = "\t", row.names = FALSE)


