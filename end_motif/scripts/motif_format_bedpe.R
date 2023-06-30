# file: runFrag.R
# author: Derek Wong, Ph.D
# date: Aug 3rd, 2022

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--bedpe"), type = "character", help = "Path to bedpe file. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=999, stringsAsFactors=F)

## Load required packages
library(tidyverse)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

## Get variables from input script
id <- opt$id
bedpe_file <- opt$bedpe
outdir <- opt$outdir

## Read in bedpe
bedpe <- read.delim(bedpe_file, header = FALSE)

### Prune and format bedpe
## Same chromosome and proper orientation
chrs <- paste0("chr", c(1:22, "X", "Y"))
bedpe <- bedpe[bedpe$V1 %in% chrs, ]
bedpe <- bedpe[bedpe$V1 == bedpe$V4, ]
bedpe <- bedpe[!(bedpe$V9 == "-"), ]

## Fragments less than 600bp
bedpe$length <- bedpe$V6 - bedpe$V2
bedpe <- bedpe[bedpe$length <= 600, ]

## Get fragment starts and ends
bedpe <- bedpe[, c("V1", "V2", "V6")]
bedpe <- bedpe[order(factor(bedpe$V1, levels = chrs),
                     bedpe$V2), ]
bedpe$Start <- ifelse(bedpe$V2 < bedpe$V6, bedpe$V2, bedpe$V6) # - 1
bedpe$End <- ifelse(bedpe$V6 > bedpe$V2, bedpe$V6, bedpe$V2) # - 1

### Get last 4 bases
bedpe$front <- bedpe$Start + 4
bedpe$back <- bedpe$End - 4

bedpe_1 <- bedpe[, c("V1", "Start", "front")]
bedpe_2 <- bedpe[, c("V1", "back", "End")]

### Write file
write.table(bedpe_1, file.path(outdir, paste0(id, "_5.bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bedpe_2, file.path(outdir, paste0(id, "_3.bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
