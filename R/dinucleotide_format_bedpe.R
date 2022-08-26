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
chrs <- paste0("chr", c(1:22))
bedpe <- bedpe[bedpe$V1 %in% chrs, ]
bedpe <- bedpe[bedpe$V1 == bedpe$V4, ]
bedpe$length <- bedpe$V6 - bedpe$V2
bedpe <- bedpe[bedpe$length == 167, ]
bedpe <- bedpe[, c("V1", "V2", "V6")]
bedpe <- bedpe[order(factor(bedpe$V1, levels = chrs),
                     bedpe$V2), ]

### Add 50bp padding
bedpe$V2 <- bedpe$V2 - 50
bedpe$V6 <- bedpe$V6 + 50

### Write file
write.table(bedpe, file.path(outdir, paste0(id, ".bed")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)






