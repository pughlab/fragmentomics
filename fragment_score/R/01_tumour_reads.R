
library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--bam"), type = "character", help = "Path to bam file. Required."),
  make_option(c("--vcf"), type = "character", help = "Path to vcf file. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

## Get variables from input script
id <- opt$id
bam <- opt$bam
vcf <- opt$vcf
libdir <- opt$libdir
outdir <- opt$outdir

## Import functions
source(paste0(libdir,"/functions.R"))

## Run script
Tumour_Reads <- GetTumourReads(id, vcf, bam)
write.table(Tumour_Reads, file.path(outdir, paste0(id, "_tumour_reads.txt")), sep = "\t", row.names = FALSE)

q('no')
