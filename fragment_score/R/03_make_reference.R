

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--tumour"), type = "character", help = "Path to tumour files. Required"),
  make_option(c("--healthy"), type = "character", help = "Path to healthy files. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

## Get variables from input script
tumour <- opt$tumour
healthy <- opt$healthy
libdir <- opt$libdir
outdir <- opt$outdir

## Import functions
source(paste0(libdir,"/functions.R"))

## Run script
Reference_set <- GenerateReferenceSet(tlen, hlen)
write.table(Reference_set, file.path(outdir, "reference_set.txt"), sep = "\t", row.names = FALSE)

q('no')
