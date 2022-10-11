

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--bam"), type = "character", help = "Path to bam file. Required."),
  make_option(c("--ref"), type = "character", help = "Path to Reference Set. Required."),
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
ref <- opt$ref
libdir <- opt$libdir
outdir <- opt$outdir

## Import functions
source(paste0(libdir,"/functions.R"))

## Read in files
reference <- read.delim(ref, header = FALSE)
reference <- as.vector(reference$V1)

## Run script
Patient_score <- GeneratePatientFS(reference, bam)
write.table(Patient_score, file.path(outdir, paste0(id, "_score.txt")), sep = "\t", row.names = FALSE)

q('no')
