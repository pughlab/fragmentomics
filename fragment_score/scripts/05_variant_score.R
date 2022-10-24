library(optparse)
library(vcfR)
library(stringr)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--bam"), type = "character", help = "Path to bam file. Required."),
  make_option(c("--ref"), type = "character", help = "Path to Reference Set. Required."),
  make_option(c("--vcf"), type = "character", help = "Path to vcf. Required."),
  make_option(c("--type"), type = "character", help = "germline or somatic. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

## Get variables from input script
id <- opt$id
bam_file <- opt$bam
ref_file <- opt$ref
vcf_file <- opt$vcf
libdir <- opt$libdir
outdir <- opt$outdir
type <- opt$type

## Import functions
source(paste0(libdir,"/functions.R"))

## Read reference
ref <- read.delim(ref_file, header = FALSE)
ref <- as.vector(ref$V1)

## Read and format VCF
vcf <- read.vcfR(vcf_file)
vcf <- cbind(as.data.frame(vcf@fix), as.data.frame(vcf@gt))
vcf$POS <- as.numeric(vcf$POS)
x <- str_split_fixed(vcf[,10], ":", n = Inf)[, 2]
v1 <- as.numeric(sapply(str_split(x, ",", n = Inf), "[[", 1))
v2 <- as.numeric(sapply(str_split(x, ",", n = Inf), "[[", 2))
vcf$vaf <- v2/(v1 + v2)
vcf$ALT_F1R2 <- str_split_fixed(vcf[,10], ":", n = Inf)[, 4]
vcf$ALT_F2R1 <- str_split_fixed(vcf[,10], ":", n = Inf)[, 5]
vcf <- vcf[vcf$ALT_F1R2 > 0 & vcf$ALT_F2R1 > 0, ]

if (type == "germline") {
  vcf <- vcf[(vcf$vaf > 0.45 & vcf$vaf < 0.55), ]
}

## Run script
variant_scores <- GenerateVariantFS(ref, vcf, bam_file, id)
write.table(variant_scores, file.path(outdir, paste0(id, "_", type, ".txt")), sep = "\t", row.names = FALSE)

q('no')
