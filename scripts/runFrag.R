# file: runFrag.R
# author: Derek Wong, Ph.D
# date: June 8th, 2021

library(optparse)

option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--bamdir"), type = "character", help = "Path to bam files. Required."),
  make_option(c("--filters"), type = "character", help = "Path to genomic blacklist regions. Required."),
  make_option(c("--gaps"), type = "character", help = "Path to genome gaps. Required."),
  make_option(c("--tiles"), type = "character", help = "Path to 100kb tiled genome. Required."),
  make_option(c("--healthy"), type = "character", help = "Path to panel of healthy controls. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

library(tidyverse)
library(multidplyr)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(GenomicRanges)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
class(Homo.sapiens)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

id <- opt$id
bamdir <- opt$bamdir
filters <- opt$filters
gaps <- opt$gaps
tiles <- opt$tiles
healthy <- opt$healthy
outdir <- opt$outdir
libdir <- opt$libdir

source(paste0(libdir,"/R/git_01-read_bam.R"))
source(paste0(libdir,"/R/git_02-frag.R"))
#source(paste0(libdir,"/R/git_03-mito_frag.R"))
#source(paste0(libdir,"/R/git_04-100kb_bins.R"))
#source(paste0(libdir,"/R/git_05-5Mb_bins.R"))
#source(paste0(libdir,"/R/git_06-summary.R"))
#source(paste0(libdir,"/R/git_07-plotting.R"))
#source(paste0(libdir,"/R/git_08-gbm_prediction.R"))








