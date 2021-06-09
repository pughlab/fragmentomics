library(optparse)

### Used for getting information from shell script
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample ID"),
  make_option(c("--fragdir"), type = "character", help = "Path to frag files"),
  make_option(c("--outdir"), type = "character", help = "Path to output"),
  make_option(c("--filters"), type = "character", help = "Genomic regions to filter out"),
  make_option(c("--gaps"), type = "character", help = "Genomic gaps to filter out"), 
  make_option(c("--tiles"), type = "character", help = "BED file of tiled genome")
)
parseobj <- OptionParser(option_list = option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen = 0, stringAsFactors = F)

###
id <- opt$id
fragpath <- opt$fragdir
outdir <- opt$outdir
filters <- opt$filters
gaps <- opt$gaps
tiles <- opt$tiles

#GC correct function
gc.correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias, na.action = na.omit)
  coverage.model <- loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <- coverage - coverage.pred + median(coverage, na.rm=TRUE)
}

gc.pred <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias, na.action = na.omit)
  coverage.model <- loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred <- predict(coverage.model, bias)
}

###
fragfile <- file.path(fragpath, paste0(id, "_frags.rds"))
filename <- file.path(outdir, paste0(id, "_bin_100kb.rds"))
filters <- file.path(filters)
gaps <- file.path(gaps)
tiles <- file.path(tiles)
if(file.exists(filename)) q('no')

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)
class(Homo.sapiens)
library(devtools)
library(biovizBase)
load(filters)
load(gaps)
hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
options(stringsAsFactors = FALSE)
options(bitmapType = 'cairo')

#Made my own hg38 tiling for this

AB <- read.table(tiles, col.names = c("chrom", "chromStart", "chromEnd", "Seqlength"))
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))

tcmeres <- gaps.hg38[grepl("centromere|telomere", gaps.hg38$type)]

arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg38))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(hsapiens, AB)

## These bins had no coverage
#AB <- AB[-c(8780, 13665)]
fragments <- readRDS(fragfile)
# 
### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg38))]
w.all <- width(fragments)

fragments <- fragments[which(w.all >= 90 & w.all <= 220)]
w <- width(fragments)

frag.list <- split(fragments, w)

counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
if(min(w) > 90) {
  m0 <- matrix(0, ncol=min(w) - 90, nrow=nrow(counts),
               dimnames=list(rownames(counts), 90:(min(w)-1)))
  counts <- cbind(m0, counts)
}

if(max(w) < 220) {
  m1 <- matrix(0, ncol=220 - max(w), nrow=nrow(counts),
               dimnames=list(rownames(counts), (max(w)+1):220))
  counts <- cbind(counts, m1)
}

olaps <- findOverlaps(fragments, AB)
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
bingc <- rep(NA, length(bin.list))
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))
gc <- as.vector(AB$gc)
bingc <- ifelse(is.na(bingc), gc, bingc)
bingc <- ifelse(bingc < min(gc), gc, bingc)
bingc <- ifelse(bingc > max(gc), gc, bingc)

### Get modes
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

short <- rowSums(counts[,1:61])
long <- rowSums(counts[,62:121])
long <- ifelse(long < mean(long)-(sd(long)*3), NA, long)
short <- ifelse(is.na(long), NA, short)
ratio <- short/long
ratio[is.nan(ratio)] <- NA
ratio[is.infinite(ratio)] <- NA
nfrags <- short+long
coverage <- (nfrags)/sum(nfrags, na.rm=TRUE)
short.corrected <- gc.correct(short, bingc)
long.corrected <- gc.correct(long, bingc)
nfrags.corrected <- gc.correct(short+long, bingc)
ratio.corrected <- gc.correct(ratio, bingc)
coverage.corrected <- gc.correct(coverage, bingc)
short.predicted <- gc.pred(short, bingc)
long.predicted <- gc.pred(long, bingc)
nfrags.predicted <- gc.pred(short+long, bingc)
coverage.predicted <- gc.pred(coverage, bingc)

AB$short <- short
AB$long <- long
AB$ratio <- ratio
AB$nfrags <- short+long
AB$coverage <- coverage
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected
AB$coverage.corrected <- coverage.corrected
AB$short.predicted <- short.predicted
AB$long.predicted <- long.predicted
AB$nfrags.predicted <- nfrags.predicted
AB$ratio.predicted <- short.predicted/long.predicted
AB$coverage.predicted <- coverage.predicted

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc
AB$id <- id

for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]

saveRDS(AB, filename)
q('no')