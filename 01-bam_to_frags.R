library(optparse)

### Used for getting information from shell script
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample ID"),
  make_option(c("--bamdir"), type = "character", help = "Path to bam files"),
  make_option(c("--outdir"), type = "character", help = "Path to output")
)
parseobj <- OptionParser(option_list = option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen = 0, stringAsFactors = F)

###
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(GenomicRanges)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
class(Homo.sapiens)
options(stringsAsFactors = FALSE)
options(bitmapType = 'cairo')

## Set Variables
id <- opt$id
bamdir <- opt$bamdir
outdir <- opt$outdir

## Read GAlignmentPairs
bamfile <- file.path(bamdir, paste0(id, ".bam"))
indexed.bam <- gsub("$", ".bai", bamfile)
if (!file.exists(indexed.bam)) {
  indexBam(bamfile)
}

param <- ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                                         isProperPair = TRUE,
                                         isDuplicate = FALSE,
                                         isSecondaryAlignment = FALSE,
                                         isUnmappedQuery = FALSE),
                      mapqFilter = 30)

galp <- readGAlignmentPairs(bamfile, param = param)

## Extract Mitochondrial reads
mito <- granges(keepSeqlevels(galp, paste0("chrM"), pruning.mode="coarse"),
                 on.discordant.seqnames="drop")
mt_nFrag <- length(mito)
mt_width <- width(mito)
rm(mito)

## Only keep Autosomes
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                 on.discordant.seqnames="drop")
rm(galp)
gcs <- GCcontent(Hsapiens, unstrand(frags))
frags$gc <- gcs
nFrags <- length(frags)

## Calculate coverage stats
bamCov <- coverage(frags)
raw_mean <- mean(bamCov)
raw_sd <- sd(bamCov)
raw_max <- max(bamCov)

## Set cutoff and filter high coverage areas (VNTRs)
cutoff <- ceiling(median(quantile(bamCov, 0.9999)))
adjust <- ceiling(10*log2(median(raw_mean)+1))
cutoff <- cutoff + ifelse(adjust > 5, adjust, 5)
high_cov_regions <- slice(bamCov, lower = cutoff)
high_cov_regions <- ranges(high_cov_regions)
rm(bamCov)

export.bed(high_cov_regions, file.path(outdir, paste0(id, "_highcovregions.bed")))

high_cov_regions <- read.table(file.path(outdir, paste0(id, "_highcovregions.bed")), col.names = c("chrom", "chromStart", "chromEnd","x", "y", "z"))
high_cov_regions <- high_cov_regions[, -c(4:6)]
high_cov_regions <- makeGRangesFromDataFrame(high_cov_regions, keep.extra.columns=TRUE)

frags_filtered <- frags[-queryHits(findOverlaps(frags, high_cov_regions))]
rm(frags)

## Calculate filtered coverage stats
bamCov_filtered <- coverage(frags_filtered)
filtered_mean <- mean(bamCov_filtered)
filtered_sd <- sd(bamCov_filtered)
filtered_max <- max(bamCov_filtered)
rm(bamCov_filtered)

saveRDS(frags_filtered, file.path(outdir, paste0(id, "_frags.rds")) )

## Calculate mitochondrial stats
mt_median <- median(mt_width)
mt_mean <- mean(mt_width)
dens <- density(mt_width)
mt_mode <- dens$x[which.max(dens$y)]
mt_lower_quartile <- as.integer(quantile(mt_width, 0.25))
mt_upper_quartile <- as.integer(quantile(mt_width, 0.75))
mt_short <- count(mt_width <= 150)
mt_long <- count(mt_width > 150)
mt_ratio <- mt_short/mt_long
mt_coverage <- sum(mt_width)/16569
mt_fraction <- mt_nFrag/nFrags

## Generate mitochondrial report
mt_matrix <- data.frame(mt_short, mt_long, mt_nFrag, mt_ratio, mt_coverage, mt_fraction, mt_median,
                       mt_mean, mt_mode, mt_lower_quartile, mt_upper_quartile)
write.table(mito.list, file.path(outdir, paste0(id, "mito_stats.txt")))

## Generate QC report
QC_matrix <- data.frame(raw_mean, raw_sd, raw_max, filtered_mean, filtered_sd, filtered_max)
write.table(QC_matrix, file.path(outdir, paste0(id, "_cov_stats.txt")))
q('no')
