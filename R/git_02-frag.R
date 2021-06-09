# file: git_02-frag.R
# author: Derek Wong, Ph.D
# date: June 9th, 2021

## Filter reads: 90-220bp on Autosomes
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                 on.discordant.seqnames="drop")
w.all <- width(frags)
frags <- frags[which(w.all >= 90 & w.all <= 220)]
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

frags <- frags[-queryHits(findOverlaps(frags, high_cov_regions))]

## Calculate filtered coverage stats
bamCov_filtered <- coverage(frags)
filtered_mean <- mean(bamCov_filtered)
filtered_sd <- sd(bamCov_filtered)
filtered_max <- max(bamCov_filtered)

## Generate QC report
QC_matrix <- data.frame(raw_mean, raw_sd, raw_max, filtered_mean, filtered_sd, filtered_max)
write.table(QC_matrix, file.path(outdir, paste0(id, "_cov_stats.txt")))