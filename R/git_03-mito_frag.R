# file: git_03-mito_frag.R
# author: Derek Wong, Ph.D
# date: June 8th, 2021

## Extract Mitochondrial reads
mito <- granges(keepSeqlevels(galp, paste0("chrM"), pruning.mode="coarse"),
                 on.discordant.seqnames="drop")
rm(galp)
mt_nFrag <- length(mito)
mt_width <- width(mito)
rm(mito)

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