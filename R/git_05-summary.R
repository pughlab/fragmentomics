# file: git_05-5Mb_bins.R
# author: Derek Wong, Ph.D
# date: June 16th, 2021

## Load in healthy controls
healthy <- file.path(healthy)
load(healthy)

sd_dist <- function(sample, mean, sd) {
  sds <- (sample - mean)/sd
}

## Generate correlations
correlations <- df.fr3 %>% ungroup() %>% group_by(sample_id, seqnames, arm, start) %>%
  dplyr::summarize(ratio_cor=cor(ratio, healthy_median$median_ratio, method="pearson", use="complete.obs"),
            ratio_corrected_cor=cor(ratio_corrected, healthy_median$median_corrected_ratio, method="pearson", use="complete.obs"),
            ratio_centered_cor=cor(ratio_centered, healthy_median$median_centered_ratio, method="pearson", use="complete.obs"),
            coverage_cor=cor(coverage, healthy_median$median_coverage, method="pearson", use="complete.obs"),
            coverage_corrected_cor=cor(coverage_corrected, healthy_median$median_coverage_corrected, method="pearson", use="complete.obs"),
            coverage_centered_cor=cor(coverage_centered, healthy_median$median_coverage_centered, method="pearson", use="complete.obs"),
            combined_cor=cor(combined, healthy_median$median_combined, method="pearson", use="complete.obs"),
            combined_centered_cor=cor(combined_centered, healthy_median$median_combined_centered, method="pearson", use="complete.obs"),
            nfrags = sum(nfrags),
            mode_size=unique(mode_size),
            mean_size=unique(mean_size),
            median_size=unique(median_size),
            q25_size=unique(q25_size),
            q75_size=unique(q75_size),
            hqbases_analyzed = 100*sum(nfrags)*2,
            depth = hqbases_analyzed/(504*5e6)
  )

## Generate distance from median
distance <- df.fr3 %>% ungroup() %>% group_by(sample_id, seqnames, arm, start) %>%
  dplyr::summarize(ratio_sd=sd_dist(ratio, healthy_median$ratio, healthy_median$ratio_sd),
                   ratio_corrected_sd=sd_dist(ratio_corrected, healthy_median$ratio_corrected, healthy_median$ratio_corrected_sd),
                   ratio_centered_sd=sd_dist(ratio_centered, healthy_median$ratio_centered, healthy_median$ratio_centered_sd),
                   coverage_sd=sd_dist(coverage, healthy_median$coverage, healthy_median$coverage_sd),
                   coverage_corrected_sd=sd_dist(coverage_corrected, healthy_median$coverage_corrected, healthy_median$coverage_corrected_sd),
                   coverage_centered_sd=sd_dist(coverage_centered, healthy_median$coverage_centered, healthy_median$coverage_centered_sd),
                   combined_sd=sd_dist(combined, healthy_median$combined, healthy_median$combined_sd),
                   combined_centered_sd=sd_dist(combined_centered, healthy_median$combined_centered, healthy_median$combined_centered_sd),
                   nfrags = sum(nfrags),
                   mode_size=unique(mode_size),
                   mean_size=unique(mean_size),
                   median_size=unique(median_size),
                   q25_size=unique(q25_size),
                   q75_size=unique(q75_size),
                   hqbases_analyzed = 100*sum(nfrags)*2,
                   depth = hqbases_analyzed/(504*5e6)
  )

## Generate summaries
summary_cor <- df.fr3 %>% ungroup() %>% group_by(sample_id) %>%
  dplyr::summarize(ratio_cor=cor(ratio, healthy_median$median_ratio, method="pearson", use="complete.obs"),
                   ratio_corrected_cor=cor(ratio_corrected, healthy_median$median_corrected_ratio, method="pearson", use="complete.obs"),
                   ratio_centered_cor=cor(ratio_centered, healthy_median$median_centered_ratio, method="pearson", use="complete.obs"),
                   coverage_cor=cor(coverage, healthy_median$median_coverage, method="pearson", use="complete.obs"),
                   coverage_corrected_cor=cor(coverage_corrected, healthy_median$median_coverage_corrected, method="pearson", use="complete.obs"),
                   coverage_centered_cor=cor(coverage_centered, healthy_median$median_coverage_centered, method="pearson", use="complete.obs"),
                   combined_cor=cor(combined, healthy_median$median_combined, method="pearson", use="complete.obs"),
                   combined_centered_cor=cor(combined_centered, healthy_median$median_combined_centered, method="pearson", use="complete.obs"),
                   nfrags = sum(nfrags),
                   mode_size=unique(mode_size),
                   mean_size=unique(mean_size),
                   median_size=unique(median_size),
                   q25_size=unique(q25_size),
                   q75_size=unique(q75_size),
                   hqbases_analyzed = 100*sum(nfrags)*2,
                   depth = hqbases_analyzed/(504*5e6)
  )

summary_sd <- df.fr3 %>% ungroup() %>% group_by(sample_id) %>%
  dplyr::summarize(ratio_sd=mean(sd_dist(ratio, healthy_median$ratio, healthy_median$ratio_sd)),
                   ratio_corrected_sd=mean(sd_dist(ratio_corrected, healthy_median$ratio_corrected, healthy_median$ratio_corrected_sd)),
                   ratio_centered_sd=mean(sd_dist(ratio_centered, healthy_median$ratio_centered, healthy_median$ratio_centered_sd)),
                   coverage_sd=mean(sd_dist(coverage, healthy_median$coverage, healthy_median$coverage_sd)),
                   coverage_corrected_sd=mean(sd_dist(coverage_corrected, healthy_median$coverage_corrected, healthy_median$coverage_corrected_sd)),
                   coverage_centered_sd=mean(sd_dist(coverage_centered, healthy_median$coverage_centered, healthy_median$coverage_centered_sd)),
                   combined_sd=mean(sd_dist(combined, healthy_median$combined, healthy_median$combined_sd)),
                   combined_centered_sd=mean(sd_dist(combined_centered, healthy_median$combined_centered, healthy_median$combined_centered_sd)),
                   nfrags = sum(nfrags),
                   mode_size=unique(mode_size),
                   mean_size=unique(mean_size),
                   median_size=unique(median_size),
                   q25_size=unique(q25_size),
                   q75_size=unique(q75_size),
                   hqbases_analyzed = 100*sum(nfrags)*2,
                   depth = hqbases_analyzed/(504*5e6)
  )

summary_df <- cbind(summary_cor, summary_df)
colnames(summary_df) <- c("correlation_to_median", "sd_from_median")
rownames(summary_df) <- c("ratio", "ratui_corrected", "ratio_centered", "coverage", "coverage_corrected",
                          "coverage_centered", "combined", "combined_centered", "nfrags" ,"mode_size",
                          "mean_size", "median_size", "q25_size", "q75_size", "hqbases_analyzed",
                          "depth")

## Write summary tables
write.table(correlations, file.path(outdir, paste0(id, "_summarycor.txt")), sep = "\t")
write.table(distance, file.path(outdir, paste0(id, "_summarysd.txt")), sep = "\t")
write.table(summary_df, file.path(outdir, paste0(id, "_summary.txt")), sep = "\t")

rm(correlations, distance, summary_df, sd_dist)
