# file: git_05-5Mb_bins.R
# author: Derek Wong, Ph.D
# date: June 8th, 2021

## Load in healthy controls
healthy <- file.path(healthy)
load(healthy)

## Generate correlations
summary_df <- df.fr %>% ungroup() %>% group_by(sample_id) %>%
  dplyr::summarize(ratio_cor=cor(ratio, healthy_median$median_ratio, method="pearson", use="complete.obs"),
            ratio_corrected_cor=cor(ratio_corrected, healthy_median$median_corrected_ratio, method="pearson", use="complete.obs"),
            ratio_centered_cor=cor(ratio_centered, healthy_median$median_centered_ratio, method="pearson", use="complete.obs"),
            coverage_ratio_cor=cor(coverage, healthy_median$median_coverage_ratio, method="pearson", use="complete.obs"),
            coverage_centered_cor=cor(coverage_centered, healthy_median$median_centered_coverage, method="pearson", use="complete.obs"),
            nfrags = sum(nfrags),
            mode_size=unique(mode_size),
            mean_size=unique(mean_size),
            median_size=unique(median_size),
            q25_size=unique(q25_size),
            q75_size=unique(q75_size),
            hqbases_analyzed = 100*sum(nfrags)*2,
            coverage = hqbases_analyzed/(504*5e6)
  )

summary_df <- inner_join(summary_df, metadata, by=c("sample_id"="TGL_ID"))
#summary_df$`type` = relevel(as_factor(summary_df$`type`), "Healthy")
write.table(summary_df, file.path(outdir, paste0("summary_correlations.txt")), sep = "\t")

# Plot profiles
df.fr %>% group_by(sample_id) %>% filter(bin==1) %>% summarize(n=n())
df.fr$cancer_status <- relevel(factor(df.fr$cancer_status), "negative")
metadata <- metadata[metadata$depth > 0.5,]
df.fr <- df.fr[df.fr$sample_id %in% metadata$TGL_ID, ]
df.fr <- df.fr[df.fr$sample_type == "plasma", ]

mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=12),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(face="bold", size=17),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr$arm <- factor(df.fr$arm, levels=armlevels)
healthy_median$arm <- factor(healthy_median$arm, levels=armlevels)

arm <- df.fr %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

# Generate Fragmentation and Coverage plots
g1 <- ggplot(df.fr, aes(x=bin, y=ratio_centered_01, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.5)
  #geom_line(data=subset(df.fr, opt_TF > 0.03), size=0.75, alpha=1)
g1 <- g1 + geom_line(data=healthy_median, size=0.75, alpha=0.5, color="black")
g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
g1 <- g1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g1 <- g1 + mytheme
g1

c1 <- ggplot(df.fr, aes(x=bin, y=coverage_centered_01, group=sample_id, color="red")) + 
  geom_line(size=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr, opt_TF == 0), size=0.75, alpha=1)
c1 <- c1 + geom_line(data=healthy_median, size=0.75, alpha=0.5, color="black")
c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
c1 <- c1 + mytheme
c1
