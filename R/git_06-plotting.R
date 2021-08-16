# file: git_05-5Mb_bins.R
# author: Derek Wong, Ph.D
# date: July 22nd, 2021

## Load in healthy controls
load(healthy)

## Set themes and plot layouts
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
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
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr3$arm <- factor(df.fr3$arm, levels=armlevels)
healthy_median$arm <- factor(healthy_median$arm, levels=armlevels)

arm <- df.fr3 %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

## Plot Fragmentation profile
f1 <- ggplot(df.fr3, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
f1 <- f1 + geom_line(data=healthy_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
f1 <- f1 + labs(x="", y="Fragmentation profile\n", color="")
f1 <- f1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
f1 <- f1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
f1 <- f1 + mytheme
ggsave(file.path(outdir, paste0(id, "_fragment.pdf")), f1, width=15, height=3, units="in")

## Plot short fragment coverage profile
c1 <- ggplot(df.fr3, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
c1 <- c1 + geom_line(data=healthy_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
c1 <- c1 + mytheme
ggsave(file.path(outdir, paste0(id, "_coverage.pdf")), c1, width=15, height=3, units="in")

## Plot combined profile
b1 <- ggplot(df.fr3, aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
b1 <- b1 + geom_line(data=healthy_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
b1 <- b1 + labs(x="", y="Combined profile\n", color="")
b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
b1 <- b1 + mytheme
ggsave(file.path(outdir, paste0(id, "_combined.pdf")), b1, width=15, height=3, units="in")