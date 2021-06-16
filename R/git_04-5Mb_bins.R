# file: git_05-5Mb_bins.R
# author: Derek Wong, Ph.D
# date: June 16th, 2021

## Generate raw 100kb bin file
tib.list <- as_tibble(AB)
rm(AB)
write.table(tib.list, file.path(outdir, paste0(id, "_100kb_bins.txt")), sep = "\t")
tib.list <- tib.list %>% dplyr::select(-matches("X"))

## Plot GC Correction metrics
pdf(file = file.path(outdir, paste0(id, "_GC_metrics.pdf")))
par(mfrow=c(2,2))
## short
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$short,
              main = "short",
              xlab = "frag_GC", 
              ylab = "short")
## short corrected
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$short.corrected,
              main = "short corrected",
              xlab = "frag_GC", 
              ylab = "short_corrected")
## short vs short predicted
smoothScatter(x = tib.list$short.predicted, 
              y = tib.list$short, 
              main = "short predicted vs actual",
              xlab = "short_predicted", 
              ylab = "short")
## short corrected vs short predicted
smoothScatter(x = tib.list$short.predicted, 
              y = tib.list$short.corrected, 
              main = "short predicted vs corrected",
              xlab = "short_predicted", 
              ylab = "short_corrected")
## long
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$long, 
              main = "long",
              xlab = "frag_GC", 
              ylab = "long")
## long corrected
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$long.corrected, 
              main = "corrected long",
              xlab = "frag_GC", 
              ylab = "long_corrected")
## long vs long predicted
smoothScatter(x = tib.list$long.predicted, 
              y = tib.list$long, 
              main = "long predicted vs actual",
              xlab = "long_predicted", 
              ylab = "long")
## long corrected vs long predicted
smoothScatter(x = tib.list$long.predicted, 
              y = tib.list$long.corrected, 
              main = "long predicted vs corrected",
              xlab = "long_predicted", 
              ylab = "long_corrected")
## total fragments
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$nfrags, 
              main = "nfrags",
              xlab = "frag_GC", 
              ylab = "nfrags")
## corrected total fragments
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$nfrags.corrected, 
              main = "corrected nfrags",
              xlab = "frag_GC", 
              ylab = "nfrags_corrected")
## fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted, 
              y = tib.list$nfrags, 
              main = "nfrags predicted vs actual",
              xlab = "nfrags_predicted", 
              ylab = "nfrags")
## corrected fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted, 
              y = tib.list$nfrags.corrected, 
              main = "nfrags predicted vs corrected",
              xlab = "nfrags_predicted", 
              ylab = "nfrags_corrected")
## ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio, 
              main = "ratios",
              xlab = "frag_gc", 
              ylab = "ratio")
## corrected ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio.corrected, 
              main = "corrected ratios",
              xlab = "frag_gc", 
              ylab = "ratio_corrected")
## predicted ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio.predicted, 
              main = "predicted ratios",
              xlab = "frag_gc", 
              ylab = "ratio_predicted")
## bin GC vs frag GC
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$C.G, 
              main = "GC content",
              xlab = "frag_GC", 
              ylab = "bin_GC")
## coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage, 
              main = "coverage",
              xlab = "frag_GC", 
              ylab = "coverage")
## corrected coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage.corrected, 
              main = "corrected coverage",
              xlab = "frag_GC", 
              ylab = "coverage_corrected")
## predicted coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage.predicted, 
              main = "predicted coverage",
              xlab = "frag_GC", 
              ylab = "coverage_predicted")
dev.off()

## Set arm levels
df.fr2 <- tib.list
rm(tib.list)
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)

## Combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric end and remove the bin closest to the centromere if it is
## smaller than 5mb.
df.fr2 <- df.fr2 %>% group_by(arm) %>%
  mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                          ceiling(rev((1:length(arm))/50) )))

df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
  dplyr::summarize(short2=sum(short, na.rm=TRUE),
            long2=sum(long, na.rm=TRUE),
            short.corrected2=sum(short.corrected, na.rm=TRUE),
            long.corrected2=sum(long.corrected, na.rm=TRUE),
            gc=mean(C.G, na.rm=TRUE),
            frag.gc2=mean(frag.gc, na.rm=TRUE),
            ratio2=mean(ratio, na.rm=TRUE),
            ratio.corrected2=mean(ratio.corrected, na.rm=TRUE),
            nfrags2=sum(nfrags, na.rm=TRUE),
            nfrags.corrected2=sum(nfrags.corrected, na.rm=TRUE),
            coverage2=mean(coverage, na.rm=TRUE),
            coverage.corrected2=mean(coverage.corrected, na.rm=TRUE),
            combined2=mean(combined, na.rm=TRUE),
            short.var=var(short.corrected, na.rm=TRUE),
            long.var=var(long.corrected, na.rm=TRUE),
            nfrags.var=var(nfrags.corrected, na.rm=TRUE),
            mode_size=unique(mode, na.rm=TRUE),
            mean_size=unique(mean, na.rm=TRUE),
            median_size=unique(median, na.rm=TRUE),
            q25_size=unique(quantile.25, na.rm=TRUE),
            q75_size=unique(quantile.75, na.rm=TRUE),
            start=start[1],
            end=rev(end)[1],
            binsize = n())

## Assign bins
df.fr3 <- df.fr3 %>% filter(binsize==50)
df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(bin = 1:length(id))

## Center and scale fragment and coverage ratios
df.fr3$ratio.centered <- ((df.fr3$ratio.corrected2 - mean(df.fr3$ratio.corrected2))/sd(df.fr3$ratio.corrected2))*0.01
df.fr3$coverage.centered <- ((df.fr3$coverage.corrected2 - mean(df.fr3$coverage.corrected2))/sd(df.fr3$coverage.corrected2))*0.01
df.fr3$combined.centered <- ((df.fr3$combined2 - mean(df.fr3$combined2))/sd(df.fr3$combined2))*0.01

### Rename and reorder dataframe
df.fr3 <- df.fr3 %>%
  dplyr::rename(
    sample_id = id,
    short = short2,
    long = long2,
    short_corrected = short.corrected2,
    long_corrected = long.corrected2,
    ratio = ratio2,
    ratio_corrected = ratio.corrected2,
    nfrags = nfrags2,
    nfrags_corrected = nfrags.corrected2,
    coverage_corrected = coverage.corrected2,
    short_var = short.var,
    long_var = long.var,
    nfrags_var = nfrags.var,
    ratio_centered = ratio.centered,
    coverage_centered = coverage.centered,
    combined = combined2,
    combined_centered = combined.centered,
    frag_gc = frag.gc2,
    coverage = coverage2
  )

df.fr3 <- df.fr3 %>%
  relocate(sample_id, seqnames, arm, start, end, gc, frag_gc, short, long, nfrags, ratio,
           short_corrected, long_corrected, nfrags_corrected, ratio_corrected, ratio_centered, 
           coverage, coverage_corrected, coverage_centered, combined, combined_centered,
           short_var, long_var, nfrags_var, mode_size,mean_size, median_size, q25_size, q75_size, 
           binsize, bin)

write.table(df.fr3, file.path(outdir, paste0(id, "_5Mb_bins.txt")), sep = "\t")