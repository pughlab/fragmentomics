# file: git_01-read_bam.R
# author: Derek Wong, Ph.D
# date: June 8th, 2021

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