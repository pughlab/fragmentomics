# Users can choose to either construct their own reference set of tumour vs healthy lengths, OR use the built-in reference set. 
# If the latter option is chosen, users can skip STEP 1 to 3 as described below. 

## SUMMARY: 
# STEP 1: Collect tumour fragment lengths
# STEP 2: Collect non-tumour fragment lengths
# STEP 3: Create fragment-level reference scores
# STEP 4: Determine Patient-level fragmentation score (FS)
# STEP 5: Determine Variant-level fragmentation score (VFS)


##############################################
################## STEP 1 ####################
##############################################
### COLLECT TUMOUR FRAGMENT LENGTHS

# Input required: VCF of tumour tissue data (germline variants and indels excluded). File name indicating corresponding plasma BAM file.

GetTumourReads <- function(patientID="patient_ID", VCF="Tissue_VCF", plasmaID="Plasma sample identifier") {
  require(Rsamtools)
  
  if (nrow(VCF) == 0) {print(paste0("No variants in ", patientID))}
  else {
    bamfile <- list.files(pattern="bam$")
    bamfileP <- bamfile[grep(plasmaID, bamfile)]
    indexx <- list.files(pattern="bai$")
    indexx <- indexx[grep(plasmaID, indexx)]
    
    if (length(bamfileP) + length(indexx) != 2) {print(paste0("BAM file ", plasmaID, " not found, proceeding without."))}
    else {
      bamfile <- BamFile(bamfileP, yieldSize=5e8)
      
      print(paste0("Processing ", plasmaID))
      ## @USER: Lines 36 and 46 refer to column names corresponding to Roche AVENIO VCFs. If the column names of your local pipeline are different, please change correspondingly
      ## @USER: Column V1 is the chromosome number, V2 is the chromosome position, V5 is the alternate base
      param <- ScanBamParam(which=GRanges(VCF$V1, IRanges(start=VCF$V2, end=VCF$V2)), what=scanBamWhat())
      bam <- scanBam(bamfile, param=param)
      
      ## Helper function to grab the lengths of the reads with tumour-derived mutations
      FragmentsWithVariants <- function(Chr, Pos, Alt, i, Bam) {
        base <- substring(as.character(Bam[[i]]$seq), Pos-Bam[[i]]$pos+1, Pos-Bam[[i]]$pos+1)
        lengths <- abs(Bam[[i]]$isize[base == Alt])
        lengths
      }
      
      lengths <- as.integer(unlist(mapply(FragmentsWithVariants, Chr=VCF$V1, Pos=as.numeric(VCF$V2), Alt=VCF$V5, i=1:nrow(VCF), MoreArgs=list(Bam=bam))))
      lengths <- table(factor(lengths, levels=1:900))
      return(lengths)
    }
  }
}

## @USER: Please write an appropriate 'wrapper' function that applies GetTumourReads() to all plasma samples in your cohort and calculates the rowSums() of the resulting tables. That way you end up with the total number of reads of each given length in the tumour reads set. 

##############################################
################## STEP 2 ####################
##############################################
### COLLECT NON-TUMOUR FRAGMENT LENGTHS

# Input required: File names indicating BAM file of healthy donor or patient with non-malignant disease. 

GetHealthyReads <- function(ID="healthy_ID") {
  require(Rsamtools)
  bamfileP <- list.files(pattern="bam$")[grep(ID, list.files(pattern="bam$"))]
  bamfile <- BamFile(bamfileP, yieldSize = 5e8)
  lengths <- scanBam(bamfile, param=ScanBamParam(what="isize"))
  lengths <- table(factor(abs(lengths[[1]]$isize), levels=1:900))
  return(lengths)
}

## @USER: Please write an appropriate 'wrapper' function that applies GetHealthyReads() to all healthy donor plasma samples and calculates the rowSums() as above. 

##############################################
################## STEP 3 ####################
##############################################
### GENERATE MUTANT:WILDTYPE LOG-RATIO REFERENCE DATASET

# Input required: Frequencies of fragment lengths in mutant and wildtype reads (collected in steps 1 and 2, respectively)

GenerateReferenceSet <- function(tlen="Tumour_lengths", hlen="Healthy_lengths") {
  scores <- NULL
  
  for (i in 1:1000) {
    hdlen <- table(factor(sample(as.numeric(names(hlen)), size=10000, prob=hlen, replace=TRUE), levels=1:900))/10000
    tumlen <- table(factor(sample(as.numeric(names(tlen)), size=10000, prob=tlen, replace=TRUE), levels=1:900))/10000
    
    tumlen[(tumlen + hdlen) < 20/10000] <- 0
    hdlen[(tumlen + hdlen) < 20/10000] <- 0
    
    score <- log(tumlen/hdlen, 2)
    score[is.na(score)] <- 0
    score[score == Inf] <- 5
    score[score == -Inf] <- -5
    
    scores <- cbind(scores, score)
  }
  score <- rowMeans(scores)
  return(score)
}

## @USER: Store the result.

##############################################
################## STEP 4 ####################
##############################################
### DETERMINE PATIENT-LEVEL FRAGMENTATION SCORES

# Input required: Reference dataset (from step 3 OR from pre-supplied data) and file names for plasma BAM files.

GeneratePatientFS <- function(ref, bamfile) {
  require(Rsamtools)
  
  bamfile <- bamfile
  indexx <- paste0(bamfile, ".bai")
  
  if (length(bamfile) + length(indexx) != 2) {print(paste0("BAM file ", bamfile, " not found, proceeding without."))}
  else {
    bamfile <- BamFile(bamfile)
    param <- ScanBamParam(what="isize",
                          flag = scanBamFlag(isPaired = TRUE,
                                             isProperPair = TRUE,
                                             isDuplicate = FALSE,
                                             isSecondaryAlignment = FALSE,
                                             isUnmappedQuery = FALSE),
                          mapqFilter = 30)
    reads <- scanBam(bamfile, param=param)
    lengths <- abs(reads[[1]]$isize)
    lengths[is.na(lengths)] <- 0
    lengths <- lengths[lengths < 900 & lengths > 0]
    PatientFS <- mean(ref[lengths])
    return(PatientFS)
  }
}

## @USER: Please write an appropriate 'wrapper' function to apply GeneratePatientFS() to all plasma samples in the cohort AND the healthy donors. 

##############################################
################## STEP 5 ####################
##############################################
### DETERMINE VARIANT-LEVEL FRAGMENTATION SCORES

# Input required: Reference dataset (from step 3 OR from pre-supplied data), list of VCFs for plasma samples of the same patient (excluding germline variants and indels), and string of file names for plasma BAM files.

GenerateVariantFS <- function(ref, vcf, bam_file, id) {
  require(dplyr)
  require(Rsamtools)
  
  ## @USER: Line 145 refers to column names corresponding to Roche AVENIO VCFs. If the column names of your local pipeline are different, please change correspondingly
  ## @USER: Column V1 is the chromosome number, V2 is the chromosome position, V5 is the alternate base
  variants <- unique(data.frame(CHR=vcf$CHROM, POS=vcf$POS, ALT=vcf$ALT, REF=vcf$REF)) # generates df of all unique variants ever seen in plasma of this patient
  
  if (nrow(variants) == 0) {print(paste0("No variants detected"))}
  else {
    all_lengths <- data.frame(NULL)
    
    bamfileP <- bam_file
    indexx <- paste0(bam_file, ".bai")
    
    if (length(bamfileP) + length(indexx) != 2) {print(paste0("BAM file:", bamfileP, "not found, proceeding without."))}
    else {
      bamfile <- BamFile(bamfileP)
      
      print(paste0("Processing ", bam_file))
      
      param <- ScanBamParam(what = scanBamWhat(),
                            which = GRanges(variants$CHR, 
                                            IRanges(start = variants$POS, 
                                                    end = variants$POS)),
                            flag = scanBamFlag(isPaired = TRUE,
                                               isProperPair = TRUE,
                                               isDuplicate = FALSE,
                                               isSecondaryAlignment = FALSE,
                                               isUnmappedQuery = FALSE),
                            mapqFilter = 30)
      bam <- scanBam(bamfile, param=param)
      
      ## Helper function to grab the lengths of the reads with tumour-derived mutations
      FragmentsWithVariants2 <- function(Chr, Pos, Alt, i, Bam) {
        base <- substring(as.character(Bam[[i]]$seq), Pos-Bam[[i]]$pos+1, Pos-Bam[[i]]$pos+1)
        lengths <- abs(Bam[[i]]$isize[base == Alt])
        res <- tibble(CHR=rep(Chr, length(lengths)), POS=rep(Pos, length(lengths)), ALT=rep(Alt, length(lengths)), LEN=lengths)
        return(res)
      }
      
      FragmentsWithoutVariants2 <- function(Chr, Pos, Ref, i, Bam) {
        base <- substring(as.character(Bam[[i]]$seq), Pos-Bam[[i]]$pos+1, Pos-Bam[[i]]$pos+1)
        lengths <- abs(Bam[[i]]$isize[base == Ref])
        res <- tibble(CHR=rep(Chr, length(lengths)), POS=rep(Pos, length(lengths)), REF=rep(Ref, length(lengths)), LEN=lengths)
        return(res)
      }
      
      var_lengths <- as.data.frame(t(mapply(FragmentsWithVariants2, 
                                            Chr=variants$CHR, 
                                            Pos=as.numeric(variants$POS), 
                                            Alt=variants$ALT,
                                            i=1:nrow(variants), 
                                            MoreArgs=list(Bam=bam))))
      var_lengths <- data.frame(cbind(CHR=unlist(var_lengths$CHR), 
                                      POS=as.numeric(unlist(var_lengths$POS)), 
                                      ALT=unlist(var_lengths$ALT), 
                                      LEN=as.numeric(unlist(var_lengths$LEN)),
                                      CNT=1))
      
      rownames(var_lengths) <- NULL
      var_lengths$LEN <- as.numeric(var_lengths$LEN)
      var_lengths$FS <- ref[var_lengths$LEN]
      
      all_lengths <- rbind(all_lengths, var_lengths)
      
      ### Repeat for wildtype fragments
      wt_lengths <- as.data.frame(t(mapply(FragmentsWithoutVariants2, 
                                           Chr=variants$CHR, 
                                           Pos=as.numeric(variants$POS), 
                                           Ref=variants$REF,
                                           i=1:nrow(variants), 
                                           MoreArgs=list(Bam=bam))))
      wt_lengths <- data.frame(cbind(CHR=unlist(wt_lengths$CHR), 
                                     POS=as.numeric(unlist(wt_lengths$POS)), 
                                     REF=unlist(wt_lengths$REF), 
                                     LEN=as.numeric(unlist(wt_lengths$LEN)),
                                     CNT=1))
      
      rownames(wt_lengths) <- NULL
      wt_lengths$LEN <- as.numeric(wt_lengths$LEN)
      wt_lengths$FS <- ref[wt_lengths$LEN]
    }
  }
  temp <- all_lengths %>%
    group_by(CHR, POS, ALT) %>%
    dplyr::summarize(VFS = mean(FS, na.rm = TRUE),
                     COUNT_VAR = n(),
                     LENGTH_VAR = mean(LEN, na.rm = TRUE),
                     var_s = sum((LEN - mean(LEN))^2))
  temp2 <- wt_lengths %>%
    group_by(CHR, POS, REF) %>%
    dplyr::summarize(WFS = mean(FS, na.rm = TRUE),
                     COUNT_WT = n(),
                     LENGTH_WT = mean(LEN, na.rm = TRUE),
                     wt_s = sum((LEN - mean(LEN))^2))
  
  var_scores <- merge(temp, temp2, by = c("CHR", "POS"))
  var_scores$s <- (var_scores$var_s + var_scores$wt_s)/(var_scores$COUNT_VAR + var_scores$COUNT_WT - 2)
  var_scores$t <- (var_scores$LENGTH_VAR - var_scores$LENGTH_WT)/sqrt(var_scores$s*((1/var_scores$COUNT_VAR) + (1/var_scores$COUNT_WT)))
  var_scores$t_test <- pt(var_scores$t, var_scores$COUNT_VAR + var_scores$COUNT_WT - 2, lower.tail = TRUE)
  ks_test <- c()
  for (i in c(1:nrow(var_scores))) {
    variant <- var_scores[i, c("CHR", "POS")]
    CHR <- variant$CHR
    POS <- variant$POS
    ks <- ks.test(all_lengths$LEN[all_lengths$CHR == CHR &
                                all_lengths$POS == POS],
                  wt_lengths$LEN[wt_lengths$CHR == CHR &
                               wt_lengths$POS == POS],
                  alternative = "less")$p.value
    ks_test <- c(ks_test, ks)
    }  
  var_scores$ks_test <- ks_test  
  var_scores <- var_scores[, c("CHR", "POS", "ALT", "REF", "COUNT_WT", "COUNT_VAR", "LENGTH_WT", "LENGTH_VAR", "WFS", "VFS", "t_test", "ks_test")]  
  return(var_scores)
}

## @USER: Please write an appropriate 'wrapper' function to apply GenerateVariantFS() to each patient (the function expects that the BAMs and VCFs from all plasma samples of that patient are provided at once)