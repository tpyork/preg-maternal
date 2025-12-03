# CREATED: 2025-02-11
# PREG MATERNAL
# RaMWAS Prediction of GAD


## PREDICT GA_AD WITH AND WITHOUT PC AND COVARIATE ADJUSTMENT FOR VISIT ONE.  INSPECT LAMBDA.

# library(foreach)
# library(doParallel)
# 
# cores <- detectCores()
# 
# cl <- makeCluster(cores-1)
# 
# registerDoParallel(cl)
# 
# n <- 1000
# 
# big_list <- list()
# 
# big_list <- foreach(i = 1:n) %dopar% {
#   biglist[i] <- i^2
# }
# 
# stopCluster(cl)





# LIBRARIES ----
library("minfi")
library("ramwas")
library("tidyverse")
# library("FlowSorted.Blood.EPIC")
# library("AnnotationHub") #for EPICv2 manifest (Peters 2024)



# LOAD DATA OBJECTS ----
load("data_objects/RGset_checkpoint2.Rdata")  #RGset2
gset.quant <- readr::read_rds("data_objects/gset.quant.2.rds")
# covariates <- read_rds("data_objects/covariates.rds")



# SUBSET DATA ----
# Visit 1 only
gset.quant.v4 <- gset.quant[, gset.quant$visit == 4]
covariates <- as.data.frame(pData(gset.quant.v4))

covariates$SampleID <- rownames(covariates)
covariates <- covariates |>
  dplyr::relocate(SampleID, .before = study_id)




# SAVE IN RAMWAS FORMAT ----
dir.create('rw_v4', showWarnings = FALSE)

# NORMALIZE DATA FIRST HERE OR WHEN TESTING NORMALIZED QC
# beta2 = getBeta(rgSetRaw)
beta <- getBeta(gset.quant.v4)

rng = granges(mapToGenome(RGset2))
chr = seqnames(rng)

# Save CpG locations
library(filematrix)
locs = cbind(chr = as.integer(chr), position = start(rng))
fmloc = fm.create.from.matrix(
  filenamebase = paste0("rw_v4/CpG_locations"),
  mat = locs,
  size = 4)
close(fmloc)
writeLines(con = 'rw_v4/CpG_chromosome_names.txt', text = levels(chr))

# Save estimates
fm = fm.create.from.matrix(
  filenamebase = paste0("rw_v4/Coverage"),
  mat = t(beta))   #normalized beta values
close(fm)


# RUN RAMWAS ----

## 1. Set Up Parameters ----
param <- ramwasParameters(
  dircoveragenorm = 'rw_v4',
  covariates      = covariates,
  modelcovariates = NULL,
  modeloutcome    = "ga_ad",
  modelPCs        = 0,
  cvnfolds        = 5)


## 2. Covariate Pruning ----


## 3. MWAS Without Covariates ----
param$modelcovariates <- NULL
param$modelPCs <- 0
# ramwas4PCA(param)
ramwas5MWAS(param)
qqPlotFast(getMWAS(param)$`p-value`)
title('QQ-plot\nNo covariates, no PCs')




## 4. MWAS With Covariates ----
param$modelcovariates <- c(
  "mat_age_years", "preg_smoke", "dm_16", "parous",
  "row", "Slide",
  "mMed", "uMed",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
param$modelPCs <- 0
# ramwas4PCA(param)
ramwas5MWAS(param)
qqPlotFast(getMWAS(param)$`p-value`)
title('QQ-plot\nwith covariates and 0 PC')



## 5. MWAS With Covariates and PCs ----
param$modelcovariates <- c(
  "mat_age_years", "preg_smoke", "dm_16", "parous",
  "row", "Slide",
  "mMed", "uMed",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran",
  "PCdata1", "PCdata2")
param$modelPCs <- 0
# ramwas4PCA(param)
ramwas5MWAS(param)
qqPlotFast(getMWAS(param)$`p-value`)
title('QQ-plot\nwith covariates and 2 PC')




# ramwas6annotateTopFindings(param)
# ramwas7riskScoreCV(param)







