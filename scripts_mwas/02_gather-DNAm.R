# CREATED: 2025-02-05
# PREG MATERNAL
# Read in and process DNAm array data



# LIBRARIES ----
library("minfi")
library("ramwas")
library("tidyverse")
# library("FlowSorted.Blood.EPIC")
# library("AnnotationHub") #for EPICv2 manifest (Peters 2024)


# For ancestry estimation.
#select which SNP set to use at "tpyork/r-functions/DNAm-pc-pop-strat/
#needs to be made into a package
githubURL <- "https://raw.githubusercontent.com/tpyork/r-functions/master/DNAm-pc-pop-strat/cpgs_with_TGP_SNP_ON_PROBE.Rdata"
load(url(githubURL))

#SNP check function
devtools::source_url("https://raw.githubusercontent.com/tpyork/r-functions/master/DNAm_450k_FXNs.R")




# LOAD RGSET ----
# Created locally (RGset_temp) in local `scripts/02_gather-DNAm`
load("data_objects/RGset_temp.Rdata")   #622399    655



# INIITIAL QC OF RAW DATA ----
qcReport(RGset, pdf= "qc/quality-control-1.pdf") 

Mset <- preprocessRaw(RGset)
qc <- getQC(Mset)

detP <- detectionP(RGset)       #identifies failed positions
failed <- detP > 0.01

pdf("qc/quality-control-2.pdf")
plotQC(qc)                       #save this separate
plot(colMeans(failed), main="Percentage of failed probes")
dev.off()


# Identify poor quality sample
x <-  cbind(1:dim(RGset)[2], as.data.frame(colData(RGset)[c(1,7)]), as.data.frame(qc))
names(x)[1] <- 'num'

omit.sample.name <- x[(x$mMed < 10.756) & (x$uMed < 10.5), 2]
omit.num <- which(pData(RGset)$Sample_Name %in% omit.sample.name)


# Reassess QC on reduced sample
RGset2 <- RGset[, -omit.num]
Mset2 <- preprocessRaw(RGset2)    #convert R/G channels to methylated and unmethylated signals
qc2 <- getQC(Mset2)
detP <- detectionP(RGset2)  #identifies failed positions
failed <- detP > 0.01


p1 <- pData(RGset) %>% 
  as.data.frame() %>% 
  mutate(sid = str_sub(Sample_Name, 1, 4))
length(unique(p1$sid))  #174

p2 <- pData(RGset2) %>% 
  as.data.frame() %>% 
  mutate(sid = str_sub(Sample_Name, 1, 4))
length(unique(p2$sid))




# updated qc with subjects removed
pdf("qc/qcReport_2.5.pdf")
plotQC(qc2)                      #save this separate
densityPlot(RGset2, main= "Beta", xlab= "Beta")   #save this separate
# densityBeanPlot(RGset2, sampGroups= pData(RGset)$timepoint)                           #save this separate
plot(colMeans(failed), main="Percentage of failed probes")
dev.off()


rm(Mset, failed, detP, failed, qc)



# FILTER PROBES AND SAMPLES ----

## Filter cross-reactive probes ----
exclude.xreactive <- maxprobes::xreactive_probes(array_type = "450K")
length(exclude.xreactive)


## Filter probes with low beadcount ----
# Needs RGChannelSetExtended; obtain by setting extended = TRUE in read.metharray.exp
lb <- getNBeads(RGset) < 3
pi1 <- getProbeInfo(RGset, type = "I")
pi2 <- getProbeInfo(RGset, type = "II")
ex1 <- pi1$Name[rowMeans(lb[pi1$AddressA,] | lb[pi1$AddressB,]) > 0.01]
ex2 <- pi2$Name[rowMeans(lb[pi2$AddressA,]) > 0.01]
exclude.bds <- unique(c(ex1, ex2))
rm(lb, pi1, pi2, ex1, ex2)



## Filter SNPs ----
# Identify probes with SNPs in critical areas
#Probe_SNPs:  Assays with SNPs present within probe > 10bp from query site
#Probe_SNPs_10:  Assays with SNPs present within probe <= 10bp from query site
#Probe_rs:  RS number for a SNP overlapping the probe
#Probe_maf:  minor allele freq; if multiple SNPs overlap only one is recorded
#CpG_rs:  SNPs overlapping the CpG site
#SBE_rs:  the single base extension of the measured methylation loci

# ann <- getAnnotation(RGset)
snps <- getSnpInfo(RGset)  #representation of db.snp
hyb.SNPs <- snps[!is.na(snps$CpG_rs) | !is.na(snps$SBE_rs), ] # & is.na(snps$Probe_rs)
exclude.snp <- rownames(hyb.SNPs)
rm(snps, hyb.SNPs)



## Filter probes and samples with low detection p-values ----
# hp <- detectionP(RGset) > 0.01
hp <- detP > 0.01
sum(rowMeans(hp) > 0.01)  #number of failed positions in > 10% of samples; 5227
exclude.hpv <- rownames(hp)[rowMeans(hp) > 0.01]
keep.samples <- colMeans(hp) < 0.01
rm(hp)
# table(keep.samples)  #N = 9



## Filter low quality samples and probes ----
# REPLACES CHECKPOINT 2
# Remove `exclude.snp` later
RGset3 <- subsetByLoci(rgSet       = RGset2[, keep.samples],
                       excludeLoci = c(exclude.xreactive, exclude.bds, exclude.hpv))
save(RGset3, file = "data_objects/RGset_checkpoint3.Rdata")

p3 <- pData(RGset3) %>% 
  as.data.frame() %>% 
  mutate(sid = str_sub(Sample_Name, 1, 4))
length(unique(p3$sid))




# CHECKPOINT 2 ----
# load("data_objects/RGset_checkpoint2.Rdata")





# ESTIMATE ANCESTRY PCs -------------------------------------------------
#done here to verify self-reported race
#will exclude.reactive later
gset.quant <- preprocessQuantile(RGset3)

# Estimate before probes with SNPs are screened
temp2 <- gset.quant
Bval <- getBeta(temp2)  #can optionally est PCs based on Mval
dim(Bval)


location.based.pc2 <- function(beta.obj, a){
  #load(file.pathname)
  a2 <- a
  
  beta.obj<-as.matrix(beta.obj)
  #gc()
  if(sum(is.na(beta.obj))>0){
    cpgmeans<-t(t(rowMeans(beta.obj,na.rm=T))) %*% matrix(1,ncol=ncol(beta.obj))
    #gc()
    missing.sites<-which(is.na(beta.obj))
    #gc()
    beta.obj[missing.sites]<-cpgmeans[missing.sites]
    rm(cpgmeans)
    #gc()
  }
  beta.obj<-beta.obj[which(rownames(beta.obj)%in%a2[,1]),]
  theresult<-princomp(x=beta.obj,cor=T)
  rm(a2)
  gc()
  return(theresult)
}

pc <- location.based.pc2(Bval, a2)


# Add PCs to gset.quant
table(row.names(pData(temp2))==row.names(pc$loadings))

head(pc$loadings[,1:10])
head(pData(temp2))

pData(temp2)$Comp.Anc.1 <- pc$loadings[,1]
pData(temp2)$Comp.Anc.2 <- pc$loadings[,2]
pData(temp2)$Comp.Anc.3 <- pc$loadings[,3]
pData(temp2)$Comp.Anc.4 <- pc$loadings[,4]
pData(temp2)$Comp.Anc.5 <- pc$loadings[,5]
pData(temp2)$Comp.Anc.6 <- pc$loadings[,6]
pData(temp2)$Comp.Anc.7 <- pc$loadings[,7]
pData(temp2)$Comp.Anc.8 <- pc$loadings[,8]
pData(temp2)$Comp.Anc.9 <- pc$loadings[,9]
pData(temp2)$Comp.Anc.10 <- pc$loadings[,10]


# plot separation by PCs
temp <- as.data.frame(pData(temp2))
g <- ggplot(temp, aes(x=Comp.Anc.3, y=Comp.Anc.2)) + geom_point(aes(colour = as.factor(temp$dm_16)))


library(scatterplot3d)
colors <- ifelse(temp2$dm_16==1, 'black', 'red')
pdf("qc/PC-ancestry.pdf")
scatterplot3d(temp2$Comp.Anc.3, temp2$Comp.Anc.1, temp2$ga_ad, main="Ancestry", color= colors)
g
dev.off()


summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.3)) #*
summary(lm(temp2$ga_ad ~ temp2$dm_16))
summary(lm(temp2$ga_ad ~ temp2$dm_16 + temp2$Comp.Anc.3))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.1)) #*
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.2))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.4))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.5))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.6))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.7))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.8))  #*
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.9))
summary(lm(temp2$ga_ad ~ temp2$Comp.Anc.10)) #*


summary(lm(temp2$Comp.Anc.1 ~ temp2$dm_16)) #*** R2 = 0.0664
summary(lm(temp2$Comp.Anc.2 ~ temp2$dm_16)) #*** R2 = 0.0291
summary(lm(temp2$Comp.Anc.3 ~ temp2$dm_16)) #*** R2 = 0.8444
summary(lm(temp2$Comp.Anc.4 ~ temp2$dm_16)) #*** R2 = 0.0731
summary(lm(temp2$Comp.Anc.5 ~ temp2$dm_16))
summary(lm(temp2$Comp.Anc.6 ~ temp2$dm_16))
summary(lm(temp2$Comp.Anc.7 ~ temp2$dm_16))
summary(lm(temp2$Comp.Anc.8 ~ temp2$dm_16))
summary(lm(temp2$Comp.Anc.9 ~ temp2$dm_16))
summary(lm(temp2$Comp.Anc.10 ~ temp2$dm_16))


#update Race for 7185
#pheno3 <- pData(temp2)
#pheno3[pheno3$Comp.Anc.2 < -.05 & pheno3$dm_16==1, ]  #find which study_id to update
#pData(temp2)[pData(temp2)$study_id==7185, 'dm_16'] <- 0


#temp <- as.data.frame(pData(temp2)[c("dm_16", "Comp.Anc.1", "Comp.Anc.2", "Comp.Anc.3")])
#g <- ggplot(temp, aes(x=Comp.Anc.3, y=Comp.Anc.2)) + geom_point(color= as.factor(temp$dm_16+1))
#g

#colors <- ifelse(temp2$dm_16==1, 'black', 'red')
#pdf("qc/PC-ancestry2.pdf")
#scatterplot3d(temp2$Comp.Anc.2, temp2$Comp.Anc.1, temp2$Comp.Anc.3, main="Ancestry", color= colors)
#g
#dev.off()

# Store ancestry PCs
readr::write_rds(pc, file = "data_objects/ancestry-pc.rds")
rm(temp2, Bval)



# SNP CHECK AND REMOVE MISLABELED SAMPLES ----
# 7115_EK_V4, `predicted_sex` == 'M'
RGset3$study_id <- str_sub(RGset3$Sample_Name, 1, 4)

out <- snpCheckMySamples(RGset_object       = RGset3,
                         twins              = FALSE,
                         repeated_measure   = TRUE,
                         common_ID          = "study_id",
                         individual_ID      = "Sample_Name",
                         highMatchThreshold = .80,
                         lowMatchThreshold  = .65)

dim(out$highMatch)
dim(out$lowMatch)

out$lowMatch    #falls below `highMatchThreshold`
out$highMatch   #falls above `lowMatchThreshold`


# Remove mislabeled samples
exclude.mislabel <- 
  !(pData(RGset3)$Sample_Name %in%
  c('7031_BE_V1', '7080_DB_V2', '7105_EA_V2', '7106_EB_V2', '7109_EE_V3', '7115_EK_V4'))
table(exclude.mislabel)  

RGset4 <- RGset3[, exclude.mislabel]



# Obtain methylation estimates for QC ----
# used for minfi::getQC() below; inspect QC reports after minimal processing
# compare to QC report after normalization
rgSetRaw <- fixMethOutliers(preprocessRaw(RGset4))  #no normalization

# Will normalize data below
# beta = BMIQ(rgSetRaw)
# beta <- getBeta(rgSetRaw)  # needed only for ramwas format


# Principal components analysis (PCA) on control probes ----
# Uses cpg-filtered data

# We extract red and green channel for the control probes.
controlType <- unique(getManifest(RGset4)@data$TypeControl$Type)
controlSet <- getControlAddress(RGset4, controlType = controlType)
probeData <- rbind(getRed(RGset4)[controlSet,], getGreen(RGset4)[controlSet,])

# Next we run principal component analysis on the data after light normalization.
data <- probeData - rowMeans(probeData)  #must be centered; scaling is optional
covmat <- crossprod(data)                #covariance matrix
eig <- eigen(covmat)

# The number of principal components included as covariates is usually determined by the scree plot.
pdf("qc/PC_controlprobes.pdf")
plotPCvalues(eig$values, n = 20)

plotPCvectors(eig$vectors[,1], i = 1, col = 'blue')
plotPCvectors(eig$vectors[,2], i = 2, col = 'red')
plotPCvectors(eig$vectors[,3], i = 3, col = 'green')
plotPCvectors(eig$vectors[,4], i = 4, col = 'purple')
plotPCvectors(eig$vectors[,5], i = 5, col = 'orange')
plotPCvectors(eig$vectors[,6], i = 6, col = 'brown')
dev.off()

# (eig$values[1]*100) / sum(eig$values)

nPCs = 5
covariates.pca <- eig$vectors[, seq_len(nPCs)]
colnames(covariates.pca) <- paste0('PCtech',seq_len(nPCs))
rm(probeData, data, covmat, eig, nPCs) 



# Cell type composition ----

# For 450k and EPICv1 arrays
covariates.cel <- estimateCellCounts(
  rgSet             = RGset4,
  compositeCellType = "Blood",
  referencePlatform = "IlluminaHumanMethylation450k",
  cellTypes         = c("CD8T", "CD4T", "NK",
                        "Bcell", "Mono", "Gran"),
  meanPlot          = TRUE,
  verbose           = TRUE)


# For EPICv2 arrays
# https://support.bioconductor.org/p/9156892/
# MSet <- preprocessNoob(RGset2)
# 
# Betas <- getBeta(MSet)
# Betas <- sesame::betasCollapseToPfx(Betas) #you can also use ENmix::rm.cgsuffix(Betas) or other function to remove replicates 
# 
# IDOLOptimizedCpGsBloodv2 <- IDOLOptimizedCpGs[which(IDOLOptimizedCpGs%in%rownames(Betas))]
# identical(rownames(IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,]), IDOLOptimizedCpGsBloodv2)
# covariates.cel <- projectCellType_CP(
#   Betas[IDOLOptimizedCpGsBloodv2, ],
#   IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,],
#   contrastWBC = NULL, nonnegative = TRUE,
#   lessThanOne = FALSE
# )



# Median methylated and unmethylated intensity ----
covariates.umm = getQC(rgSetRaw)



# Phenotypic covariates from the sample sheet ----
# Separate for ramwas; but add to RGset2 for other
covariates.phe <- pData(RGset4)


covariates.anc <- as.data.frame(pc$loadings[,1:5])
names(covariates.anc) <- paste0("PCanc", 1:5)
covariates.anc$samples <- rownames(covariates.anc)

covariates <- data.frame(
  samples = rownames(covariates.umm),
  covariates.umm,
  covariates.pca,
  covariates.cel,
  covariates.phe)

covariates <- left_join(covariates, covariates.anc, by = "samples")
# skimr::skim(covariates)

# CHECKPOINT 4 ----
save(RGset4, covariates, file = "data_objects/checkpoint4.Rdata")
# load("data_objects/checkpoint4.Rdata")

p4 <- pData(RGset4) %>% 
  as.data.frame() %>% 
  mutate(sid = str_sub(Sample_Name, 1, 4))
length(unique(p4$sid))



# Get study_id and visit
covariates$visit <- as.numeric(substr(covariates$Sample_Name, 10, 10))
# covariates$study_id <- as.numeric(substr(covariates$Sample_Name, 1, 4))

# Get array row
covariates$row <- as.numeric(substr(covariates$Array, 3, 3))

# Force slide as categorical
covariates$Slide <- paste0("x", covariates$Slide)

# Finalize covariates
covariates <- covariates %>% 
  relocate(c(study_id, Sample_Name), .before = samples) %>% 
  relocate(c(Basename, filenames), .after = last_col())


# table(row.names(covariates) == covariates$samples)
write_rds(covariates, file = "data_objects/covariates.rds")
# covariates <- readr::read_rds(file = "data_objects/covariates.rds")


# Add all covariates to pData
# table(rownames(covariates) == rownames(pData(RGset2)))
pData(RGset4) <- as(covariates, "DataFrame")



# NORMALIZE ----
gset.quant <- preprocessQuantile(RGset4)
# table(rownames(covariates) == rownames(pData(gset.quant)))


# Remove chrY, chrX probes
gset.quant <- gset.quant[seqnames(gset.quant) != 'chrY']  #removes 3 chrY probes
# temp.quant <- temp.quant[seqnames(temp.quant) != 'chrX']  #NOT RUN; would remove 10397 chrX probes



# CHECKPOINT 4 ----
readr::write_rds(gset.quant, file = "data_objects/gset.quant.4.rds")
# gset.quant <- read_rds("data_objects/gset.quant.rds")
# covariates <- read_rds("data_objects/covariates.rds")



# QC COMPARISON ----
pdf("qc/quality-control-3.pdf")
densityPlot(RGset4, main= "Raw Data", xlab= "Beta")   #save this separate
densityPlot(getBeta(gset.quant), main= "Quant Normalized", xlab= "Beta")
dev.off()







