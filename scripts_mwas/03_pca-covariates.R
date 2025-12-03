# CREATED: 2025-02-17
# PREG MATERNAL
# Principal Component Analysis and Covariate Selection

# https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_07_Functions_PCA.pdf



# LIBRARIES ----
library(minfi)
library(polycor)
library(tidyverse)
library(lme4)
library(ramwas)



# LOAD DATA OBJECTS ----
gset.quant <- readr::read_rds("data_objects/gset.quant.rds")

gset.temp <- gset.quant
ptemp <- as.data.frame(pData(gset.temp))
rm(gset.quant)



# LOAD AND PROCESS TIMEPOINT OBJECT ----
# Created locally in `05_sample-timepoints.R`; keeps all timepoints instead of those used for DNAm analysis
ga_timepoint_long_tbl <-
  readr::read_rds("data_objects/ga_timepoint_long_AllVisitsDates.rds") %>%
  mutate(timepoint = as.numeric(stringr::str_sub(visit, 2, 2))) %>%
  mutate(study_id = as.numeric(as.character(study_id))) %>% 
  dplyr::rename(timepoint_day = days)



# 1. ESTIMATE INITIAL PCA ----
# Not to be confused with the PCA of control probes in script `02_gather-DNAm.R`
# Compare initial PCA with residual PCA from covariate adjustment

# Non-Adjusted DNAm values
Mval <- getM(gset.temp)


Mval0 <- Mval - rowMeans(Mval)  #must be centered; scaling is optional
covmat <- crossprod(Mval0)                #covariance matrix
eig <- eigen(covmat)


# The number of principal components included as covariates is usually determined by the scree plot.
pdf("qc/PC_Mval-nonadjusted.pdf")
plotPCvalues(eig$values, n = 20)

plotPCvectors(eig$vectors[,1], i = 1, col = 'blue')
plotPCvectors(eig$vectors[,2], i = 2, col = 'red')
plotPCvectors(eig$vectors[,3], i = 3, col = 'green')
plotPCvectors(eig$vectors[,4], i = 4, col = 'purple')
plotPCvectors(eig$vectors[,5], i = 5, col = 'orange')
plotPCvectors(eig$vectors[,6], i = 6, col = 'brown')
plotPCvectors(eig$vectors[,7], i = 7, col = 'yellow')
plotPCvectors(eig$vectors[,8], i = 8, col = 'pink')

dev.off()

# (eig$values[1]*100) / sum(eig$values)
nPCs = 8
covariates.pca <- eig$vectors[, seq_len(nPCs)]
colnames(covariates.pca) <- paste0('PCdata',seq_len(nPCs))
# rm(Mval0, covmat, eig, nPCs)

ptemp <- bind_cols(ptemp, covariates.pca)
write_rds(ptemp, file = "data_objects/ptemp.rds")


# Use heatmeap to explore correlations of PCs with covariates
df1 <- ptemp %>%
  select(study_id, visit, Slide, row, 
         mMed, uMed, PC1, PC2, PC3, PC4, PC5,
         CD4T, CD8T, NK, Bcell, Mono, Gran,
         ga_ad, dm_16, mat_age_years, parous, preg_smoke, baby_a_gender)

pc1 <- ptemp %>%
  select(paste0("PCdata", 1:8))

# df1 <- as.data.frame(df1)
# pc1 <- as.data.frame(pc1)


out <- matrix(NA, (dim(df1)[2]*dim(pc1)[2]), 4)
count <- 0
for (i in 1:dim(df1)[2]) {
  for (j in 1:dim(pc1)[2]) {
    count <- count + 1
    a <- hetcor(df1[,i], pc1[,j])
    out[count,1] <- names(df1)[i]
    out[count,2] <- names(pc1)[j]
    out[count,3] <- round(a$corr[2,1], 4)
    out[count,4] <- round(a$std.errors[1, 2], 4)
  }
}

out <- as.data.frame(out)
names(out) <- c('Factor','PC','Corr', 'std_error')
#out$Factor <- as.character(out$Factor)
#out$PC <- as.character(out$PC)
# out$Corr <- as.character(out$Corr)
out$Corr <- as.numeric(out$Corr)
out$std_error <- as.numeric(out$std_error)

out$test_statistic <- out$Corr / out$std_error
out$p_value <- 2 * (1 - pnorm(abs(out$test_statistic)))


# Identify covariates significantly correlated with major PCs
out %>% 
  filter(PC == "PCdata8" & p_value < 0.05)


summary(lm(pc1$PCdata1 ~ df1$row + df1$mMed + df1$uMed + 
             df1$PC1 + df1$PC2 + df1$PC3 + df1$PC4 + df1$PC5 +
             df1$CD8T + df1$NK + df1$Bcell + df1$Mono + df1$Gran))



# Heatmap of experimental factors by PCs
hm <- ggplot(out, aes(PC, Factor )) +
  geom_tile(aes(fill = Corr), color = "white") +
  scale_fill_gradient(low = "yellow", high = "blue") +
  ylab("Experimental Factors") +
  xlab("PC") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Correlation")
hm

pdf("qc/PC-correlation-before-adjustment.pdf")
hm
dev.off()



# 2. Residualize signficant covariates ----
table(colnames(Mval) == rownames(ptemp))

ptemp.ortho <- ramwas::orthonormalizeCovariates(
  ptemp %>% 
  select(Slide, row, mMed, uMed, PC1, PC2, PC3, PC4, PC5, 
         CD8T, NK, Bcell, Mono, Gran,
         dm_16, parous, preg_smoke, mat_age_years),
  modelhasconstant = FALSE)



# Random selection of Mvals to test PCs of covariate adjusted mvals
iters <- 50000
mvals <- sample(1:dim(Mval)[1], iters, replace = FALSE)
rows <- list()
for (i in 1:iters) {     #dim(Mval)[1]) {
  z <- xpectr::suppress_mw(residuals(lmer(Mval[mvals[i],] ~ ptemp.ortho + (1|ptemp$study_id))))
  # z <- summary(lm(Mval[i,] ~ 1))$residuals
  rows[[i]] <- z
  if (i %% 1000 == 0) print(i)
}

Mval_adj_1 <- do.call(rbind, rows)
dim(Mval_adj_1)




# 3. ESTIMATE ADJUSTED PCA ----
# This code snippet has been added to `03_pca-data.R`
Mval1 <- Mval_adj_1 - rowMeans(Mval_adj_1)  #must be centered; scaling is optional
covmat <- crossprod(Mval1)                #covariance matrix
eig <- eigen(covmat)

pdf("qc/PC_Mval-adjusted.pdf")
plotPCvalues(eig$values, n = 20)

plotPCvectors(eig$vectors[,1], i = 1, col = 'blue')
plotPCvectors(eig$vectors[,2], i = 2, col = 'red')
plotPCvectors(eig$vectors[,3], i = 3, col = 'green')
plotPCvectors(eig$vectors[,4], i = 4, col = 'purple')
plotPCvectors(eig$vectors[,5], i = 5, col = 'orange')
plotPCvectors(eig$vectors[,6], i = 6, col = 'brown')
plotPCvectors(eig$vectors[,7], i = 7, col = 'yellow')
plotPCvectors(eig$vectors[,8], i = 8, col = 'pink')

dev.off()


# ADD TIMEPOINT_DAY TO PTEMP OBJECT ----
ptemp <- ptemp %>% 
  mutate(timepoint = visit,
         array_name = rownames(ptemp))


pheno_tbl <- 
  left_join(ptemp, ga_timepoint_long_tbl, by = c("study_id", "timepoint")) %>%
    mutate(timepoint_day = as.numeric(timepoint_day))

dim(pheno_tbl)  #646 46
table(is.na(pheno_tbl$timepoint_day))
table(colnames(gset.temp) == rownames(ptemp))
table(colnames(gset.temp) == pheno_tbl$array_name)



# OUTPUT FINAL GSET OBJECT ----
gset.quant.2 <- gset.temp

table(pheno_tbl$array_name == colnames(gset.quant.2))
pData(gset.quant.2) <- as(pheno_tbl, "DataFrame")

write_rds(gset.quant.2, file = "data_objects/gset.quant.2.rds")
# gset.quant.2 <- readr::read_rds("data_objects/gset.quant.2.rds")



