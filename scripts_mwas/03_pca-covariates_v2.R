# CREATED: 2025-02-17
# PREG MATERNAL
# Principal Component Analysis and Covariate Selection

# https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_07_Functions_PCA.pdf

# USE THIS: _v2 updates the data PCs estimated on all residualized mvals



# LIBRARIES ----
library(minfi)
library(polycor)
library(tidyverse)
library(lme4)
library(ramwas)



# LOAD DATA OBJECTS ----
gset.quant <- readr::read_rds("data_objects/gset.quant.4.rds")

gset.temp <- gset.quant
ptemp <- as_tibble(pData(gset.temp)) %>% 
  mutate(study_id = as.factor(study_id))
rm(gset.quant)



# LOAD DATA FOR PROCESSED.DATE ----
# Processed_date
load("data_objects/so.Rdata")

so_tbl <- so_tbl %>%
  mutate(study_id       = as.factor(study_id),
         processed_date = as.factor(processed_date),
         visit          = timepoint) %>%
  select(study_id, visit, processed_date)


ptemp <- left_join(ptemp, so_tbl, by= c('study_id','visit')) %>% 
  select(study_id, ga_ad,
         mMed, uMed, row,
         processed_date, visit,
         PCtech1, PCtech2, PCtech3, PCtech4, PCtech5,
         PCanc1, PCanc2, PCanc3, PCanc4, PCanc5,
         CD8T, CD4T, NK, Bcell, Mono, Gran,
         mat_age_years, preg_smoke, dm_16) 



# LOAD AND PROCESS TIMEPOINT OBJECT ----
# Created locally in `05_sample-timepoints.R`; keeps all timepoints instead of those used just for DNAm analysis
ga_timepoint_long_tbl <-
  readr::read_rds("data_objects/ga_timepoint_long_AllVisitsDates.rds") %>%
  mutate(timepoint = as.numeric(stringr::str_sub(visit, 2, 2))) %>%
  # mutate(study_id = as.numeric(as.character(study_id))) %>% 
  dplyr::rename(timepoint_day = days) %>% 
  select(-visit) %>% 
  dplyr::rename(visit = timepoint)

ptemp <- 
  left_join(ptemp, ga_timepoint_long_tbl, by = c('study_id', 'visit'))

rm(ga_timepoint_long_tbl, so_tbl)




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
# nPCs = 10
# covariates.pca <- eig$vectors[, seq_len(nPCs)]
# colnames(covariates.pca) <- paste0('PCdata',seq_len(nPCs))
# rm(Mval0, covmat, eig, nPCs)
# ptemp2 <- bind_cols(ptemp, covariates.pca)
# write_rds(ptemp, file = "data_objects/ptemp.rds")




# Use heatmeap to explore correlations of PCs with covariates
df1 <- ptemp2 %>%
  select(ga_ad,
         mMed, uMed, row,
         processed_date, visit,
         PCtech1, PCtech2, PCtech3, PCtech4, PCtech5,
         PCanc1, PCanc2, PCanc3, PCanc4, PCanc5,
         CD8T, CD4T, NK, Bcell, Mono, Gran,
         mat_age_years, preg_smoke, dm_16)

# pc1 <- ptemp %>%
#   select(paste0("PCdata", 1:8))
pc1 <- covariates.pca
df1 <- as.data.frame(df1)
pc1 <- as.data.frame(pc1)


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
  filter(PC == "PCdata1" & p_value < 0.05)

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

rm(hm, out, pc1, a, df1, covariates.pca, Mval)



# 2. Residualize mvals by measured covariates ----

# This code run in `03_pca-data.R` to get residualized mvals
table(ptemp$study_id == pData(gset.temp)$study_id)

#Send data for `03_pca-data.R`
save(gset.temp, ptemp, file = "data_objects/residuals-input.Rdata")

#Get data PCs to add to ptemp
#Get residualized mvals to plot adjusted PCs
load("data_objects/residuals-input.Rdata")
load("data_objects/residuals-output.Rdata")



# 3. ESTIMATE ADJUSTED PCA ----
# This code snippet has been added to `03_pca-data.R`
Mval1 <- Mval_adj_1 - rowMeans(Mval_adj_1)  #must be centered; scaling is optional
covmat <- crossprod(Mval1)                #covariance matrix
eig <- eigen(covmat)

pdf("qc/PC_Mval-adjusted-2.pdf")
plotPCvalues(eig$values, n = 30)

plotPCvectors(eig$vectors[,1], i = 1, col = 'blue')
plotPCvectors(eig$vectors[,2], i = 2, col = 'red')
plotPCvectors(eig$vectors[,3], i = 3, col = 'green')
plotPCvectors(eig$vectors[,4], i = 4, col = 'purple')
plotPCvectors(eig$vectors[,5], i = 5, col = 'orange')
plotPCvectors(eig$vectors[,6], i = 6, col = 'brown')
plotPCvectors(eig$vectors[,7], i = 7, col = 'yellow')
plotPCvectors(eig$vectors[,8], i = 8, col = 'pink')

dev.off()

# Add data PCs to ptemp
nPCs = 20
covariates.pca <- eig$vectors[, seq_len(nPCs)]
colnames(covariates.pca) <- paste0('PCdata',seq_len(nPCs))
# rm(Mval0, covmat, eig, nPCs)

ptemp <- bind_cols(ptemp, covariates.pca)


# Get sample loadings
eignenvalues <- eig$values
eigenvectors <- eig$vectors
loadings <- eigenvectors %*% diag(sqrt(eignenvalues))

loadings_df <- data.frame(id  = paste(ptemp$study_id, ptemp$visit, sep = "_"),
                          PC1outlier = loadings[,1], PC2outlier = loadings[,2], 
                          PC3outlier = loadings[,3], PC4outlier = loadings[,4], 
                          PC5outlier = loadings[,5], PC6outlier = loadings[,6], 
                          PC7outlier = loadings[,7], PC8outlier = loadings[,8], 
                          PC9outlier = loadings[,9], PC10outlier = loadings[,10])

# Plot the loadings using ggplot2
# ggplot(loadings_df, aes(x = PC1, y = PC2, label = id)) +
#   geom_point() +
#   geom_text(vjust = -0.5) +
#   labs(title = "Loadings Plot", x = "Principal Component 1", y = "Principal Component 2") +
#   theme_minimal()


z_scores <- scale(loadings[,20])
temp_tbl <- tibble(sample = 1:625, scaled_loading = z_scores)
ggplot(temp_tbl, aes(x = sample, y = scaled_loading)) +
  geom_point() +
  theme_minimal()

# Identify outliers
z_scores <- scale(loadings[, 3])
outliers <- which(abs(z_scores) > 3)
outliers


# Add scaled loadings to ptemp
temp_tbl <- loadings_df %>% 
  mutate(across(starts_with("PC"), scale))

which(abs(temp_tbl %>% select(PC3outlier)) > 3)

ptemp <- cbind(ptemp, temp_tbl %>% select(-id))




# OUTPUT FINAL GSET OBJECT ----
gset.quant.5 <- gset.temp

pData(gset.quant.5) <- as(ptemp, "DataFrame")

write_rds(gset.quant.5, file = "data_objects/gset.quant.5.rds")



