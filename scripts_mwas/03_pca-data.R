# CREATED: 2025-03-12
# PREG MATERNAL
# Data PCs

# Estimate and output data PCs from covariate adjusted mvals
# These data PCs will be added to pData in '03_pca-covariates.R`



# LIBRARIES ----
library(minfi)
library(tidyverse)
library(lme4)
library(tictoc)



# LOAD DATA FROM `03_pca-covariates_v2.R` ----
load("/lustre/home/tpyork/projects/preg-maternal/data_objects/residuals-input.Rdata")



# PROCESS DATA ----
Mval <- getM(gset.temp)

# Verify missing rates
# skimr::skim(
#   ptemp %>% 
#     select(study_id, ga_ad,
#            mMed, uMed, row,
#            processed_date, visit,
#            PCtech1, PCtech2, PCtech3, PCtech4, PCtech5,
#            PCanc1, PCanc2, PCanc3, PCanc4, PCanc5,
#            CD8T, CD4T, NK, Bcell, Mono, Gran,
#            mat_age_years, preg_smoke, dm_16))



# ORTHONORMALIZE COVARIATES ----
# Done for covariates that absolutely will be considered
ptemp.ortho <- ramwas::orthonormalizeCovariates(
  ptemp %>% 
    select(processed_date, row, mMed, uMed, 
           PCtech1, PCtech2, PCtech3, PCtech4, PCtech5, 
           PCanc3,
           CD8T, NK, Bcell, Mono, Gran,
           preg_smoke, mat_age_years),
  modelhasconstant = FALSE)



# RESIDUALIZE MVALS ----
# Regress out measured covariates that will be used in analysis

# Random selection of Mvals to test PCs of covariate adjusted mvals
iters <- dim(Mval)[1]
# mvals <- sample(1:dim(Mval)[1], iters, replace = FALSE)
rows <- list()
tic("residualize")
for (i in 1:iters) {     
  z <- xpectr::suppress_mw(residuals(lmer(Mval[i,] ~ ptemp.ortho + (1|ptemp$study_id))))
  # z <- summary(lm(Mval[i,] ~ 1))$residuals
  rows[[i]] <- z
  if (i %% 1000 == 0) print(i)
}
toc()

Mval_adj_1 <- do.call(rbind, rows)
dim(Mval_adj_1)



# OUTPUT ----
save(Mval_adj_1, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects/residuals-output.Rdata")

