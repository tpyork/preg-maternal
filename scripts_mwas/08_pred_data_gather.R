# CREATED: 2025-05-01
# Prediction Data Gather


# OUTLINE ----

# Predict GAD using v1 methylation levels
# Select CpG based on whether their levels change across pregnancy
# Vary CpG inclusion by 1%, 5%, 10% FDR
# Consider robustness of selectin CpGs based on change by permuting with randomly selected CpGs
# Compare to existing clinical risk predictors; maternal age, parity, smoking
# Combine CpG and clinical risk predictors
# Test using V1, V2, V3 CpG levels
# Try XGBoost; Ridge; LASSO; MARS; ensemble
#
# OTHER:
# Add back removed subjects based on PC outliers and reidentify CpGs



# SET UP WORKSPACE AND LOAD PACKAGES ----
# DNAm object
library(minfi)

# Core
library(tidyverse)


# GATHER DATA ----

#* GATHER CPG CHANGE RESULTS ----
gout <- readr::read_rds("data_objects/Gout-25cov.rds")

cpg_keep <- gout[gout$mix.qval.bac < 0.01, ] %>% 
  as.data.frame() %>% rownames()
cpg_keep %>% length()


#* GATHER CPG SUBJECT TIMEPOINT DATA ----
gset.quant <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects/gset.quant.578.rds")
# length(unique(pData(gset.quant)$study_id)) # 171


cpg_id <- granges(gset.quant) %>% 
  names()

mval <- getM(gset.quant) %>% 
  as.data.frame() %>% 
  t()

names(mval) <- cpg_id

mval <- mval %>% 
  as_tibble()



#* GATHER PHENOTYPIC DATA ----
pdata <- as_tibble(pData(gset.quant)) %>%
  select(study_id, ga_ad, timepoint_day, visit,
         preg_smoke, mat_age_years, dm_16,
         CD8T, NK, Bcell, Mono, Gran) %>% 
  mutate(study_id   = as.factor(study_id),
         preg_smoke = as.factor(preg_smoke),
         dm_16      = as.factor(dm_16))

# table(pdata$visit)



#* CREATE DATA FOR MODEL ----
mval <- mval[, cpg_keep]

preg_dat <- bind_cols(pdata, mval)
# length(unique(preg_dat$study_id)) #171

preg_dat <- preg_dat %>% 
  filter(visit == 1)
# length(unique(preg_dat$study_id)) #151


# Clean up
rm(gset.quant, gout, mval, pdata, cpg_id, cpg_keep)



# OUTPUT DATA ----
readr::write_rds(preg_dat, file = "data_objects_pred/preg_dat_25cov_qval01.rds")























