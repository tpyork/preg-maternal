# CREATED: 2025-05-01
# Prediction Data Gather

# updated for 10-fold cv

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

#* 1. GATHER CPG CHANGE RESULTS ----
gout <- readr::read_rds("data_objects_cv/mwas_cv_train_7.rds")


# Filter CpG results by pval

# cpg_keep <- gout %>% as.data.frame() %>% 
#   arrange(mix.pval) %>% slice(1:2000) %>% rownames()
# cpg_keep %>% length()


# Filter CpG results by qval
cpg_keep <- gout[gout$mix.qval.bac < 0.05, ] %>%
  as.data.frame() %>% rownames()
cpg_keep %>% length()


#* 2. GATHER CPG TRAIN and TEST TIMEPOINT V1 DATA ----

#train mval
gset.train <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/cv_train_7.rds")

cpg_id <- granges(gset.train) %>% 
  names()

mval.train <- getM(gset.train) %>% 
  as.data.frame() %>% 
  t()

names(mval.train) <- cpg_id

mval.train <- mval.train %>% 
  as_tibble()

mval.train <- mval.train[, cpg_keep]


#test mval
gset.test <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/cv_test_7.rds")

cpg_id <- granges(gset.test) %>% 
  names()

mval.test <- getM(gset.test) %>% 
  as.data.frame() %>% 
  t()

names(mval.test) <- cpg_id

mval.test <- mval.test %>% 
  as_tibble()

mval.test <- mval.test[, cpg_keep]



#* 3. GATHER PHENOTYPIC DATA ----
# train
pdata.train <- as_tibble(pData(gset.train)) %>%
  select(study_id, ga_ad, timepoint_day, visit,
         preg_smoke, mat_age_years, dm_16,
         CD8T, NK, Bcell, Mono, Gran) %>% 
  mutate(study_id   = as.factor(study_id),
         preg_smoke = as.factor(preg_smoke),
         dm_16      = as.factor(dm_16))

# test
pdata.test <- as_tibble(pData(gset.test)) %>%
  select(study_id, ga_ad, timepoint_day, visit,
         preg_smoke, mat_age_years, dm_16,
         CD8T, NK, Bcell, Mono, Gran) %>% 
  mutate(study_id   = as.factor(study_id),
         preg_smoke = as.factor(preg_smoke),
         dm_16      = as.factor(dm_16))


#* 4. CREATE DATA FOR MODEL ----
preg_dat_train <- bind_cols(pdata.train, mval.train) %>% 
  filter(visit == 1)
# length(unique(preg_dat_train$study_id))


preg_dat_test <- bind_cols(pdata.test, mval.test) %>% 
  filter(visit == 1)
# length(unique(preg_dat_test$study_id))



# OUTPUT DATA ----
readr::write_rds(preg_dat_train, file = "data_objects_cv/preg_dat_train_cv7_qval05.rds")
readr::write_rds(preg_dat_test, file = "data_objects_cv/preg_dat_test_cv7_qval05.rds")



















