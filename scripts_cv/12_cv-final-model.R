# CREATED: 2025-10-15
# Establish Final Model on Full Data


# Get full data MWAS results; bacon corrected 
# Select top 3000 CpGs based on CV results
# Fit elastic net model; re-estimate hyperparameters via 10-fold CV; select best and refit on entire dataset
# Show correlation / scatterplot of observed vs predicted GAD
# Use VIP to identify main drivers of predictor

# Data Summary Flowchart
# FILE  | OBJECT  |  COMMENT

# data_objects/gset.quant.578.rds          | gset.quant | processed methylation data     |
# data_objects/dmp-timepoint_day-25cov.rds | Gout.dmp   | MWAS results on full data      |
# data_objects/Gout-25cov.rds              | cov25      | Gout.dmp with bacon adjustment |
# sdf | sdf | sdf |
# sdf | sdf | sdf |



# LIBRARIES ----
library(tidymodels)
tidymodels_prefer()
library(minfi)
library(yardstick)
library(vip) 




# This is the variables used
# data is orthonormalized

# pdata <- as_tibble(pData(gset.quant)) %>%
#   select(study_id, ga_ad, timepoint_day,
#          processed_date, row, mMed, uMed, 
#          PCtech1, PCtech2, PCtech3, PCtech4, PCtech5, 
#          PCdata1, PCdata2, PCdata3, PCdata4, PCdata5, PCdata6, PCdata7, PCdata8, PCdata9, PCdata10,
#          PCanc3,
#          CD8T, NK, Bcell, Mono, Gran,
#          preg_smoke, mat_age_years) %>% 
#   mutate(study_id = as.factor(study_id),
#          ga_ad    = ga_ad/10)


# PARAMETER SETTINGS ----
cpg_max <- 3000
tune_size <- 500




# LOAD MWAS RESULTS FROM FULL SAMPLE ----
gout <- readr::read_rds("data_objects/dmp-timepoint_day-25cov.rds")


# Filter CpG results by pval
cpg_keep <- gout %>% as.data.frame() %>%
  arrange(mix.pval) %>% slice(1:cpg_max) %>% rownames()
cpg_keep %>% length()


# LOAD SUBJECT-LEVEL PHENO AND CpG VALUES FOR TOP 3000 FEATURES ----
# 2. GATHER CPG TRAIN and TEST TIMEPOINT V1 DATA ----
gset.quant <- readr::read_rds("data_objects/gset.quant.578.rds")


cpg_id <- granges(gset.quant) %>% 
  names()

mval <- getM(gset.quant) %>% 
  as.data.frame() %>% 
  t()

names(mval) <- cpg_id

mval <- mval %>% 
  as_tibble()

mval <- mval[, cpg_keep]


# 3. GATHER PHENOTYPIC DATA ----
# train
pdata <- as_tibble(pData(gset.quant)) %>%
  select(study_id, ga_ad, timepoint_day, visit,
         preg_smoke, mat_age_years, dm_16,
         CD8T, NK, Bcell, Mono, Gran) %>% 
  mutate(study_id   = as.factor(study_id),
         preg_smoke = as.factor(preg_smoke),
         dm_16      = as.factor(dm_16))


preg <- bind_cols(pdata, mval) %>% 
  filter(visit == 1)



# APPLY MODEL TO FULL DATA ----

# 5. MODEL SETUP ----
# Data splitting and resampling
set.seed(123)
preg_rs <- vfold_cv(preg, v = 10)  #preg_tr


# Model specification
# Elastic Net
elastic_mod <- linear_reg(
  penalty = tune(),
  mixture = tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")

param_set <- extract_parameter_set_dials(elastic_mod)
# dials::penalty()
# dials::mixture()

# Recipe
lr_rec <-
  recipe(ga_ad ~ ., data = preg) %>% 
  # step_rm(study_id, timepoint_day, visit) %>%
  step_rm(study_id, timepoint_day, visit, preg_smoke, mat_age_years, dm_16,
          CD8T, NK, Bcell, Mono, Gran) %>%
  step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors())

# Workflow
elastic_wrk <-
  workflow() %>% 
  add_model(elastic_mod) %>% 
  add_recipe(lr_rec)


# 6. TUNE ----
elastic_grid_random <- grid_random(
  penalty(),
  mixture(),
  size = tune_size)

set.seed(234)
elastic_res <- tune_grid(
  elastic_wrk,
  resamples = preg_rs,
  grid      = elastic_grid_random
)
readr::write_rds(elastic_res, file = "data_objects_final/elastic_res_final.rds")


#Choose model parameter set
elastic_best <- elastic_res %>% 
  select_best(metric = "rsq")

# elastic_best <- elastic_res %>% 
#     select_best(metric = "rmse")



# 7. FINAL MODEL ----
# Fit best model on training
elastic_workflow_final <- finalize_workflow(
  elastic_wrk,
  elastic_best
)



#Fit on training with selected parameters on training
fit_final <- 
  elastic_workflow_final %>% 
  fit(data = preg)

readr::write_rds(fit_final, file = "data_objects_final/fit_final.rds")


# variable importance plot
fit_final %>% 
  # pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip(num_features = 50, 
      aesthetics = list(color = 'blue', fill = "grey50", size = 0.3))

fit_final %>% 
  # pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip::vi() %>% 
  readr::write_rds("data_objects_final/vip.rds")












