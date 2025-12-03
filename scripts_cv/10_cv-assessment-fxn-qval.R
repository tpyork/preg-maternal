# CREATED: 2025-09-04
# Prediction Model Assessment Function

# Combines: `08_pred_data_gather.R` and `09_pred_asssessment.R`

### VARIATION THAT BUILDS MODELS BASED ON QVAL INSTEAD OF TOP CpGs




# SET UP WORKSPACE AND LOAD PACKAGES ----
library(minfi)
library(tidyverse)
library(tidymodels)
tidymodels_prefer()
library(yardstick)

setwd("/lustre/home/tpyork/projects/preg-maternal")

pred_assess <- function(cv, qval_thresh, tune_size, model_best = "rsq") {
  
  # 1.GATHER CPG CHANGE RESULTS ----
  gout <- readr::read_rds(paste0("data_objects_cv/mwas_cv_train_", cv, ".rds"))
  
  # Filter CpG results by qval
  cpg_keep <- gout %>% as.data.frame() %>%
    filter(mix.qval.bac < qval_thresh) %>% rownames()
  # cpg_keep %>% length()

  
  # 2. GATHER CPG TRAIN and TEST TIMEPOINT V1 DATA ----
  #train mval
  gset.train <- readr::read_rds(paste0("data_objects_cv/cv_train_", cv, ".rds"))
  
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
  gset.test <- readr::read_rds(paste0("data_objects_cv/cv_test_", cv, ".rds"))
  
  cpg_id <- granges(gset.test) %>% 
    names()
  
  mval.test <- getM(gset.test) %>% 
    as.data.frame() %>% 
    t()
  
  names(mval.test) <- cpg_id
  
  mval.test <- mval.test %>% 
    as_tibble()
  
  mval.test <- mval.test[, cpg_keep]
  
  
  # 3. GATHER PHENOTYPIC DATA ----
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
  
  
  # 4. CREATE DATA FOR MODEL ----
  preg_train <- bind_cols(pdata.train, mval.train) %>% 
    filter(visit == 1)
  # length(unique(preg_dat_train$study_id))
  
  
  preg_test <- bind_cols(pdata.test, mval.test) %>% 
    filter(visit == 1)
  # length(unique(preg_dat_test$study_id))
  
  
  # 5. MODEL SETUP ----
  # Data splitting and resampling
    set.seed(123)
    preg_rs <- vfold_cv(preg_train, v = 10)  #preg_tr
  
  
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
    recipe(ga_ad ~ ., data = preg_train) %>% 
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
  
  #Choose model parameter set
  
  if (model_best == "rsq") {
  elastic_best <- elastic_res %>% 
    select_best(metric = "rsq")
  }
  
  if (model_best == "rmse") {
  elastic_best <- elastic_res %>% 
    select_best(metric = "rmse")
  }
  
  
  # 7. FINAL MODEL ----
  # Fit best model on training
  elastic_workflow_final <- finalize_workflow(
    elastic_wrk,
    elastic_best
  )
  
  #Fit on training with selected parameters on training
  fit_train <- 
    elastic_workflow_final %>% 
    fit(data = preg_train)
  
  
  # 8. ASSESS ON TEST ----
  test_pred <- predict(fit_train, new_data = preg_test) %>% 
    bind_cols(preg_test)
  
  metrics_out <- metrics(test_pred, truth = ga_ad, estimate = .pred)
  
  return(metrics_out)
  
  }


# RUN FUNCTION AND SAVE RESULTS ----


#* Q-VALUE 0.01 ----
# out <- list()
# 
# for (i in 1:10) {
#   out[[i]] <- pred_assess(
#                cv          = i,
#                qval_thresh = 0.01,
#                tune_size   = 500,
#                model_best  = "rsq")
# }
# 
# # Output
# readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_01-tune_size500-model_bestRSQ.rds")


#* Q-VALUE 0.05 ----
# out <- list()
# 
# for (i in 1:10) {
#   out[[i]] <- pred_assess(
#     cv          = i,
#     qval_thresh = 0.05,
#     tune_size   = 500,
#     model_best  = "rsq")
# }
# 
# # Output
# readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_05-tune_size500-model_bestRSQ.rds")


#* Q-VALUE 0.10 ----
# out <- list()
# 
# for (i in 1:10) {
#   out[[i]] <- pred_assess(
#     cv          = i,
#     qval_thresh = 0.10,
#     tune_size   = 500,
#     model_best  = "rsq")
# }
# 
# # Output
# readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_10-tune_size500-model_bestRSQ.rds")



#* Q-VALUE 0.20 ----
# out <- list()
# 
# for (i in 1:10) {
#   out[[i]] <- pred_assess(
#     cv          = i,
#     qval_thresh = 0.20,
#     tune_size   = 500,
#     model_best  = "rsq")
# }
# 
# # Output
# readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_20-tune_size500-model_bestRSQ.rds")



#* Q-VALUE 0.25 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.25,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_25-tune_size500-model_bestRSQ.rds")


#* Q-VALUE 0.30 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.30,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_30-tune_size500-model_bestRSQ.rds")



#* Q-VALUE 0.35 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.35,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_35-tune_size500-model_bestRSQ.rds")




#* Q-VALUE 0.40 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.40,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_40-tune_size500-model_bestRSQ.rds")


#* Q-VALUE 0.45 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.45,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_45-tune_size500-model_bestRSQ.rds")


#* Q-VALUE 0.50 ----
out <- list()

for (i in 1:10) {
  out[[i]] <- pred_assess(
    cv          = i,
    qval_thresh = 0.50,
    tune_size   = 500,
    model_best  = "rsq")
}

# Output
readr::write_rds(out, file = "/lustre/home/tpyork/projects/preg-maternal/data_objects_cv/assess-qval_50-tune_size500-model_bestRSQ.rds")



