# CREATED: 2025-07-17
# Create rset object for 10-fold cross-validation



# LIBRARIES ----
library(tidymodels)
library(minfi)




# LOAD DATA ----
# Baseline data from final workflow (all data)
# preg_dat <- readr::read_rds("data_objects/preg_dat_25cov_qval01.rds")
# length(unique(preg_dat$study_id))


# Processed methylation and pheno data
gset.quant <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects/gset.quant.578.rds")
pdata <- as.data.frame(pData(gset.quant))




# DATA SPLITTING AND RESAMPLING ----
set.seed(123456)

#resampling
preg_rs <- group_vfold_cv(pdata, 
                          group = study_id,
                          balance = "groups",
                          v = 10)


# Verify the grouping
# fold_01 <- preg_rs$splits[[1]]
# train_x <- analysis(fold_01)          #Baseline CpG values
# test_x <- assessment(fold_01)      #Baseline CpG values
# train_subjects <- unique(train_x$study_id)
# test_subjects <- unique(test_x$study_id)
# intersect(train_subjects, test_subjects)


# Function to create CV training and test sets
create_cv_dat <- function(rs, fold, GRset) {
  fold_x <- preg_rs$splits[[fold]]
  train_x <- analysis(fold_x)
  test_x <- assessment(fold_x)
  
  train_x_keep <- as.numeric(as.character(unique(train_x$study_id)))
  train_gset <- GRset[, pdata$study_id %in% train_x_keep]
  
  test_x_keep <- as.numeric(as.character(unique(test_x$study_id)))
  test_gset <- GRset[, pdata$study_id %in% test_x_keep]
  

  return(list(train_gset, test_gset))
  }


# Loop to generate all cv objects
# nb: not feasible to put in a single list since these object are large
for (i in 1:10) {
  cv_train_name <- paste0("cv_train_", i)
  cv_test_name <- paste0("cv_test_", i)
  
  cv_x <- create_cv_dat(rs    = preg_rs, 
                        fold  = i,
                        GRset = gset.quant)

  assign(cv_train_name, cv_x[[1]])
  assign(cv_test_name, cv_x[[2]])
}


# Write objects to file
for (i in 1:10) {
  readr::write_rds(get(paste0("cv_train_", i)), file = paste0("data_objects_pred/cv_train_", i, ".rds"))
  readr::write_rds(get(paste0("cv_test_", i)), file = paste0("data_objects_pred/cv_test_", i, ".rds"))
}

