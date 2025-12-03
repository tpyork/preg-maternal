# CREATED: 2025-09-12
# Gather and Evaluate Prediction Model Assessments




# SET UP WORKSPACE AND LOAD PACKAGES ----
library(tidyverse)
library(readr)



qval01 <- read_rds(file = "data_objects_cv/assess-qval_01-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval01, ~ .x[[2, 3]])
mean(qval_vector)  # 0.02032292


qval05 <- read_rds(file = "data_objects_cv/assess-qval_05-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval05, ~ .x[[2, 3]])
mean(qval_vector)  # 0.053848


qval10 <- read_rds(file = "data_objects_cv/assess-qval_10-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval10, ~ .x[[2, 3]])
mean(qval_vector)  # 0.05968664


qval20 <- read_rds(file = "data_objects_cv/assess-qval_20-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval20, ~ .x[[2, 3]])
mean(qval_vector)  # 0.1075071


qval25 <- read_rds(file = "data_objects_cv/assess-qval_25-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval25, ~ .x[[2, 3]])
mean(qval_vector)  # 0.07869189


qval30 <- read_rds(file = "data_objects_cv/assess-qval_30-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval30, ~ .x[[2, 3]])
mean(qval_vector)  # 0.1034531


qval35 <- read_rds(file = "data_objects_cv/assess-qval_35-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval35, ~ .x[[2, 3]])
mean(qval_vector)  # 0.04303049


qval40 <- read_rds(file = "data_objects_cv/assess-qval_40-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval40, ~ .x[[2, 3]])
mean(qval_vector)  # 0.02888662


qval45 <- read_rds(file = "data_objects_cv/assess-qval_45-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval45, ~ .x[[2, 3]])
mean(qval_vector)  # 0.01944246


qval50 <- read_rds(file = "data_objects_cv/assess-qval_50-tune_size500-model_bestRSQ.rds")
qval_vector <- map_dbl(qval50, ~ .x[[2, 3]])
mean(qval_vector)  # 0.01664624














