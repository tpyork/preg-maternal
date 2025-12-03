# CREATED: 2025-08-21
# Collect mwas results






# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# library(minfi)
library(minfi)
library(ramwas)
library(bacon)


# MWAS-TIMEPOINT-DAY-16COV ----
# mwas1 <- readr::read_rds("data_objects_pred/mwas_cv_train_1.rds")
# 
# table(mwas1$mix.qval < 0.01)
# table(mwas1$mix.qval < 0.0001)
# 
# table(is.na(mwas1$mix.qval))



mwas_bacon <- function(mwas, fold) {
  p_values <- mwas$mix.pval
  
  z_scores <- qnorm(p_values / 2, lower.tail = FALSE) * sign(rnorm(length(p_values)))
  bc <- bacon(z_scores)
  corrected_pvals <- pval(bc)
  
  qval <- qvalue::qvalue(corrected_pvals)$qvalue
  
  mwas$mix.pval.bac <- corrected_pvals
  mwas$mix.qval.bac <- qval
  
  # Output bacon corrected values
  pdf(paste0("/lustre/home/tpyork/projects/preg-maternal/figures_pred/qqplots_cv", fold, ".pdf"))
    ramwas::qqPlotFast(mwas$mix.pval, ylim = c(0,6))
    ramwas::qqPlotFast(mwas$mix.pval.bac, ylim = c(0,6))
  dev.off()
  
  # Replace object
  readr::write_rds(mwas, file= paste0("/lustre/home/tpyork/projects/preg-maternal/data_objects_pred/mwas_cv_train_", fold, ".rds"))
}



mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_1.rds"), fold = 1)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_2.rds"), fold = 2)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_3.rds"), fold = 3)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_4.rds"), fold = 4)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_5.rds"), fold = 5)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_6.rds"), fold = 6)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_7.rds"), fold = 7)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_8.rds"), fold = 8)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_9.rds"), fold = 9)
mwas_bacon(mwas = readr::read_rds("data_objects_pred/mwas_cv_train_10.rds"), fold = 10)





