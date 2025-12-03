# CREATED: 2025-08-18
# Mixed-model assessment of change in DNAm over pregnancy
# Uses continuous ga for `time`

# CONFIGURED TO RUN ON ATHENA SLURM
# Builds upon `09_DNAm-timepoint_day-ewas.R` on local
# Copied from: `05_mwas-timepoint-day_25cov_578n.R` in Athena preg-maternal

# Includes minimal (N = ) number of covariates plus the data PCs.
# Data PCs estimated from complete residulized data.
# Not includes processing_date and row. Omits slide.



# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# DNAm object
library(minfi)

# Wrangle
library(dplyr)

# Manage env
library(tictoc)

# Model dev
library(lme4)
library(lmerTest)
library(broom.mixed)
library(qvalue)


# SET CV FOLD ----
fold <- 8


# LOAD DATA ----
gset.quant <- readr::read_rds(paste0("/lustre/home/tpyork/projects/preg-maternal/data_objects_pred/cv_train_", fold, ".rds"))



# PROCESS DATA ----
Mval <- getM(gset.quant)


pdata <- as_tibble(pData(gset.quant)) %>%
  select(study_id, ga_ad, timepoint_day,
         processed_date, row, mMed, uMed, 
         PCtech1, PCtech2, PCtech3, PCtech4, PCtech5, 
         PCdata1, PCdata2, PCdata3, PCdata4, PCdata5, PCdata6, PCdata7, PCdata8, PCdata9, PCdata10,
         PCanc3,
         CD8T, NK, Bcell, Mono, Gran,
         preg_smoke, mat_age_years) %>% 
  mutate(study_id = as.factor(study_id),
         ga_ad    = ga_ad/10)

Gout.dmp <- granges(gset.quant)      #[1:2000]

rm(gset.quant)


# ORTHONORMALIZE ----
# Use ramwas to orthonormalize predictors
ptemp.ortho <- ramwas::orthonormalizeCovariates(
  pdata %>% 
    select(-study_id, -ga_ad),
  modelhasconstant = FALSE)

pdata <- cbind(pdata[, c('study_id', 'ga_ad')], ptemp.ortho)
names(pdata) <- c('study_id', 'ga_ad', paste0("c", 1:dim(ptemp.ortho)[2]))
names(pdata)[3] <- "timepoint_day"






# FUNCTION: mixed model analysis
model_nocov <- function(mval, pheno) {
  
  pheno2 <- cbind(mval, pheno)
  
  mod1 <- lmer(mval ~ timepoint_day + 
                 c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 +
                 c11 + c12 + c13 + c14 + c15 + c16 + c17 + c18 + c19 +
                 c20 + c21 + c22 + c23 + c24 + c25 + c26 + c27 + c28 +
                 c29 + c30 + c31 + 
                 (1 + timepoint_day | study_id), data = pheno2)
  mod1.tidy <- tidy(mod1)
  
  return(format(unlist(mod1.tidy[2, c(4, 5, 6, 8)]), scientific = TRUE))
}



# RUN MODEL LOOP ----

# ONE
start <- 1
stop <- 100000
out.1 <- matrix(NA, stop-start+1, 4)
count <- 0
tic("DNAm_1")
for (i in start:stop) {
  count <- count + 1
  if(i %% 1000 == 0) print(i)
  out.1[count,] <- xpectr::suppress_mw(model_nocov(mval = Mval[i,], pheno = pdata))
}
toc()


# TWO
start <- 100001
stop <- 200000
out.2 <- matrix(NA, stop-start+1, 4)
count <- 0
tic("DNAm_2")
for (i in start:stop) {
  count <- count + 1
  if(i %% 1000 == 0) print(i)
  out.2[count,] <- xpectr::suppress_mw(model_nocov(Mval[i,], pdata))
}
toc()


#THREE
start <- 200001
stop <- 300000
out.3 <- matrix(NA, stop-start+1, 4)
count <- 0
tic("DNAm_3")
for (i in start:stop) {
  count <- count + 1
  if(i %% 1000 == 0) print(i)
  out.3[count,] <- xpectr::suppress_mw(model_nocov(Mval[i,], pdata))
}
toc()


#FOUR
start <- 300001    #31
stop <-  435982
out.4 <- matrix(NA, stop-start+1, 4)
count <- 0
tic("DNAm_4")
for (i in start:stop) {
  count <- count + 1
  if(i %% 1000 == 0) print(i)
  out.4[count,] <- xpectr::suppress_mw(model_nocov(Mval[i,], pdata))
}
toc()


out <- rbind(out.1, out.2, out.3, out.4)


# create output object
Gout.dmp$mix.coef <- as.numeric(out[, 1])   #DNAm change
Gout.dmp$mix.se   <- as.numeric(out[, 2])
Gout.dmp$mix.stat <- as.numeric(out[, 3])
Gout.dmp$mix.pval <- as.numeric(out[, 4])




# INCLUDE FDR -------------------------------------------------------------
Gout.dmp$mix.qval <- qvalue(Gout.dmp$mix.pval)$qvalue




# OUTPUT ------------------------------------------------------------------
readr::write_rds(Gout.dmp, file= paste0("/lustre/home/tpyork/projects/preg-maternal/data_objects_pred/mwas_cv_train_", fold, ".rds"))




