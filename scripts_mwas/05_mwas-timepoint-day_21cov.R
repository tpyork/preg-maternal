# CREATED: 2025-03-12
# Mixed-model assessment of change in DNAm over pregnancy
# Uses continuous ga for `time`

# CONFIGURED TO RUN ON ATHENA SLURM
# Builds upon `09_DNAm-timepoint_day-ewas.R` on local

# Includes minimal (N = 21) number of covariates plus the data PCs.
# Includes ga_ad as covariate.
# Data PCs estimated from complete residulized data.
# Not includes processing_date and row. Omits slide.



# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# DNAm object
library(minfi)

# Wrangle
library(dplyr)

# Manage env
# library(doParallel)
# library(foreach)
library(tictoc)

# Model dev
library(lme4)
library(lmerTest)
library(broom.mixed)
library(qvalue)



# LOAD DATA ----
gset.quant <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects/gset.quant.2.rds")



# pheno$timepoint_day + pheno$Comp.Anc.3 + pheno$mat_age_years + pheno$parous +
#   pheno$preg_smoke + pheno$Gran + pheno$Bcell + pheno$Mono +
#   (1 + pheno$timepoint_day | pheno$study_id))


# PROCESS DATA ----
Mval <- getM(gset.quant)

pdata <- as_tibble(pData(gset.quant)) %>%
  select(study_id, timepoint_day,
         ga_ad, 
         mMed, uMed, 
         PC1, PC2, PC3, PC4, PC5,
         PCdata1, PCdata2, PCdata3,
         CD8T, NK, Bcell, Mono, Gran,
         processed_date, row,
         mat_age_years, preg_smoke, dm_16) %>% 
  mutate(study_id = as.factor(study_id),
         ga_ad    = ga_ad/10)

Gout.dmp <- granges(gset.quant)     #[1:400]

rm(gset.quant)


# ORTHONORMALIZE ----
# Use ramwas to orthonormalize predictors
ptemp.ortho <- ramwas::orthonormalizeCovariates(
  pdata %>% 
    select(-study_id),
  modelhasconstant = FALSE)

pdata <- cbind(pdata[, c('study_id')], ptemp.ortho)
names(pdata) <- c('study_id', paste0("c", 1:dim(ptemp.ortho)[2]))
names(pdata)[2] <- "timepoint_day"
names(pdata)[3] <- "ga_ad"


# temp <- 
# as_tibble(pData(gset.quant)) %>%
#   select(study_id, ga_ad, timepoint_day,
#          mMed, uMed, PC1, PC2, PC3, PC4, PC5,
#          CD4T, NK, Bcell, Mono, Gran,
#          Slide,
#          mat_age_years, preg_smoke, dm_16)
# cor(temp$timepoint_day, pdata$c2)
# cor(pdata$c2, pdata$c3)



# 1. MODEL DNAm CHANGE - 19 COVARIATES ----


# mval <- Mval[1,]
# pheno2 <- cbind(mval, pdata)
# 
# mod1 <- lmer(mval ~ timepoint_day +
#                c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 +
#                c21 + c22 + c23 + c24 +
#                (1 + timepoint_day | study_id), data = pheno2)
# summary(mod1)

#### *PICK UP HERE ----
# try without slide but include data PCs



# FUNCTION: 20 covariates
model_nocov <- function(mval, pheno) {
  
  pheno2 <- cbind(mval, pheno)
  
  mod1 <- lmer(mval ~ timepoint_day + 
                 ga_ad + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 +
                 c11 + c12 + c13 + c14 + c15 + c16 + c17 + c18 + c19 +
                 c20 + c21 + c22 + c23 + c24 + c25 +
                 (1 + timepoint_day | study_id), data = pheno2)
  mod1.tidy <- tidy(mod1)
  
  # b. capture effect of change on ga_ad ----
  mod1.rand <- coef(mod1)$study_id[, 1:2]  %>%
    mutate(study_id = as.factor(rownames((coef(mod1)$study_id)))) %>%
    as_tibble() %>%
    rename(rand.slp = `timepoint_day`) %>%
    rename(rand.int = `(Intercept)`)
  
  pheno.sub <- pheno2 %>%
    select(study_id, ga_ad) %>%
    filter(!duplicated(study_id))
  
  dat.rand <-
    left_join(mod1.rand, pheno.sub, by = 'study_id')
  
  mod2.int <- lm(ga_ad ~ rand.int, data = dat.rand)
  mod2.int.tidy <- tidy(mod2.int)
  
  mod2.slp <- lm(ga_ad ~ rand.slp, data = dat.rand)
  mod2.slp.tidy <- tidy(mod2.slp)
  
  return(round(unlist(c(mod1.tidy[2, c(4, 6, 8)], mod2.int.tidy[2, c(2, 4, 5)], mod2.slp.tidy[2, c(2, 4, 5)] )), 5))
}



# RUN MODEL LOOP ----

# ONE
start <- 1
stop <- 100000
out.1 <- matrix(NA, stop-start+1, 9)
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
out.2 <- matrix(NA, stop-start+1, 9)
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
out.3 <- matrix(NA, stop-start+1, 9)
count <- 0
tic("DNAm_3")
for (i in start:stop) {
  count <- count + 1
  if(i %% 1000 == 0) print(i)
  out.3[count,] <- xpectr::suppress_mw(model_nocov(Mval[i,], pdata))
}
toc()


#FOUR
start <- 300001
stop <-  413238
out.4 <- matrix(NA, stop-start+1, 9)
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
Gout.dmp$mix.coef <- round(out[, 1], 5)   #DNAm change
Gout.dmp$mix.stat <- round(out[, 2], 5)
Gout.dmp$mix.pval <- round(out[, 3], 5)
Gout.dmp$int.coef <- round(out[, 4], 5)   #ga_ad on int
Gout.dmp$int.stat <- round(out[, 5], 5)
Gout.dmp$int.pval <- round(out[, 6], 5)
Gout.dmp$slp.coef <- round(out[, 7], 5)   #ga_ad on slp
Gout.dmp$slp.stat <- round(out[, 8], 5)
Gout.dmp$slp.pval <- round(out[, 9], 5)



# INCLUDE FDR -------------------------------------------------------------
Gout.dmp$mix.qval <- qvalue(Gout.dmp$mix.pval)$qvalue
Gout.dmp$int.qval <- qvalue(Gout.dmp$int.pval)$qvalue
Gout.dmp$slp.qval <- qvalue(Gout.dmp$slp.pval)$qvalue


table(Gout.dmp$mix.qval < .01)
table(Gout.dmp$int.qval < .01)
table(Gout.dmp$slp.qval < .01)




# OUTPUT ------------------------------------------------------------------
readr::write_rds(Gout.dmp, file= "/lustre/home/tpyork/projects/preg-maternal/data_objects/dmp-timepoint_day-21cov.rds")




