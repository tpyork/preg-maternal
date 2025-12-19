# CREATED: 2025-12-19
# PREG Maternal DNAm Demographic Table



# LIBRARIES ----
library(tidyverse)
library(minfi)

library(gt)
library(gtExtras)
library(gtsummary)


# SET TABLE THEME ----
theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()





# PREPARE DATA ----
gset.quant <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects/gset.quant.578.rds")
pdata <- as_tibble(pData(gset.quant))

rm(gset.quant)

# few PTBs
table(
  pdata[pdata$visit==1, "ga_ad"] < 7*37
)

pdata_tbl <- pdata %>% 
  mutate(preg_smoke = factor(preg_smoke, labels = c("No", "Yes"))) %>% 
  mutate(dm_16 = factor(dm_16, labels = c("White", "Black"))) %>% 
  mutate(visit = factor(visit, labels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))) %>% 
  mutate(timepoint_day = as.numeric(timepoint_day)) %>% 
  select(visit, timepoint_day, ga_ad, mat_age_years, preg_smoke, dm_16)


table_1 <- 
  pdata_tbl %>%
  tbl_summary(
    by = visit,                                    #split table by group
    missing_text = "Prefer not to answer/Unknown",
    label       = list(ga_ad                ~ 'Gestational age at delivery (days)',
                       mat_age_years        ~ 'Maternal age',
                       preg_smoke           ~ 'Smoking during pregnancy',
                       dm_16                ~ 'Self-reported Race',
                       timepoint_day        ~ 'Gestational age at assessment (days)'

    )
  ) %>%
  # add_n() %>%                                       #add column with total number of non-missing observations
  # add_overall() %>% 
  # add_p() %>%                                         #test for a difference between groups
  # separate_p_footnotes() %>% 
  modify_header(label  = "**Variable**") %>%           #update the column header
  bold_labels() 

table_1


# EXPORT ----

# Word format
library(flextable)
library(officer)

table_1 %>% 
  as_flex_table() %>%
  flextable::save_as_docx(path = "figures/table_1.docx")

















