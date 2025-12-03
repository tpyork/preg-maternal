# CREATED: 2025-03-21
# PREG MATERNAL
# Remove sample outliers loading high on data PCs



# LIBRARIES ----
library(minfi)
library(tidyverse)


# LOAD DATA ----
gset.temp <- readr::read_rds("/lustre/home/tpyork/projects/preg-maternal/data_objects/gset.quant.5.rds")


ptemp <- as_tibble(pData(gset.temp))

out1 <- which(abs(ptemp %>% select(PC1outlier)) > 3)
out2 <- which(abs(ptemp %>% select(PC2outlier)) > 3)
out3 <- which(abs(ptemp %>% select(PC3outlier)) > 3)
out4 <- which(abs(ptemp %>% select(PC4outlier)) > 3)
out5<- which(abs(ptemp %>% select(PC5outlier)) > 3)
out6 <- which(abs(ptemp %>% select(PC6outlier)) > 3)
out7 <- which(abs(ptemp %>% select(PC7outlier)) > 3)
out8 <- which(abs(ptemp %>% select(PC8outlier)) > 3)
out9 <- which(abs(ptemp %>% select(PC9outlier)) > 3)
out10 <- which(abs(ptemp %>% select(PC10outlier)) > 3)

outlier <- c(out1, out2, out3, out4, out5, out6, out7, out8, out9, out10)
table(outlier)
length(unique(outlier))
outlier47 <-  unique(outlier)

gset.quant.578 <- gset.temp[, -outlier47]

write_rds(gset.quant.578, file = "data_objects/gset.quant.578.rds")




