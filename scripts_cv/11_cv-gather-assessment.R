# CREATED: 2025-09-12
# Gather and Evaluate Prediction Model Assessments




# SET UP WORKSPACE AND LOAD PACKAGES ----
library(tidyverse)
library(readr)



cpg5 <- read_rds(file = "data_objects_cv/assess-cpg_max5-tune_size500-model_bestRSQ.rds")
cpg10 <- read_rds(file = "data_objects_cv/assess-cpg_max10-tune_size500-model_bestRSQ.rds")
cpg20 <- read_rds(file = "data_objects_cv/assess-cpg_max20-tune_size500-model_bestRSQ.rds")
cpg50 <- read_rds(file = "data_objects_cv/assess-cpg_max50-tune_size500-model_bestRSQ.rds")
cpg100 <- read_rds(file = "data_objects_cv/assess-cpg_max100-tune_size500-model_bestRSQ.rds")
cpg200 <- read_rds(file = "data_objects_cv/assess-cpg_max200-tune_size500-model_bestRSQ.rds")
cpg500 <- read_rds(file = "data_objects_cv/assess-cpg_max500-tune_size500-model_bestRSQ.rds")
cpg1000 <- read_rds(file = "data_objects_cv/assess-cpg_max1000-tune_size500-model_bestRSQ.rds")
cpg1500 <- read_rds(file = "data_objects_cv/assess-cpg_max1500-tune_size500-model_bestRSQ.rds")
cpg2000 <- read_rds(file = "data_objects_cv/assess-cpg_max2000-tune_size500-model_bestRSQ.rds")
cpg3000 <- read_rds(file = "data_objects_cv/assess-cpg_max3000-tune_size500-model_bestRSQ.rds")
cpg4000 <- read_rds(file = "data_objects_cv/assess-cpg_max4000-tune_size500-model_bestRSQ.rds")
cpg5000 <- read_rds(file = "data_objects_cv/assess-cpg_max5000-tune_size500-model_bestRSQ.rds")





cpg5_vector <- map_dbl(cpg5, ~ .x[[2, 3]])
cpg10_vector <- map_dbl(cpg10, ~ .x[[2, 3]])
cpg20_vector <- map_dbl(cpg20, ~ .x[[2, 3]])
cpg50_vector <- map_dbl(cpg50, ~ .x[[2, 3]])
cpg100_vector <- map_dbl(cpg100, ~ .x[[2, 3]])
cpg200_vector <- map_dbl(cpg200, ~ .x[[2, 3]])
cpg500_vector <- map_dbl(cpg500, ~ .x[[2, 3]])
cpg1000_vector <- map_dbl(cpg1000, ~ .x[[2, 3]])
cpg1500_vector <- map_dbl(cpg1500, ~ .x[[2, 3]])
cpg2000_vector <- map_dbl(cpg2000, ~ .x[[2, 3]])
cpg3000_vector <- map_dbl(cpg3000, ~ .x[[2, 3]])
cpg4000_vector <- map_dbl(cpg4000, ~ .x[[2, 3]])
cpg5000_vector <- map_dbl(cpg5000, ~ .x[[2, 3]])



out_tbl <- tibble(
    rsq     = c(cpg5_vector, cpg10_vector, cpg20_vector, cpg50_vector, cpg100_vector, 
                cpg200_vector, cpg500_vector ,cpg1000_vector, cpg1500_vector, cpg2000_vector, cpg3000_vector, cpg4000_vector, cpg5000_vector),
    markers = rep(c("5", "10", "20", "50", "100", "200", "500", "1000", "1500", "2000", "3000", "4000", "5000"), each = 10)
  )


a <- 
out_tbl %>% 
  group_by(markers) %>% 
  summarise(mean_rsq = mean(rsq, na.rm = TRUE)) %>% 
  mutate(markers = factor(markers, level = c("5", "10", "20", "50", "100", "200", "500", "1000", "1500", "2000", "3000", "4000", "5000"))) %>% 
  ggplot(aes(x = markers, y = sqrt(mean_rsq))) +
  geom_point(size = 3) +
  labs(title = "Prediction of GAD", x = "Number of Markers", y = "Correlation") +
  theme_minimal()

ggsave(path     ="figures_cv",
       filename = "CV_number-of-markers.jpeg",
       device   = "jpeg")

pdf("figures_cv/CV_number-of-markers.pdf")
  print(a)
dev.off()

















