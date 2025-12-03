# CREATED: 2025-03-12
# Collect mwas results






# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# library(minfi)
library(minfi)
library(ramwas)




# MWAS-TIMEPOINT-DAY-0COV ----
cov0 <- readr::read_rds("data_objects/dmp-timepoint_day-0cov.rds")

table(cov0$mix.qval < 0.05)
table(cov0$int.qval < 0.05)
table(cov0$slp.qval < 0.05)

table(is.na(cov0$mix.qval))
table(is.na(cov0$int.qval))
table(is.na(cov0$slp.qval))

qqPlotFast(cov0$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov0$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov0$slp.pval), ylim = c(0,6))



# MWAS-TIMEPOINT-DAY-16COV ----
cov16 <- readr::read_rds("data_objects/dmp-timepoint_day-16cov.rds")

table(cov16$mix.qval < 0.01)
table(cov16$int.qval < 0.01)
table(cov16$slp.qval < 0.01)

table(cov16$mix.qval < 0.0001)
table(cov16$int.qval < 0.01)
table(cov16$slp.qval < 0.01)

table(is.na(cov16$mix.qval))
table(is.na(cov16$int.qval))
table(is.na(cov16$slp.qval))

qqPlotFast(cov16$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov16$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov16$slp.pval), ylim = c(0,6))




# MWAS-TIMEPOINT-DAY-19COV ----
cov19 <- readr::read_rds("data_objects/dmp-timepoint_day-19cov.rds")

table(cov19$mix.qval < 0.01)
table(cov19$int.qval < 0.01)
table(cov19$slp.qval < 0.01)

table(cov19$mix.qval < 0.0001)
table(cov19$int.qval < 0.01)
table(cov19$slp.qval < 0.01)

table(is.na(cov19$mix.qval))
table(is.na(cov19$int.qval))
table(is.na(cov19$slp.qval))

qqPlotFast(cov19$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov19$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov19$slp.pval), ylim = c(0,6))




# MWAS-TIMEPOINT-DAY-20COV ----
cov20 <- readr::read_rds("data_objects/dmp-timepoint_day-20cov.rds")

table(cov20$mix.qval < 0.01)
table(cov20$int.qval < 0.01)
table(cov20$slp.qval < 0.01)

table(cov20$mix.qval < 0.0001)
table(cov20$int.qval < 0.01)
table(cov20$slp.qval < 0.01)

table(is.na(cov20$mix.qval))
table(is.na(cov20$int.qval))
table(is.na(cov20$slp.qval))

qqPlotFast(cov20$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov20$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov20$slp.pval), ylim = c(0,6))

table(cov20$mix.qval[1:400] < .01)


# MWAS-TIMEPOINT-DAY-21COV ----
cov21 <- readr::read_rds("data_objects/dmp-timepoint_day-21cov.rds")

table(cov21$mix.qval < 0.01)
table(cov21$int.qval < 0.01)
table(cov21$slp.qval < 0.01)

table(cov21$mix.qval < 0.0001)
table(cov21$int.qval < 0.01)
table(cov21$slp.qval < 0.01)

table(is.na(cov21$mix.qval))
table(is.na(cov21$int.qval))
table(is.na(cov21$slp.qval))

qqPlotFast(cov21$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov21$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov21$slp.pval), ylim = c(0,6))



# MWAS-TIMEPOINT-DAY-22COV ----
cov22 <- readr::read_rds("data_objects/dmp-timepoint_day-22cov.rds")

table(cov22$mix.qval < 0.01)
table(cov22$int.qval < 0.01)
table(cov22$slp.qval < 0.01)

table(cov22$mix.qval < 0.0001)
table(cov22$int.qval < 0.01)
table(cov22$slp.qval < 0.01)

table(is.na(cov22$mix.qval))
table(is.na(cov22$int.qval))
table(is.na(cov22$slp.qval))

qqPlotFast(cov22$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov22$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov22$slp.pval), ylim = c(0,6))



# MWAS-TIMEPOINT-DAY-23COV ----
cov23 <- readr::read_rds("data_objects/dmp-timepoint_day-23cov.rds")

table(cov23$mix.qval < 0.01)
table(cov23$int.qval < 0.01)
table(cov23$slp.qval < 0.01)

table(cov23$mix.qval < 0.0001)
table(cov23$int.qval < 0.01)
table(cov23$slp.qval < 0.01)

table(is.na(cov23$mix.qval))
table(is.na(cov23$int.qval))
table(is.na(cov23$slp.qval))

qqPlotFast(cov23$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov23$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov23$slp.pval), ylim = c(0,6))


# MWAS-TIMEPOINT-DAY-24COV ----
cov24 <- readr::read_rds("data_objects/dmp-timepoint_day-24cov.rds")

table(cov24$mix.qval < 0.01)
table(cov24$int.qval < 0.01)
table(cov24$slp.qval < 0.01)

table(cov24$mix.qval < 0.0001)
table(cov24$int.qval < 0.01)
table(cov24$slp.qval < 0.01)

table(is.na(cov24$mix.qval))
table(is.na(cov24$int.qval))
table(is.na(cov24$slp.qval))

qqPlotFast(cov24$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov24$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov24$slp.pval), ylim = c(0,6))


# MWAS-TIMEPOINT-DAY-25COV ----
cov25 <- readr::read_rds("data_objects/dmp-timepoint_day-25cov.rds")

min(cov25$mix.pval)
max(cov25$mix.pval)
hist(cov25$mix.pval)


table(cov25$mix.qval < 0.01)
table(cov25$int.qval < 0.01)
table(cov25$slp.qval < 0.01)

table(cov25$mix.qval < 0.05)
table(cov25$int.qval < 0.01)
table(cov25$slp.qval < 0.01)

table(is.na(cov25$mix.qval))
table(is.na(cov25$int.qval))
table(is.na(cov25$slp.qval))

qqPlotFast(cov25$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov25$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov25$slp.pval), ylim = c(0,6))


# MWAS-TIMEPOINT-DAY-25COV ----
cov26 <- readr::read_rds("data_objects/dmp-timepoint_day-26cov.rds")

min(cov26$mix.pval)
max(cov26$mix.pval)
hist(cov26$mix.pval)


table(cov26$mix.qval < 0.01)
table(cov26$int.qval < 0.01)
table(cov26$slp.qval < 0.01)

table(cov26$mix.qval < 0.05)
table(cov26$int.qval < 0.01)
table(cov26$slp.qval < 0.01)

table(is.na(cov26$mix.qval))
table(is.na(cov26$int.qval))
table(is.na(cov26$slp.qval))

qqPlotFast(cov26$mix.pval, ylim = c(0,6))
qqPlotFast(na.omit(cov26$int.pval), ylim = c(0,6))
qqPlotFast(na.omit(cov26$slp.pval), ylim = c(0,6))






# bacon ----
library(bacon)
p_values <- cov25$mix.pval
table(cov26$mix.pval == 0)

### This is incorrect - need actual 'sign' information for each coefficient
z_scores <- qnorm(p_values / 2, lower.tail = FALSE) * sign(rnorm(length(p_values)))
# z_scores <- qnorm(p_values / 2, lower.tail = FALSE)
# signed_z <- sign(effect_sizes) * z_score



bc <- bacon(z_scores)
corrected_pvals <- pval(bc)

ramwas::qqPlotFast(corrected_pvals, ylim = c(0,6))

qval <- qvalue::qvalue(corrected_pvals)$qvalue
table(qval < 0.01)

cov25$mix.pval.bac <- corrected_pvals
cov25$mix.qval.bac <- qval
table(cov25$mix.qval.bac < 0.01)
 
# Output bacon corrected values
# readr::write_rds(cov25, file= "/lustre/home/tpyork/projects/preg-maternal/data_objects/Gout-25cov.rds")
# cov25 <- readr::read_rds(file= "/lustre/home/tpyork/projects/preg-maternal/data_objects/Gout-25cov.rds")
# ramwas::qqPlotFast(cov25$mix.pval, ylim = c(0,6))
# ramwas::qqPlotFast(cov25$mix.pval.bac, ylim = c(0,6))
# table(cov25$mix.qval.bac < 0.01)

# Better coding needed here:
readr::write_rds(cov25, file= "/lustre/home/tpyork/projects/preg-maternal/data_objects/Gout-25cov.rds")
readr::write_rds(cov26, file= "/lustre/home/tpyork/projects/preg-maternal/data_objects/Gout-26cov.rds")






# test_statistic <- 0.516
# df <- 475.23296
# p_value <- 2 * pt(-abs(test_statistic), df)
# 
# 
# test_statistic <- 10
# standard_error <- 0.5
# z_score <- test_statistic/standard_error
# p_value <- 2 * (1 - pnorm(abs(z_score)))
# 
# significant_digits <- 4
# signif(p_value, digits = 12)
# 
# format(p_value, scientific = TRUE)
# 
# 
# 
# 
# which(cov25$mix.pval == 0)
# 

