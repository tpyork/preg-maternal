



library(ramwas)


load("initial-analysis-temp/dmp-timepoint_day-mixed_standardized.Rdata")
pval <- unname(Gout.dmp$mix.pval)
qqPlotFast(pval, ylim = c(0, 6))
title('QQ-plot\nwith covariates and 2 PC')

table(Gout.dmp$mix.qval < 0.01)



load("initial-analysis-temp/dmp-timepoint_day-mixed.Rdata")
pval <- unname(Gout.dmp$mix.pval)
qqPlotFast(pval, ylim = c(0, 6))
title('QQ-plot\nwith covariates and 2 PC')

table(Gout.dmp$mix.qval < 0.01)






