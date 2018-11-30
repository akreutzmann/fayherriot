#===================================================================================================
# SUBJECT:  Example code for Fay-Herriot Function with arcsin transformation
# AUTHOR:   Timo Schmid, Fabian Bruckschen, Nicola Salvati and Till Zbiranski, June 06, 2017
# TITLE:    Example code for the paper "Constructing socio-demographic indicators for National
#           Statistical Institutes using mobile phone data: estimating literacy rates in Senegal"
#===================================================================================================

# Load Packages ------------------------------------------------------------------------------------
require(MASS)
require(nlme)
require(formula.tools)

# Load Data ----------------------------------------------------------------------------------------
setwd("C:/ here your / directory path / for the working example/")
load("testdata.RData")

# Source Function ----------------------------------------------------------------------------------
source("FH-function-arcsin.R")

# Run Function -------------------------------------------------------------------------------------


fit<-FH_arcsin(formula = y~x1+x2,dataframe_sample = data_sample, saind = data_sample$id,
                              dataframe_pop_aux = data_pop,x.total_saind = data_pop$id,
                              method=c("AP"), area_count = data_sample$effsample, Boot =T,
                              B = 20,ratio = weights, bench_value = national_mean)

# Show Estimates -----------------------------------------------------------------------------------
# Small area point estimates
fit$est_mean_FH
fit$sigma2v
fit$CI_FH$Low
fit$CI_FH$Up

combined_data <- combine_data(pop_data = data_pop, pop_domains = "id",
                              smp_data = data_sample, smp_domains = "id",
                              vardir = NULL)

FH_test <- fh(fixed = y ~ x1.y + x2.y, vardir = "vardir",
         combined_data = combined_data, domains = "idD",
         method = "ampl", interval = c(0, 100), transformation = "arcsin",
         eff_smpsize = "effsample", alpha = 0.05)

FH_test$sigmau2
FH_test$MSE$Li
FH_test$MSE$Ui

all.equal(fit$est_mean_FH$pred, FH_test$ind$EBLUP)
all.equal(fit$CI_FH$Low, FH_test$MSE$Li)
all.equal(fit$CI_FH$Up, FH_test$MSE$Ui)

# Corresponding confidence intervals
fit$CI_FH

# Benchmarked small area point estimates
fit$est_mean_FH_bench
