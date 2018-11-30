data(eusilcA_smpAgg)
data(eusilcA_popAgg)


combined_data <- combine_data(pop_data = eusilcA_popAgg, pop_domains = "Domain",
                              smp_data = eusilcA_smpAgg, smp_domains = "Domain",
                              vardir = NULL)
combined_data <- combined_data[!is.na(combined_data$Mean),]

library(sae)
sae_FH <- mseFH(Mean ~ eqsize + cash, vardir = Var_Mean, data = combined_data)
sae_FH$est$fit$refvar

FH_test <- fh(fixed = Mean ~ eqsize + cash, vardir = "Var_Mean",
              combined_data = combined_data, domains = "idD",
              method = "reml", interval = c(0, 2800000), transformation = "no",
              eff_smpsize = "effsample", alpha = 0.05)

FH_test$sigmau2
all.equal(FH_test$sigmau2, sae_FH$est$fit$refvar)

all.equal((FH_test$sigmau2 - sae_FH$est$fit$refvar)/sae_FH$est$fit$refvar < 1.2e-06,
          TRUE)




sae_FH_ml <- mseFH(Mean ~ eqsize + cash, vardir = Var_Mean, data = combined_data,
                method = "ML")

FH_test_ml <- fh(fixed = Mean ~ eqsize + cash, vardir = "Var_Mean",
              combined_data = combined_data, domains = "idD",
              method = "ml", interval = c(0, 2800000), transformation = "no",
              eff_smpsize = "effsample", alpha = 0.05)

all.equal((sae_FH_ml$est$fit$refvar - FH_test_ml$sigmau2)/FH_test_ml$sigmau2 < 2.98e-06,
          TRUE)


sae_estimsigma <- data.frame(REML = sae_FH$est$fit$refvar,
                             ML = sae_FH_ml$est$fit$refvar)

write.csv(sae_estimsigma, file = "./tests/testthat/sae_estimsigma.csv")
