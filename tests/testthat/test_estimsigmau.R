# Test if the estimation of the variance of the random effects corresponds between
# the sae package and the fh function for reml and ml estimation up to a specific
# level


# Get the benchmark results
benchmark <- read.csv("sae_estimsigma.csv")

# Get data for the comparison
data(eusilcA_smpAgg)
data(eusilcA_popAgg)

# Combine the data sets
combined_data <- combine_data(pop_data = eusilcA_popAgg, pop_domains = "Domain",
                              smp_data = eusilcA_smpAgg, smp_domains = "Domain",
                              vardir = NULL)
combined_data <- combined_data[!is.na(combined_data$Mean),]


test_that("Does the FH_eblup function return the same point and mse estimates like package sae?", {

  # Estimation with FH_eblup function
  fh_test <-  fh(fixed = Mean ~ eqsize + cash, vardir = "Var_Mean",
                 combined_data = combined_data, domains = "idD",
                 method = "reml", interval = c(0, 2800000), transformation = "no")

  # Compare point estimates from FH_AK and benchmark
  expect_equal((fh_test$sigmau2 - benchmark$REML)/benchmark$REML < 1.2e-06,
               TRUE)

  expect_equal((fh_test$sigmau2 - benchmark$REML)/benchmark$REML < 1.2e-06,
               TRUE)

})
