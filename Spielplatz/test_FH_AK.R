################################################################################
# Get a new FH function
################################################################################

# Load packages
library(sae)
library(MASS)
library(nlme)
library(formula.tools)


# Load data
combined_data <- read.csv("C:/Users/akreutzmann/Dropbox/QUESSAMI/Kooperation mit DIW/Stata Code für FH/SmallAreaEstimation/IN/c_data.csv")



# Have example with sae code
sae_FH <- mseFH(yi ~ MajorArea + newVar, SD2, data = combined_data)

My_FH <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
      combined_data = combined_data, domains = "SmallArea",
      method = "nicola_reml", interval = c(0, 1000), transformation = "no")

all.equal(as.numeric(sae_FH$est$eblup), as.numeric(My_FH$ind$EBLUP))
all.equal(as.numeric(sae_FH$mse), as.numeric(My_FH$MSE$MSE))

c(sae_FH$est$eblup)
as.numeric(My_FH$ind$EBLUP)
c(sae_FH$mse)
as.numeric(My_FH$MSE$MSE)

# The results using the code from the sae package and my code return the same
# results using REML estimation

My_FHCrude <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
               combined_data = combined_data, domains = "SmallArea",
               method = "sae_reml", interval = c(0, 1000),
               transformation = "log_crude")


My_FHSM <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
                    combined_data = combined_data, domains = "SmallArea",
                    method = "sae_reml", interval = c(0, 1000),
                    transformation = "log_SM")


My_FHBC2 <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
                 combined_data = combined_data, domains = "SmallArea",
                 method = "sae_reml", interval = c(0, 1000),
                 transformation = "log_BC2")


My_FHCrude$ind$EBLUP
My_FHSM$ind$EBLUP
My_FHBC2$ind$EBLUP

My_FHCrude$MSE$MSE
My_FHSM$MSE$MSE
My_FHBC2$MSE$MSE



# Adjusted methods
My_FHREML <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
                    combined_data = combined_data, domains = "SmallArea",
                    method = "nicola_reml", interval = c(0, 1000),
                    transformation = "no")


My_FHAMPL <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
                 combined_data = combined_data, domains = "SmallArea",
                 method = "AMPL", interval = c(0, 1000),
                 transformation = "no")


My_FHAMRL <- FH_eblup(formula = yi ~ MajorArea + newVar, vardir = "SD2",
                  combined_data = combined_data, domains = "SmallArea",
                  method = "AMRL", interval = c(0, 1000),
                  transformation = "no")


My_FHREML$sigmau2
My_FHAMPL$sigmau2
My_FHAMRL$sigmau2


My_FHREML$ind$EBLUP
My_FHAMPL$ind$EBLUP
My_FHAMRL$ind$EBLUP

My_FHREML$MSE$MSE
My_FHAMPL$MSE$MSE
My_FHAMRL$MSE$MSE

