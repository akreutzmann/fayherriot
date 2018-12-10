################################################################################
## Title: Estimation of the variance of the random effects
## - Model based simulation study
## Author: Sylvia Harmening


# clear workspace
rm(list = ls(all = TRUE))

# load required packages
library(fayherriot)

# Set seed
set.seed(-22123)

# Number of simulation runs
NoSim <- 50

# number of areas
#m = 15
#m = 45
m <- 100


# Set sigma_u to 1 (for both variance patterns)
sigma.u.sim <- 0.004

# Implement sigma_e for variance patterns a) and b)
#vardir  <- as.vector(rep(0, m))
#vardir <- ((0.00768 * (1:m - 1)) / (m - 1)) + 0.001


eff_smpsize <- ((150 * (1:m - 1)) / (m - 1)) + 10
vardir <- 1/ (4 * eff_smpsize)

# pattern a)
#pat1=4
#pat2=0.6
#pat3=0.5
#pat4=0.4
#pat5=0.1

#pattern b)
#pat1 <- 3.5
#pat2 <- 3
#pat3 <- 2.5
#pat4 <- 2
#pat5 <- 1.5

# for m = 15
#m1=3
#m2=6
#m3=9
#m4=12

# for m = 45
#m1=9
#m2=18
#m3=27
#m4=36

# for m = 100
#m1 <- 20
#m2 <- 40
#m3 <- 60
#m4 <- 80

#vardir[1:m1] <- pat1
#vardir[(m1 + 1):m2] <- pat2
#vardir[(m2 + 1):m3] <- pat3
#vardir[(m3 + 1):m4] <- pat4
#vardir[(m4 + 1):m] <- pat5



# General notation:
# saeREML REML variance estimation like in sae package
# REML REML variance estimation how Nicola implemented it
# pml profile ML variance estimation
# AMRL adjusted REML (Li and Lahiri)
# AMPL adjusted ML (Li and Lahiri)
# AMRL_YL adjusted REML (Yoshimori and Lahiri)
# AMPL_YL adjusted ML (Yoshimori and Lahiri)
# True True population values

# Create matrices for results

# Point estimates

True.ratio	<- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

sigmau2.arcsinBoot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.arcsinBoot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.arcsinBoot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
Li.arcsinBoot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
Ui.arcsinBoot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

sigmau2.arcsinJack <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.arcsinJack <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.arcsinJack <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.arcsinJack <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))




# Gamma estimates
gamma.arcsin <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
True.gamma <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

# MSE estimates
MSE.boot <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.jackknife <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

# Percentages of zero estimates of sigma_u
sigmau2.arcsin <- matrix(0, NoSim, m)


# beta is set 0 -> common mean model
# beta <- 0
beta <- c(0.001, 1)

below_zero <- NULL
above_one <- NULL

# Simulation (10000 runs)

for (h in 1:NoSim) {

  cat(date(),"Iteration number",h,"starting","\n",fill = T)

  #Xpop <- matrix(1, m, 1) # auxiliary variable (m * 1)
  Xpop <- matrix(c(rep(1, m), ((1:m) / (2 * m)) + 0.25), nrow = m, ncol = 2)
  # random effects randomly drawn from a normal distribution
  u <- matrix(rnorm(m, 0, sqrt(sigma.u.sim)), m, 1)
  # sampling errors randomly drawn from a normal distribution
  e <- rnorm(m, 0, sqrt(vardir))

  theta <- Xpop%*%beta + u # true means
  True.ratio[h,] <- theta
  yi <- theta  + e # direct estimates
  True.gamma[h,] <- vardir / (sigma.u.sim + vardir) # true gammas

  below_zero <- c(below_zero, any(theta < 0), any(yi < 0))
  above_one <- c(above_one, any(theta > 1), any(yi > 1))


  dom <- as.vector(1:m) # domain indicator for population
  area <- as.vector(1:m) # domain indicator for sample
  admin_data <- data.frame(Xpop, dom) # data frame with covariates and domains
  sd2 <- vardir
  # data frame with direct estimates, variances and domains
  sample_data <- data.frame(yi, area, sd2, eff_smpsize)

  # combine admin und sample data
  data_comb <- combine_data(pop_data = admin_data,"dom",
                            smp_data = sample_data, "area", "sd2")

  # Estimate FH models

  # REML variance estimation Nicolas implementation
  Boot <- fh(fixed = yi ~ 1, vardir = "sd2",
             combined_data = data_comb, domains = "area",
             method = "reml", interval = c(0, 100),
             eff_smpsize = "eff_smpsize",
             transformation = "arcsin_boot", alpha = 0.05)
  sigmau2.arcsinBoot[h,] = Boot$sigmau2
  gamma.arcsinBoot[h,] <- vardir / (sigmau2.arcsinBoot[h,] + vardir)
  EBLUP.arcsinBoot[h,] <- Boot$ind$EBLUP
  Li.arcsinBoot[h,]		<- Boot$MSE$Li
  Ui.arcsinBoot[h,]		<- Boot$MSE$Ui

  # REML variance estimation Nicolas implementation
  Jack <- fh(fixed = yi ~ 1, vardir = "sd2",
             combined_data = data_comb, domains = "area",
             method = "reml", interval = c(0, 100),
             eff_smpsize = "eff_smpsize",
             transformation = "arcsin_jack", alpha = 0.05)
  sigmau2.arcsinJack[h,] = Jack$sigmau2
  gamma.arcsinJack[h,] <- vardir / (sigmau2.arcsinJack[h,] + vardir)
  EBLUP.arcsinJack[h,] <- Jack$ind$EBLUP
  MSE.arcsinJack[h,]		<- Jack$MSE$MSE

}





################################################################################
# Comparison of different estimators of sigmau2
source("./Spielplatz/QualityMeasure.R")

True.sigma <- matrix(sigma.u.sim, NoSim, m)

#RB_sigmaREML <- mean((sigmau2.REML[, 1] - True.sigma[ ,1]) / True.sigma[ ,1])

quality_sigmaArcsinBoot <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.arcsinBoot[, 1]),
                                    MSETF = FALSE)
quality_sigmaArcsinBoot$RB * 100

quality_sigmaArcsinJack <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                          Est.mean = t(sigmau2.arcsinJack[, 1]),
                                          MSETF = FALSE)
quality_sigmaArcsinJack$RB * 100



# Comparison of different estimators of Gamma

# Relative Bias (%)


# REML
quality_gammaBoot <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.arcsinBoot),
                               MSETF = FALSE)
summary(quality_gammaBoot$RB * 100)

# ML
quality_gammaJack <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.arcsinJack),
                              MSETF = FALSE)
summary(quality_gammaJack$RB * 100)



################################################################################
# Comparison of different MSE estimators


quality_arcsinBoot <- QualityMeasure(True.mean = t(True.ratio), Est.mean = t(EBLUP.arcsinBoot),
                               MSE = t(MSE.REML), MSETF = FALSE)
summary(quality_arcsinBoot$RB * 100)

quality_arcsinJack <- QualityMeasure(True.mean = t(True.ratio), Est.mean = t(EBLUP.arcsinJack),
                              MSE = t(MSE.arcsinJack), MSETF = TRUE)
summary(quality_arcsinJack$RB * 100)
summary(quality_arcsinJack$RB_RMSE * 100)


