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
m = 45
#m <- 100


# Set sigma_u to 1 (for both variance patterns)
#sigma.u.sim <- 1
sigma.u.sim <- 2

# Implement sigma_e for variance patterns a) and b)
vardir  <- as.vector(rep(0, m))
#vardir <- ((4 * (1:m - 1)) / (m - 1)) + 2

# pattern a)
#pat1=4
#pat2=0.6
#pat3=0.5
#pat4=0.4
#pat5=0.1

#pattern b)
pat1 <- 3.5
pat2 <- 3
pat3 <- 2.5
pat4 <- 2
pat5 <- 1.5

# for m = 15
#m1=3
#m2=6
#m3=9
#m4=12

# for m = 45
m1=9
m2=18
m3=27
m4=36

# for m = 100
#m1 <- 20
#m2 <- 40
#m3 <- 60
#m4 <- 80

vardir[1:m1] <- pat1
vardir[(m1 + 1):m2] <- pat2
vardir[(m2 + 1):m3] <- pat3
vardir[(m3 + 1):m4] <- pat4
vardir[(m4 + 1):m] <- pat5



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
EBLUP.REML <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.pml <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.AMRL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.AMPL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.AMRL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
EBLUP.AMPL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
True.mean	<- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

# Gamma estimates
gamma.pml <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.REML <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.AMRL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.AMPL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.AMRL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
gamma.AMPL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
True.gamma <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

# MSE estimates
MSE.REML <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.pml <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.AMRL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.AMPL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.AMRL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))
MSE.AMPL_YL <- array(c(rep(NA, NoSim * m)), dim = c(NoSim, m))

# Percentages of zero estimates of sigma_u
sigmau2.REML <- matrix(0, NoSim, m)
sigmau2.pml <- matrix(0, NoSim, m)
sigmau2.AMRL <- matrix(0, NoSim, m)
sigmau2.AMPL <- matrix(0, NoSim, m)
sigmau2.AMRL_YL <- matrix(0, NoSim, m)
sigmau2.AMPL_YL <- matrix(0, NoSim, m)

# beta is set 0 -> common mean model
beta <- 0
#beta <- c(100, 5)

# Simulation (10000 runs)

for (h in 1:NoSim) {

  cat(date(),"Iteration number",h,"starting","\n",fill = T)

  Xpop <- matrix(1, m, 1) # auxiliary variable (m * 1)
  #Xpop <- matrix(c(rep(1, m), ((1:m) / (2 * m)) + 1), nrow = m, ncol = 2)
  # random effects randomly drawn from a normal distribution
  u <- matrix(rnorm(m, 0, sqrt(sigma.u.sim)), m, 1)
  # sampling errors randomly drawn from a normal distribution
  e <- rnorm(m, 0, sqrt(vardir))

  theta <- Xpop%*%beta + u # true means
  True.mean[h,] <- theta
  yi <- theta  + e # direct estimates
  True.gamma[h,] <- vardir / (sigma.u.sim + vardir) # true gammas


  dom <- as.vector(1:m) # domain indicator for population
  area <- as.vector(1:m) # domain indicator for sample
  admin_data <- data.frame(Xpop, dom) # data frame with covariates and domains
  sd2 <- vardir
  # data frame with direct estimates, variances and domains
  sample_data <- data.frame(yi, area, sd2)

  # combine admin und sample data
  data_comb <- combine_data(pop_data = admin_data,"dom",
                            smp_data = sample_data, "area", "sd2")

  # Estimate FH models

  # REML variance estimation Nicolas implementation
  Reml <- fh(fixed = yi ~ 1, vardir = "sd2",
                   combined_data = data_comb, domains = "area",
                   method = "reml", interval = c(0, 100),
                   transformation = "no")
  sigmau2.REML[h,] = Reml$sigmau2
  gamma.REML[h,] <- vardir / (sigmau2.REML[h,] + vardir)
  EBLUP.REML[h,] <- Reml$ind$EBLUP
  MSE.REML[h,]		<- Reml$MSE$MSE

  # ML variance estimation
  pml <- fh(fixed = yi ~ 1, vardir = "sd2",
                      combined_data = data_comb, domains = "area",
                      method = "ml", interval = c(0, 100),
                      transformation = "no")
  sigmau2.pml[h,] = pml$sigmau2
  gamma.pml[h,] <- vardir / (sigmau2.pml[h,] + vardir)
  EBLUP.pml[h,] <- pml$ind$EBLUP
  MSE.pml[h,]		<- pml$MSE$MSE

  # adjusted REML variance estimation (Li and Lahiri)
  adj_AMRL <- fh(fixed = yi ~ 1, vardir = "sd2",
                   combined_data = data_comb, domains = "area",
                   method = "amrl", interval = c(0, 100),
                   transformation = "no")
  sigmau2.AMRL[h,] = adj_AMRL$sigmau2
  gamma.AMRL[h,] <- vardir / (sigmau2.AMRL[h,] + vardir)
  EBLUP.AMRL[h,] <- adj_AMRL$ind$EBLUP
  MSE.AMRL[h,]		<- adj_AMRL$MSE$MSE

  # adjusted ML variance estimation (Li and Lahiri)
  adj_AMPL <- fh(fixed = yi ~ 1, vardir = "sd2",
                   combined_data = data_comb, domains = "area",
                   method = "ampl", interval = c(0, 100),
                   transformation = "no")
  sigmau2.AMPL[h,] = adj_AMPL$sigmau2
  gamma.AMPL[h,] <- vardir / (sigmau2.AMPL[h,] + vardir)
  EBLUP.AMPL[h,] <- adj_AMPL$ind$EBLUP
  MSE.AMPL[h,]		<- adj_AMPL$MSE$MSE

  # adjusted REML variance estimation (Yoshimori and Lahiri)
  YL_AMRL <- fh(fixed = yi ~ 1, vardir = "sd2",
                    combined_data = data_comb, domains = "area",
                    method = "amrl_yl", interval = c(0, 100),
                    transformation = "no")
  sigmau2.AMRL_YL[h,] = YL_AMRL$sigmau2
  gamma.AMRL_YL[h,] <- vardir / (sigmau2.AMRL_YL[h,] + vardir)
  EBLUP.AMRL_YL[h,] <- YL_AMRL$ind$EBLUP
  MSE.AMRL_YL[h,]		<- YL_AMRL$MSE$MSE

  # adjusted ML variance estimation (Yoshimori and Lahiri)
  YL_AMPL <- fh(fixed = yi ~ 1, vardir = "sd2",
                   combined_data = data_comb, domains = "area",
                   method = "ampl_yl", interval = c(0, 100),
                   transformation = "no")
  sigmau2.AMPL_YL[h,] = YL_AMPL$sigmau2
  gamma.AMPL_YL[h,] <- vardir / (sigmau2.AMPL_YL[h,] + vardir)
  EBLUP.AMPL_YL[h,] <- YL_AMPL$ind$EBLUP
  MSE.AMPL_YL[h,]		<- YL_AMPL$MSE$MSE

}



################################################################################
### Results

################################################################################
# Percentage of zero estimates of sigma_u
sigmau2.REML.zero <- (length(sigmau2.REML[sigmau2.REML <= 1e-4]) /
                        length(sigmau2.REML))*100
sigmau2.pml.zero <- (length(sigmau2.pml[sigmau2.pml <= 1e-4]) /
                       length(sigmau2.pml))*100
sigmau2.AMRL.zero <- (length(sigmau2.AMRL[sigmau2.AMRL <= 1e-12]) /
                        length(sigmau2.AMRL))*100
sigmau2.AMPL.zero <- (length(sigmau2.AMPL[sigmau2.AMPL <= 1e-12]) /
                        length(sigmau2.AMPL))*100
sigmau2.AMRL_YL.zero <- (length(sigmau2.AMRL_YL[sigmau2.AMRL_YL <= 1e-12]) /
                           length(sigmau2.AMRL_YL))*100
sigmau2.AMPL_YL.zero <- (length(sigmau2.AMPL_YL[sigmau2.AMPL_YL <= 1e-12]) /
                           length(sigmau2.AMPL_YL))*100

################################################################################
# Comparison of different estimators of sigmau2
source("./Spielplatz/QualityMeasure.R")

True.sigma <- matrix(sigma.u.sim, NoSim, m)

#RB_sigmaREML <- mean((sigmau2.REML[, 1] - True.sigma[ ,1]) / True.sigma[ ,1])

quality_sigmaREML <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.REML[, 1]),
                                    MSETF = FALSE)
quality_sigmaREML$RB * 100


quality_sigmapml <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.pml[, 1]),
                                    MSETF = FALSE)
quality_sigmapml$RB * 100

quality_sigmaAMRL <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.AMRL[, 1]),
                                    MSETF = FALSE)
quality_sigmaAMRL$RB * 100

quality_sigmaAMPL <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.AMPL[, 1]),
                                    MSETF = FALSE)
quality_sigmaAMPL$RB * 100


quality_sigmaAMRL_YL <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.AMRL_YL[, 1]),
                                    MSETF = FALSE)
quality_sigmaAMRL_YL$RB * 100

quality_sigmaAMPL_YL <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(sigmau2.AMPL_YL[, 1]),
                                    MSETF = FALSE)
quality_sigmaAMPL_YL$RB * 100


# Comparison of different estimators of Gamma

# Relative Bias (%)


# REML
quality_REML <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.REML),
                               MSETF = FALSE)
summary(quality_REML$RB * 100)

# ML
quality_pml <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.pml),
                               MSETF = FALSE)
summary(quality_pml$RB * 100)

# AMRL
quality_AMRL <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.AMRL),
                              MSETF = FALSE)
summary(quality_AMRL$RB * 100)

# AMPL
quality_AMPL <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.AMPL),
                              MSETF = FALSE)
summary(quality_AMPL$RB * 100)

# AMRL_YL
quality_AMRL_YL <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.AMRL_YL),
                               MSETF = FALSE)
summary(quality_AMRL_YL$RB * 100)

# AMPL
quality_AMPL_YL <- QualityMeasure(True.mean = t(True.gamma), Est.mean = t(gamma.AMPL_YL),
                               MSETF = FALSE)
summary(quality_AMPL_YL$RB * 100)



################################################################################
# Comparison of different MSE estimators


quality_REML <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.REML),
                               MSE = t(MSE.REML), MSETF = FALSE)
summary(quality_REML$RB * 100)
summary(quality_REML$RB_RMSE * 100)

plot(quality_REML$True.RMSE, type = "l")
lines(quality_REML$Est.RMSE, col = "red")


quality_REMLPoint <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.REML),
                               MSE = t(MSE.REML), MSETF = FALSE)
plot(quality_REMLPoint$True, type = "l")
lines(quality_REMLPoint$Dir, col = "red")
summary(quality_REMLPoint$Bias)
summary(quality_REMLPoint$True)


quality_pml <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.pml),
                               MSE = t(MSE.pml), MSETF = TRUE)
summary(quality_pml$RB * 100)
summary(quality_pml$RB_RMSE * 100)


quality_AMRL <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.AMRL),
                               MSE = t(MSE.AMRL), MSETF = TRUE)
summary(quality_AMRL$RB * 100)
summary(quality_AMRL$RB_RMSE * 100)

quality_AMPL <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.AMPL),
                               MSE = t(MSE.AMPL), MSETF = TRUE)
summary(quality_AMPL$RB * 100)
summary(quality_AMPL$RB_RMSE * 100)


quality_AMRL_YL <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.AMRL_YL),
                               MSE = t(MSE.AMRL_YL), MSETF = TRUE)
summary(quality_AMRL_YL$RB * 100)
summary(quality_AMRL_YL$RB_RMSE * 100)

quality_AMPL_YL <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(EBLUP.AMPL_YL),
                               MSE = t(MSE.AMPL_YL), MSETF = TRUE)
summary(quality_AMPL_YL$RB * 100)
summary(quality_AMPL_YL$RB_RMSE * 100)






