################################################################################
## Title: Log transformed FH models - Model based simulation study
## Author: Sylvia Harmening


# clear workspace
rm(list = ls(all = TRUE))

# load required packages
library(fayherriot)

# Set seed
set.seed(100)

# Number of simulation runs
NoSim <- 50

# Number of areas
#m <- 30
#m = 50
m <- 100

# Sigma_u for variance patterns a) and b)
# pattern a)
sigma.u.sim <- 0.2

# pattern b)
# sigma.u.sim = 1

# Implement sigma_e for variance patterns a) and b)
vardir <- as.vector(rep(0, m))

# pattern a)
pat1 <- 0.5
pat2 <- .33
pat3 <- 0.25
pat4 <- 0.2
pat5 <- 0.17

# pattern b)
#pat1=1
#pat2=0.67
#pat3=0.5
#pat4=0.4
#pat5=0.33

# for m = 30
#m1=6
#m2=12
#m3=18
#m4=24

# for m = 50
#m1 = 10
#m2 = 20
#m3 = 30
#m4 = 40

# for m = 100
m1 <- 20
m2 <- 40
m3 <- 60
m4 <- 80

vardir[1:m1] <- pat1
vardir[(m1 + 1):m2] <- pat2
vardir[(m2 + 1):m3] <- pat3
vardir[(m3 + 1):m4] <- pat4
vardir[(m4 + 1):m] <- pat5

# General notation:
# SM back transformation following Slud and Maiti
# SMBC1 back transformation following Chandra et al. bias correction 1
# SMBC2 back transformation following Chandra et al. bias correction 2
# n/naiveML crude back transformation (ML variance estimation)
# nR/naiveREML crude back transformation (REML variance estimation)
# crude extreme naive back transformation
# True True population values


# Create matrices for results

# Point estimates
SMEBP <- matrix(NA, NoSim, m)
SMEBP.n	<- matrix(NA, NoSim, m)
SMEBP.nR <- matrix(NA, NoSim, m)
SMEBP.crude <- matrix(NA, NoSim, m)
True.mean	<- matrix(NA, NoSim, m)

# MSE estimates
MSE.SM <- matrix(NA, NoSim, m)
MSE.n <- matrix(NA, NoSim, m)
MSE.nR <- matrix(NA, NoSim, m)

# sigma_u estimates
estsigma2u_SM <- rep(NA, NoSim)
estsigma2u_n <- rep(NA, NoSim)
estsigma2u_nR <- rep(NA, NoSim)

# gamma estimates (SM)
gamma.SM <- matrix(0, NoSim, m)

# common mean model
alpha <- 0
Beta <- 1

for (h in 1:NoSim) {

  cat(date(),"Iteration number",h,"starting","\n",fill = T)

  Xpop <- matrix(1, m, 1) # auxiliary variable
  Xpop[,1] <- runif(m, 0, 1)
  Xpop[,1] <- log(Xpop)

  # random effects randomly drawn from a normal distribution
  u <- matrix(rnorm(m,0,sqrt(sigma.u.sim)), m, 1)
  # sampling errors randomly drawn from a normal distribution
  e <- rnorm(m, 0, sqrt(vardir))
  theta <- exp(alpha + Beta*Xpop[,1]  + u[,1])
  True.mean[h,] <- theta # true means
  log.yi <- alpha + Beta*Xpop[,1] + u + e # direct estimates


  dom <- as.vector(1:m) # domain indicator for population
  area <- as.vector(1:m) # domain indicator for sample
  admin_data <- data.frame(Xpop, dom) # data frame with covariates and domains
  sd2 <- vardir
  # data frame with direct estimates, variances and domains
  sample_data <- data.frame(log.yi = exp(log.yi), area, sd2 = sd2 / (1/exp(log.yi))^2)

  # combine admin und sample data
  data_comb <- combine_data(pop_data = admin_data,"dom",
                            smp_data = sample_data, "area", "sd2")

  # Estimate FH models

  # Slud and Maiti back transformation
  SM <- fh(fixed = log.yi ~ Xpop, vardir = "sd2", combined_data = data_comb,
                 domains = "area", method = "ml", interval = c(0, 1000),
                 transformation = "log_SM")
  estsigma2u_SM[h] <- SM$sigmau2
  SMEBP[h,] <- SM$ind$EBLUP
  MSE.SM[h,] <- SM$MSE$MSE
  gamma.SM[h,] <- SM$gamma$Gamma

  # Crude back transformation (ML variance estimation)
  SMnaive <- fh(fixed = log.yi ~ Xpop, vardir = "sd2",
                      combined_data = data_comb, domains = "area", method = "ml",
                      interval = c(0, 1000),
                      transformation = "log_crude")
  estsigma2u_n[h] <- SMnaive$sigmau2
  SMEBP.n[h,] <- SMnaive$ind$EBLUP
  MSE.n[h,] <- SMnaive$MSE$MSE

  # Crude back transformation (REML variance estimation)
  SMnaivereml <- fh(fixed = log.yi ~ Xpop, vardir = "sd2",
                      combined_data = data_comb, domains = "area",
                      method = "reml", interval = c(0, 1000),
                      transformation = "log_crude")
  estsigma2u_nR[h] <- SMnaivereml$sigmau2
  SMEBP.nR[h,] <- SMnaivereml$ind$EBLUP
  MSE.nR[h,] <- SMnaivereml$MSE$MSE
}


################################################################################
### Results
################################################################################
# Comparison of different EBLUPS


# Relative Bias (%)
source("./Spielplatz/QualityMeasure.R")


quality_SMEBP <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(SMEBP),
                               MSE = t(MSE.SM), MSETF = TRUE)
summary(quality_SMEBP$RB * 100)
summary(quality_SMEBP$RB_RMSE * 100)

quality_crudeREML <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(SMEBP.nR),
                                MSE = t(MSE.nR), MSETF = TRUE)
summary(quality_crudeREML$RB * 100)
summary(quality_crudeREML$RB_RMSE * 100)


quality_crudeML <- QualityMeasure(True.mean = t(True.mean), Est.mean = t(SMEBP.n),
                                    MSE = t(MSE.n), MSETF = TRUE)
summary(quality_crudeML$RB * 100)
summary(quality_crudeML$RB_RMSE * 100)


True.sigma <- matrix(sigma.u.sim, NoSim, m)
quality_sigmaSM <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                    Est.mean = t(estsigma2u_SM),
                                    MSETF = FALSE)
quality_sigmaSM$RB * 100

quality_sigmaCrudeREML <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                  Est.mean = t(estsigma2u_nR),
                                  MSETF = FALSE)
quality_sigmaCrudeREML$RB * 100

quality_sigmaCrudeML <- QualityMeasure(True.mean = t(True.sigma[, 1]),
                                         Est.mean = t(estsigma2u_n),
                                         MSETF = FALSE)
quality_sigmaCrudeML$RB * 100
