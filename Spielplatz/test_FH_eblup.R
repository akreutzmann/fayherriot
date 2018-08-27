# Load data
library(readstata13)
milk <- read.dta13("H:/STATA_SAE/R Code/milk.dta")

# Install source package
install.packages("H:/fayherriot_0.1.0.tar.gz", repos = NULL, type = "source")
library(fayherriot)

# 1. REML, without transformation
FH_REML <- FH_eblup(formula = yi ~ as.factor(MajorArea), vardir = "SD2",
                 combined_data = milk, domains = "SmallArea",
                 method = "nicola_reml", interval = c(0, 1000))
FH_REML$sigmau2
FH_REML$coefficients

shapiro.test(FH_REML$real_residuals)
shapiro.test(FH_REML$random_effects)
FH_REML$model_select

# 2. AMPL, without transformation
FH_AMPL <- FH_eblup(formula = yi ~ as.factor(MajorArea), vardir = "SD2",
                 combined_data = milk, domains = "SmallArea",
                 method = "AMPL", interval = c(0, 1000))

FH_AMPL$sigmau2
FH_AMPL$coefficients

shapiro.test(FH_AMPL$real_residuals)
shapiro.test(FH_AMPL$random_effects)
FH_AMPL$model_select

# 3. AMRL_YL, without transformation
FH_AMRL_YL <- FH_eblup(formula = yi ~ as.factor(MajorArea), vardir = "SD2",
                    combined_data = milk, domains = "SmallArea",
                    method = "AMRL_YL", interval = c(0, 1000))

FH_AMRL_YL$sigmau2
FH_AMRL_YL$coefficients

shapiro.test(FH_AMRL_YL$real_residuals)
shapiro.test(FH_AMRL_YL$random_effects)
FH_AMRL_YL$model_select


# 4. ML, log transformation, SM bias correction
FH_SM <- FH_eblup(formula = yi ~ as.factor(MajorArea), vardir = "SD2",
               combined_data = milk, domains = "SmallArea",
               method = "ML", transformation = "log_SM",
               interval = c(0, 1000))

FH_SM$sigmau2
FH_SM$coefficients

shapiro.test(FH_SM$real_residuals)
shapiro.test(FH_SM$random_effects)
FH_SM$model_select


# 5. REML, log transformation, crude bias correction
FH_Crude <- FH_eblup(formula = yi ~ as.factor(MajorArea), vardir = "SD2",
                  combined_data = milk, domains = "SmallArea",
                  method = "nicola_reml", transformation = "log_crude",
                  interval = c(0, 1000))

FH_Crude$sigmau2
FH_Crude$coefficients

shapiro.test(FH_Crude$real_residuals)
shapiro.test(FH_Crude$random_effects)
FH_Crude$model_select

