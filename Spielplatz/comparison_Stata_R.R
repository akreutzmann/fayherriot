# load packages
library(readstata13)
library(sae)

################################################################################
##### Bundeslandebene
# load data
estimationData <- read.dta13("H:/STATA SAE/OUT/Bundeslaender/estimationData.dta")

# Var = SE^2
estimationData$SD2 <- (estimationData$seRGroup_i1hinc)^2

# fehlende Werte?
summary(is.na(estimationData))

# Have example with sae code
sae_FH <- mseFH(pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice,
                SD2, data = estimationData)

My_FH <- FH_AK(formula = pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice,
               vardir = "SD2",
               combined_data = estimationData, domains = "kkzID",
               method = "nicola_reml", interval = c(0, 100000))


My_FH$MSE

all.equal(as.numeric(My_FH$ind$EBLUP), as.numeric(sae_FH$est$eblup))
all.equal(as.numeric(My_FH$MSE$PR_MSE), as.numeric(sae_FH$mse))


# Check model selection
step(lm(pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice,
        data = estimationData))

summary(lm(pointDirect_i1hinc ~ unemplQuota + GDPperCapita,
           data = estimationData))

# Have example with sae code
sae_FH <- mseFH(pointDirect_i1hinc ~ unemplQuota + GDPperCapita,
                SD2, data = estimationData)

My_FH <- FH_AK(formula = pointDirect_i1hinc ~ unemplQuota + GDPperCapita,
               vardir = "SD2",
               combined_data = estimationData, domains = "kkzID",
               method = "nicola_reml", interval = c(0, 100000))


My_FH$MSE

all.equal(as.numeric(My_FH$ind$EBLUP), as.numeric(sae_FH$est$eblup))
all.equal(as.numeric(My_FH$MSE$PR_MSE), as.numeric(sae_FH$mse))


# Check Normality
qqnorm(My_FH$random_effects)
qqline(My_FH$random_effects)
shapiro.test(My_FH$random_effects)
plot(density(My_FH$random_effects))


qqnorm(My_FH$real_residuals)
qqline(My_FH$real_residuals)
shapiro.test(My_FH$real_residuals)
plot(density(My_FH$real_residuals))


# CV plots
CV_plots(My_FH)


################################################################################
##### ROR - Ebene
# load data
estimationDataROR <- read.dta13("H:/STATA SAE/OUT/estimationDataROR.dta")

# Var = SE^2
estimationDataROR$SD2 <- (estimationDataROR$seRGroup_i1hinc)^2

# fehlende Werte?
summary(is.na(estimationDataROR))

estimationDataRORsae <- estimationDataROR[complete.cases(estimationDataROR),]
estimationDataRORsae$SD2 <- (estimationDataRORsae$seRGroup_i1hinc)^2

# Have example with sae code
sae_FH <- mseFH(pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice, SD2,
                data = estimationDataRORsae)

My_FH <- FH_AK(formula = pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice,
               vardir = "SD2",
               combined_data = estimationDataROR, domains = "kkzID",
               method = "nicola_reml", interval = c(0, 100000))

sae_FH$mse
My_FH$MSE

all.equal(as.numeric(My_FH$ind$EBLUP), as.numeric(sae_FH$est$eblup))
all.equal(as.numeric(My_FH$MSE$PR_MSE), as.numeric(sae_FH$mse))

# Check model
step(lm(pointDirect_i1hinc ~ unemplQuota + GDPperCapita + landPrice,
        data = estimationData))

summary(lm(pointDirect_i1hinc ~ unemplQuota + GDPperCapita,
           data = estimationData))

# Have example with sae code
sae_FH <- mseFH(pointDirect_i1hinc ~ unemplQuota + GDPperCapita, SD2,
                data = estimationDataRORsae)

My_FH <- FH_AK(formula = pointDirect_i1hinc ~ unemplQuota + GDPperCapita,
               vardir = "SD2",
               combined_data = estimationDataROR, domains = "kkzID",
               method = "nicola_reml", interval = c(0, 100000))

# Check normality
qqnorm(My_FH$random_effects)
qqline(My_FH$random_effects)
shapiro.test(My_FH$random_effects)
plot(density(My_FH$random_effects))


qqnorm(My_FH$real_residuals)
qqline(My_FH$real_residuals)
shapiro.test(My_FH$real_residuals)
plot(density(My_FH$real_residuals))


# CV plots
CV_plots(My_FH)

