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


# Set parameters
maxiter <- 100
precision <- 0.0001
interval <- c(0, 1000)
vardir <- "SD2"
formula <- yi ~ MajorArea_new + SecondVar

milk$newVar <- rnorm(43) * milk$MajorArea

milk_smp <- milk[1:40,]
milk_admin <- milk[,c("SmallArea","MajorArea", "newVar")]
names(milk_admin) <- c("Domain", "MajorArea_new", "SecondVar") 

smp_domains <- "SmallArea"
admin_domains <- "Domain"

milk_smp$SD2 <- milk_smp$SD^2
smp_data <- milk_smp
admin_data <- milk_admin

combined_data <- combine_data(admin_data, "Domain", smp_data, "SmallArea",
                              "SD2")
combined_data <- combined_data[-which(is.na(combined_data$SmallArea)),]


# Have example with sae code
sae_FH <- mseFH(yi ~ MajorArea + newVar, SD2, data = combined_data)

My_FH <- FH_AK(formula = yi ~ MajorArea + newVar, vardir = "SD2", 
      combined_data = combined_data, domains = "SmallArea", 
      method = "sae_reml", interval = c(0, 1000), back_transformation = NULL)

all.equal(as.numeric(sae_FH$est$eblup), as.numeric(My_FH$ind$EBLUP))


combined_data$logyi <- log(combined_data$yi)
combined_data$logSD2 <- (1 / combined_data$yi)^2 * combined_data$SD2

My_FHNaive <- FH_AK(formula = logyi ~ MajorArea + newVar, vardir = "logSD2", 
               combined_data = combined_data, domains = "SmallArea", 
               method = "sae_reml", interval = c(0, 1000), 
               back_transformation = "naive")


My_FHSM <- FH_AK(formula = logyi ~ MajorArea + newVar, vardir = "logSD2", 
                    combined_data = combined_data, domains = "SmallArea", 
                    method = "sae_reml", interval = c(0, 1000), 
                    back_transformation = "SM")


My_FHBC2 <- FH_AK(formula = logyi ~ MajorArea + newVar, vardir = "logSD2", 
                 combined_data = combined_data, domains = "SmallArea", 
                 method = "sae_reml", interval = c(0, 1000), 
                 back_transformation = "BC2")


My_FHNaive$ind$EBLUP
My_FHSM$ind$EBLUP
My_FHBC2$ind$EBLUP




