library(emdi)

data(eusilcA_smp)
data(eusilcA_pop)


# Get weighted direct estimates and the corresponding variance
direct_estim <- direct("eqIncome", smp_data = eusilcA_smp, smp_domains = "district",
       weights = "weight", var = TRUE)

# Reduce the data set to mean and headcount ratio
eusilcA_smpAgg <- direct_estim$ind[, 1:3]
eusilcA_smpAgg <- data.frame(eusilcA_smpAgg, direct_estim$MSE[, 2:3])
names(eusilcA_smpAgg) <- c("Domain", "Mean", "Head_Count", "Var_Mean", "Var_Head_Count")

# Add sample sizes for effective sample size
eusilcA_smpAgg$n <- table(direct_estim$framework$smp_domains_vec)


# Get covariate information from the population data
pop_eqsize <- direct("eqsize", smp_data = eusilcA_pop, smp_domains = "district")
pop_cash <- direct("cash", smp_data = eusilcA_pop, smp_domains = "district")
pop_self_empl <- direct("self_empl", smp_data = eusilcA_pop, smp_domains = "district")
pop_unempl_ben <- direct("unempl_ben", smp_data = eusilcA_pop, smp_domains = "district")
pop_age_ben <- direct("age_ben", smp_data = eusilcA_pop, smp_domains = "district")
pop_surv_ben <- direct("surv_ben", smp_data = eusilcA_pop, smp_domains = "district")
pop_sick_ben <- direct("sick_ben", smp_data = eusilcA_pop, smp_domains = "district")
pop_dis_ben <- direct("dis_ben", smp_data = eusilcA_pop, smp_domains = "district")
pop_rent <- direct("rent", smp_data = eusilcA_pop, smp_domains = "district")
pop_fam_allow <- direct("fam_allow", smp_data = eusilcA_pop, smp_domains = "district")
pop_house_allow <- direct("house_allow", smp_data = eusilcA_pop, smp_domains = "district")
pop_cap_inv <- direct("cap_inv", smp_data = eusilcA_pop, smp_domains = "district")
pop_tax_adj <- direct("tax_adj", smp_data = eusilcA_pop, smp_domains = "district")

eusilcA_popAgg <- data.frame(Domain = pop_eqsize$ind$Domain,
                             eqsize = pop_eqsize$ind$Mean,
                             cash = pop_cash$ind$Mean,
                             self_empl = pop_self_empl$ind$Mean,
                             unempl_ben = pop_unempl_ben$ind$Mean,
                             age_ben = pop_age_ben$ind$Mean,
                             surv_ben = pop_surv_ben$ind$Mean,
                             sick_ben = pop_sick_ben$ind$Mean,
                             dis_ben = pop_dis_ben$ind$Mean,
                             rent = pop_rent$ind$Mean,
                             fam_allow = pop_fam_allow$ind$Mean,
                             house_allow = pop_house_allow$ind$Mean,
                             cap_inv = pop_cap_inv$ind$Mean,
                             tax_adj = pop_tax_adj$ind$Mean)

save(eusilcA_smpAgg, file = "./data/eusilcA_smpAgg.rda")
save(eusilcA_popAgg, file = "./data/eusilcA_popAgg.rda")
