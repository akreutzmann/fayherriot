################################################################################
# This R script is used to try all the possibilities implemented in fayherriot.
################################################################################

# Load aggregated data ---------------------------------------------------------
data("eusilcA_popAgg")
data("eusilcA_smpAgg")

# Combine sample and population data -------------------------------------------
combined_data <- combine_data(pop_data = eusilcA_popAgg, pop_domains = "Domain",
                              smp_data = eusilcA_smpAgg, smp_domains = "Domain")

# Estimation of EBLUP means without transformation -----------------------------

# REML
fh_reml <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "reml", interval = c(0, 100000000))
print(fh_reml)
summary(fh_reml)
plot(fh_reml)
library(emdi)
load_shapeaustria()
map_plot(fh_reml, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_reml, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_reml, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))

# ML
fh_ml <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "ml", interval = c(0, 100000000))
print(fh_ml)
summary(fh_ml)
plot(fh_ml)
map_plot(fh_ml, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_ml, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_ml, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))

# AMRL following Li and Lahiri
fh_amrl <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "amrl", interval = c(0, 100000000))
print(fh_amrl)
summary(fh_amrl)
plot(fh_amrl)
map_plot(fh_amrl, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_amrl, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_amrl, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))

# AMPL following Li and Lahiri
fh_ampl <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "ampl", interval = c(0, 100000000))
print(fh_ampl)
summary(fh_ampl)
plot(fh_ampl)
map_plot(fh_ampl, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_ampl, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_ampl, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))


# AMRL following Yoshimori and Lahiri
fh_amrl_yl <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "amrl_yl", interval = c(0, 100000000))
print(fh_amrl_yl)
summary(fh_amrl_yl)
plot(fh_amrl_yl)
map_plot(fh_amrl_yl, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_amrl_yl, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_amrl_yl, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))


# AMPL following Yoshimori and Lahiri
fh_ampl_yl <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
              combined_data = combined_data, domains = "Domain",
              method = "ampl_yl", interval = c(0, 100000000))
print(fh_ampl_yl)
summary(fh_ampl_yl)
plot(fh_ampl_yl)
map_plot(fh_ampl_yl, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_ampl_yl, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_ampl_yl, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))


# Estimation of EBLUP means with log transformation -----------------------------

# REML and crude back-transformation
fh_reml_crude <- fh(fixed = Mean ~ eqsize + cash + self_empl, vardir = "Var_Mean",
                    combined_data = combined_data, domains = "Domain",
                    method = "reml", interval = c(0, 100000000), transformation = "log_crude")
print(fh_reml_crude)
summary(fh_reml_crude)
plot(fh_reml_crude)
library(emdi)
load_shapeaustria()
map_plot(fh_reml_crude, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
         indicator = "EBLUP", map_dom_id = "PB")
compare_fh(model = fh_reml_crude, label = "orig",
           color = c("blue", "lightblue3"),
           shape = c(16, 16), line_type = c("solid", "solid"),
           gg_theme = NULL)
CV_plots(fh_reml_crude, label_direct = "Direct", label_FH = "FH", line = 10,
         colline = "red", color = c("blue", "lightblue3"))











