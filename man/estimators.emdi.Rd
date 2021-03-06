% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.emdi.R
\name{estimators.emdi}
\alias{estimators.emdi}
\title{Presents point, MSE and/or CV estimates of an emdiObject}
\usage{
\method{estimators}{emdi}(object, indicator = "all", MSE = FALSE,
  CV = FALSE, ...)
}
\arguments{
\item{object}{an object of type "emdi", representing point and,
if chosen, MSE estimates.}

\item{indicator}{optional character vector that selects which indicators
shall be returned: (i) all calculated indicators ("all");
(ii) each indicator name: "Mean", "Quantile_10", "Quantile_25", "Median",
"Quantile_75", "Quantile_90", "Head_Count",
"Poverty_Gap", "Gini", "Quintile_Share" or the function name/s of
"custom_indicator/s"; (iii) groups of indicators: "Quantiles", "Poverty",
"Inequality" or "Custom". If two of these groups are selected, only the first
one is returned. Defaults to "all". Note, additional custom indicators can be
defined as argument for model-based approaches (see also \code{\link{ebp}})
and do not appear in groups of indicators even though these might belong to
one of the groups.}

\item{MSE}{optional logical. If \code{TRUE}, MSE estimates for selected indicators
per domain are added to the data frame of point estimates. Defaults to
\code{FALSE}.}

\item{CV}{optional logical. If \code{TRUE}, coefficients of variation for selected
indicators per domain are added to the data frame of point estimates.
Defaults to \code{FALSE}.}

\item{...}{other parameters that can be passed to function \code{estimators}.}
}
\value{
an object of type "estimators.emdi" with point and/or MSE
estimates and/or calculated CV's per domain obtained from
\code{emdiObject$ind} and, if chosen, \code{emdiObject$MSE}. These objects
contain two elements, one data frame \code{ind} and a character naming the
indicator or indicator group \code{ind_name}.
}
\description{
Method \code{estimators.emdi} presents point and MSE estimates for regional
disaggregated indicators. Coefficients of variation are calculated
using these estimators. This method enables to select for which indicators
the estimates shall be returned. The returned object is suitable for
printing with the \code{print.estimators.emdi} method.
}
\examples{
\dontrun{
# Loading data - population and sample data
data("eusilcA_pop")
data("eusilcA_smp")

# Generate emdi object with additional indicators; here via function ebp()
emdi_model <- ebp(fixed = eqIncome ~ gender + eqsize + cash +
self_empl + unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
fam_allow + house_allow + cap_inv + tax_adj, pop_data = eusilcA_pop,
pop_domains = "district", smp_data = eusilcA_smp, smp_domains = "district",
threshold = 11064.82, transformation = "box.cox",
L = 50, MSE = TRUE, B = 50, custom_indicator =
list(my_max = function(y, threshold){max(y)},
my_min = function(y, threshold){min(y)}), na.rm = TRUE, cpus = 1)

# Example 1: Choose Gini coefficient, MSE and CV
estimators(emdi_model, indicator = "Gini", MSE = TRUE, CV = TRUE)

# Example 2: Choose custom indicators without MSE and CV
estimators(emdi_model, indicator = "Custom")
}
}
\seealso{
\code{\link{emdiObject}}, \code{\link{direct}}, \code{\link{ebp}}
}
