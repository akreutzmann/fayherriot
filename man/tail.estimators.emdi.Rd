% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.emdi.R
\name{tail.estimators.emdi}
\alias{tail.estimators.emdi}
\title{Returns the last part of predicted indicators and, if chosen, of MSE and
CV estimators.}
\usage{
\method{tail}{estimators.emdi}(x, n = 6L, addrownums = NULL, ...)
}
\arguments{
\item{x}{an object of type "estimators.emdi", representing
point estimators and, if chosen, MSE and/or CV estimates for selected
indicators.}

\item{n}{a single integer. If positive, it determines the number of rows for
the data frame. If negative, all but the n first rows of
elements of the object.}

\item{addrownums}{if there are no row names, create them from the row numbers.}

\item{...}{arguments to be passed to or from other methods.}
}
\value{
Selected rows of the object of type "estimators.emdi".
}
\description{
Returns the last part of predicted indicators and, if chosen, of MSE and
CV estimators.
}
\examples{
\dontrun{
# Loading data - population and sample data
data("eusilcA_pop")
data("eusilcA_smp")

# Generate emdi object with deleting missing values; here via function ebp()
emdi_model <- ebp(fixed = eqIncome ~ gender + eqsize + cash +
self_empl + unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
fam_allow + house_allow + cap_inv + tax_adj,
pop_data = eusilcA_pop, pop_domains = "district",
smp_data = eusilcA_smp, smp_domains = "district",
na.rm = TRUE)

# Example: Choose last lines of the Gini coefficient, MSE and CV
tail(estimators(emdi_model, indicator = c("Gini", "Head_Count")))
}
}
\seealso{
\code{\link{estimators.emdi}}
}
