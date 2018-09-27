#' EBLUPs and their corresponding MSE estimates for the Fay-Herriot model
#'
#' This function conducts the estimation of EBLUPS and MSEs based on a Fay-Herriot
#' model. Extended fitting methods can be chosen. And different back-transformations
#' of log transformed Fay-Herriot models can be selected.
#'
#' @param formula a formula object of the form: target variable ~ explanatory
#' variable 1 + explanatory variable 2 + ... .
#' @param vardir a character string indicating the name of the variable containing
#' the domain-specific sampling variances of the direct estimators that are
#' included in \code{combined_data} set.
#' @param combined_data a combined data set containing the direct estimates,
#' the sampling variances, the explanatory variables and the domains.
#' @param domains a character string indicating the domain variable that is
#' included in \code{combined_data}. If \code{NULL}, the domains are numbered
#' consecutively.
#' @param method a character string describing the method for the estimation of
#' the variance of the random effects. Methods that can be chosen (i) restricted
#' maximum likelihood (REML) method following the \pkg{sae} package ("\code{sae_reml}"),
#' (ii) restricted maximum likelihood method ("\code{nicola_reml}"), (iii) maximum
#' likelihood method ("\code{ML}"), (iv) adjusted REML following
#' \cite{Li and Lahiri (2010)} ("\code{AMRL}"), adjusted ML following
#' \cite{Li and Lahiri (2010)} ("\code{AMPL}"), (v) adjusted REML following
#' \cite{Yoshimori and Lahiri (2014)} ("\code{AMRL_YL}"), (vi) adjusted ML
#' following \cite{Yoshimori and Lahiri (2014)} ("\code{AMRL_YL}"). Defaults to
#' "\code{nicola_reml}".
#' @param transformation a character string. No transformation and log
#' transformation of the dependent variable with different back transformations
#' can be chosen (i) no transformation ("\code{no}"), (ii) crude back-transformation
#' ("\code{log_crude}"), (iii) back transformation following
#' \cite{Slud and Maiti (2006)} ("\code{log_SM}"),
#' (iv) back transformation following \cite{Chandra et al. (2017)} (Taylor series
#' approximation bias correction) ("\code{log_BC2}"). Defaults to "\code{no}".
#' @param interval a numeric vector containing a lower and upper limit for the
#' estimation of \code{sigmau2}. Default is set to \code{c(0,1000)}. In some cases
#' it may be more suitbale to choose a larger interval.
#' @param precision a precision criteria for the estimation of \code{sigmau2}.
#' Defaults to \code{0.0001}.
#' @param maxiter maximum of iterations for the estimation of \code{sigmau2}.
#' Defaults to \code{100}.
#' @details When choosing the back transformation "\code{log_SM}" \code{log_BC2}",
#' the variance estimation method must set to "\code{ML}" and the EBLUP and MSE
#' estimates are only returned for in-sample domains.
#' @return The function \code{FH_eblup} returns a list containing the following
#' objects:
#' \item{ind}{a data frame containing the domains, the direct estimates, the EBLUP
#' estimates and a variable indicating wheter the domain is an in- or out-of-sample
#' domain (0/1).}
#' \item{MSE}{a data frame containing the domains, the sampling variances, the MSE
#' estimates and a variable indicating wheter the domain is an in- or out-of-sample
#' domain (0/1).}
#' \item{method}{applied variance estimation method ("\code{sae_reml}",
#' "\code{nicola_reml}", "\code{ML}", "\code{AMRL}", "\code{AMPL}", "\code{AMRL_YL}"
#' or "\code{AMPL_YL}").}
#' \item{MSE_method}{applied MSE estimation method ("\code{Prasad_Rao}",
#' "\code{Datta_Lahiri}", "\code{AMRL_corrected}", "\code{AMPL_corrected}",
#' "\code{AMRL_YL_corrected}","\code{AMPL_YL_corrected}").}
#' \item{transformation}{a character string indicating the applied transformation
#' and back transformation method.}
#' \item{coefficients}{a data frame containing the estimated model coefficients
#' (betas), the standard errors,
#' \code{t}- and \code{p}-values of the explanatory variables.}
#' \item{sigmau2}{estimated variance of the random effects.}
#' \item{random_effects}{a matrix containing the random effects per domain.}
#' \item{real_residuals}{a matrix containing the real residuals per domain.}
#' \item{gamma}{a data frame containing the shrinkage factors per domain.}
#' \item{model_select}{a data frame containing different model selection and
#' accuracy criteria: loglikelihood,
#' Akaike information criterion (AIC), Bayesian information criterion (BIC),
#' R2 and adjusted R2 following \cite{Lahiri and Suntornchost (2015)}.}
#' @references
#' Chandra, H., Aditya, K. and Kumar, S. (2017), Small-area estimation under a
#' log-transformed area-level model, Journal of Statistical Theory and Practice
#' 12(3), 497-505. \cr \cr
#' Datta, G. S. and Lahiri, P. (2000), A unified measure of uncertainty of
#' estimated best linear unbiased predictors in small area estimation problems,
#' Statistica Sinica 10(2), 613-627. \cr \cr
#' Fay, R. E. and Herriot, R. A. (1979), Estimates of income for small places:
#' An application of James-Stein procedures to census data, Journal of the
#' American Statistical Association 74(366), 269-277. \cr \cr
#' Lahiri, P. and Suntornchost, J. (2015), Variable selection for linear mixed
#' models with applications in small area estimation, The Indian Journal of
#' Statistics 77-B(2), 312-320.
#' Li, H. and Lahiri, P. (2010), An adjusted maximum likelihood method for solving
#' small area estimation problems, Journal of Multivariate Analyis 101, 882-902. \cr \cr
#' Prasad, N. and Rao, J. (1990), The estimation of the mean squared error of
#' small-area estimation, Journal of the American Statistical Association 85(409),
#' 163-171. \cr \cr
#' Rao, J. N. K. and Molina, I. (2015), Small area estimation', New York: Wiley. \cr \cr
#' Slud, E. and Maiti, T. (2006), Mean-squared error estimation in transformed
#' Fay-Herriot models, Journal of the Royal Statistical Society: Series B 68(2),
#' 239-257.\cr \cr
#' Yoshimori, M. and Lahiri, P. (2014), A new adjusted maximum likelihood method
#' for the Fay-Herriot small area model, Journal of Multivariate Analysis
#' 124, 281-294.
#' @import formula.tools
#' @importFrom stats median model.frame model.matrix model.response optimize
#' @importFrom stats pnorm rnorm
#' @export


FH_eblup <- function(formula, vardir, combined_data, domains = NULL, method,
                     transformation = "no", interval = c(0, 1000), precision = 0.0001,
                     maxiter = 100) {


  # Notational framework
  framework <- framework_FH(combined_data = combined_data, formula = formula,
                            vardir = vardir, domains = domains,
                            transformation = transformation)


  # Estimate sigma u
  sigmau2 <- wrapper_estsigmau2(framework = framework, method = method,
                                precision = precision, maxiter = maxiter,
                                interval = interval)


  # Standard EBLUP
  eblup <- eblup_FH(framework = framework, sigmau2 = sigmau2,
                    combined_data = combined_data)


  # Criteria for model selection
  criteria <- model_select(framework = framework, sigmau2 = sigmau2,
                           real_res = eblup$real_res)


  if (transformation == "no") {

    # Analytical MSE
    MSE_data <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                               combined_data = combined_data,
                               method = method)


    # Shrinkage factor
    Gamma <- data.frame(Domain = framework$data[[framework$domains]],
                        Gamma = eblup$gamma)

    out <- list(ind = eblup$EBLUP_data,
                MSE = MSE_data$MSE_data,
                method = method,
                MSE_method = MSE_data$MSE_method,
                transformation = transformation,
                coefficients = eblup$coefficients,
                sigmau2 = sigmau2,
                random_effects = eblup$random_effects,
                real_residuals = eblup$real_res,
                gamma = Gamma,
                model_select = criteria)
  } else if (transformation != "no") {

    # Shrinkage factor
    Gamma <- data.frame(Domain = framework$data[[framework$domains]],
                        Gamma = eblup$gamma)

    # Back-transformation
    result_data <- backtransformed(framework = framework,
                                   sigmau2 = sigmau2, eblup = eblup,
                                   transformation = transformation,
                                   combined_data = combined_data,
                                   method = method)

    out <- list(ind = result_data$EBLUP_data,
                MSE = result_data$MSE_data,
                method = method,
                MSE_method = NULL,
                transformation = transformation,
                coefficients = eblup$coefficients,
                sigmau2 = sigmau2,
                random_effects = eblup$random_effects,
                real_residuals = eblup$real_res,
                gamma = Gamma,
                model_select = criteria)

  }

  class(out) <- "FH_eblup"

  return(out)

}
