#' Function for Fay-Herriot model
#'
#' This function conducts the estimation of a Fay-Herriot model.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the nested error linear regression model with the
#' dependent variable on the left of a ~ operator and the explanatory
#' variables on the right, separated by + operators.
#' @param vardir a character string indicating the name of the variable containing
#' the domain-specific sampling variances of the direct estimators that are
#' included in \cr \code{combined_data}.
#' @param combined_data a data set containing the direct estimates,
#' the sampling variances, the explanatory variables and the domains.
#' @param domains a character string indicating the domain variable that is
#' included in \cr \code{combined_data}. If \code{NULL}, the domains are numbered
#' consecutively.
#' @param method a character string describing the method for the estimation of
#' the variance of the random effects. Methods that can be chosen
#' (i) restricted maximum likelihood (REML) method ("\code{reml}"),
#' (iii) maximum likelihood method ("\code{ml}"),
#' (iv) adjusted REML following \cite{Li and Lahiri (2010)} ("\code{amrl}"),
#' (v) adjusted ML following \cite{Li and Lahiri (2010)} ("\code{ampl}"),
#' (vi) adjusted REML following \cite{Yoshimori and Lahiri (2014)}
#' ("\code{amrl_yl}"),
#' (vii) adjusted ML following \cite{Yoshimori and Lahiri (2014)}
#' ("\code{ampl_yl}").
#'  Defaults to "\code{reml}".
#' @param MSE a character string determining the estimation of the MSE estimates.
#' Analytical or jackknife MSEs can be chosen.
#' @param interval interval for the estimation of sigmau2.
#' @param precision precision criteria for the estimation of sigmau2.
#' @param maxiter maximum of iterations for the estimation of sigmau2.
#' @param transformation choses type of transformation and back-transformation.
#' @param eff_smpsize Effective sample size.
#' @return fitted FH model.
#' @import formula.tools
#' @importFrom stats median model.frame model.matrix model.response optimize
#' @importFrom stats pnorm rnorm
#' @export


fh <- function(fixed, vardir, combined_data, domains = NULL, method = "reml",
                     MSE = "analytical", transformation = "no", eff_smpsize = NULL,
                     interval = c(0, 1000),
                     precision = 0.0001, maxiter = 100, alpha = alpha) {


  # Notational framework
  framework <- framework_FH(combined_data = combined_data, fixed = fixed,
                            vardir = vardir, domains = domains,
                            transformation = transformation,
                            eff_smpsize = eff_smpsize)


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
    if (MSE == "analytical") {
      MSE_data <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                                 combined_data = combined_data,
                                 method = method)
    } else if (MSE == "jackknife") {
      MSE_data <- jiang_jackknife(framework = framework,
                                  combined_data = combined_data, sigmau2 = sigmau2,
                                  vardir = vardir, eblup = eblup, transformation = transformation,
                                  method = method, interval = interval)
    }



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
                                   method = method,
                                   precision = precision, maxiter = maxiter,
                                   interval = interval, alpha = alpha)

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
