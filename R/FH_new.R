#' Function for Fay-Herriot model
#'
#' This function conducts the estimation of a Fay-Herriot model.
#'
#' @param formula formula object.
#' @param vardir direct variance.
#' @param combined_data combined data.
#' @param domains domain level.
#' @param method method for the estimation of sigmau2.
#' @param interval interval for the estimation of sigmau2.
#' @param precision precision criteria for the estimation of sigmau2.
#' @param maxiter maximum of iterations for the estimation of sigmau2.
#' @param transformation choses type of transformation and back-transformation.
#' @return fitted FH model.
#' @import formula.tools
#' @importFrom stats median model.frame model.matrix model.response optimize
#' @importFrom stats pnorm rnorm
#' @export


FH_eblup <- function(formula, vardir, combined_data, domains = NULL, method,
                     transformation = NULL, interval = c(0, 1000), precision = 0.0001,
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
                               eblup = eblup, combined_data = combined_data,
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
