#' Summarizes an emdiObject
#'
#' Additional information about the data and model in small area estimation
#' methods and components of an emdi object are extracted. The returned object
#' is suitable for printing  with the \code{print.summary.emdi} method.
#' @param object an object of type "emdi", representing point and MSE
#' estimates. Objects differ depending on the estimation method: direct
#' vs. model-based.
#' @param ... additional arguments that are not used in this method.
#' @return an object of type "summary.emdi" with following
#' @importFrom moments skewness kurtosis
#' @export


summary.fh <- function (object) {

  N_domain_out_of_sample <- object$framework$M - object$framework$m
  N_domain_in_sample <- object$framework$m

  # Normaly checks for standardized realized residuals
  skewness_stdres <- skewness(object$model$std_real_residuals)
  kurtosis_stdres <- kurtosis(object$model$std_real_residuals)
    if (length(object$model$std_real_residuals) >= 3 & length(object$model$std_real_residuals) <
        5000) {
      shapiro_stdres_W <- shapiro.test(object$model$std_real_residuals)[[1]]
      shapiro_stdres_p <- shapiro.test(object$model$std_real_residuals)[[2]]
    }
    else {
      warning("Number of domains must be between 3 and 5000, otherwise the\n Shapiro-Wilk test is not applicable.")
      shapiro_stdres_W <- NA
      shapiro_stdres_p <- NA
    }

  # Normality checks for random effects
  skewness_random <- skewness(object$model$random_effects)
  kurtosis_random <- kurtosis(object$model$random_effects)
  if (length(object$model$random_effects) >= 3 & length(object$model$random_effects) <
      5000) {
    shapiro_random_W <- shapiro.test(object$model$random_effects)[[1]]
    shapiro_random_p <- shapiro.test(object$model$random_effects)[[2]]
  }
  else {
    shapiro_random_W <- NA
    shapiro_random_p <- NA
  }

  normality <- data.frame(Skewness = c(skewness_stdres, skewness_random),
                          Kurtosis = c(kurtosis_stdres, kurtosis_random),
                          Shapiro_W = c(shapiro_stdres_W,
                                        shapiro_random_W),
                          Shapiro_p = c(shapiro_stdres_p,
                                        shapiro_random_p),
                          row.names = c("Standardized_Residuals",
                                        "Random_effects"))

  if (object$transformation == "no") {
    transform_data <- NULL
  } else if (object$transformation == "log_crude") {
    transformation <- "log"
    backtransformation <- "crude"
    transform_data <- data.frame(Transformation  = transformation,
                                 Back_transformation = backtransformation,
                                 row.names       = ""
    )
  } else if (object$transformation == "log_SM") {
    transformation <- "log"
    backtransformation <- "Slud_Maiti"
    transform_data <- data.frame(Transformation  = transformation,
                                 Back_transformation = backtransformation,
                                 row.names       = ""
    )
  } else if (object$transformation == "arcsin") {
    transformation <- "arcsin"
    backtransformation <- "naive"
    transform_data <- data.frame(Transformation  = transformation,
                                 Back_transformation = backtransformation,
                                 row.names       = ""
    )
  }


  sum_FH_eblup <- list(call = object$call,
                       out_of_sample = N_domain_out_of_sample,
                       in_sample = N_domain_in_sample,
                       variance_estimation_method = object$method$method,
                       transformation = transform_data,
                       MSE_method = object$method$MSE_method,
                       coefficients = object$model$coefficients,
                       legend = object$legend,
                       sigmau2 = object$model$sigmau2,
                       model_select = object$model$model_select,
                       normality = normality)
  class(sum_FH_eblup) <- "summary.fh"
  sum_FH_eblup
}


#' Prints a summary.emdi object
#'
#' The elements described in summary.emdi are printed.
#' @param x an object of type "summary.emdi", generally resulting
#' from applying summary to an object of type "emdi"
#' @param ... optional arguments passed to print.default; see the documentation on
#' that method functions.
#' @export

print.summary.fh <- function(x)
{
  cat("Call:\n ")
  print(x$call)
  cat("\n")
  cat("In-sample domains: ", x$in_sample, "\n")
  cat("Out-of-sample domains: ", x$out_of_sample,
      "\n")
  cat("\n")
  cat("Variance estimation method: ", x$variance_estimation_method,
      "\n")
  cat("Estimated variance of random effects: ", x$sigmau2,
      "\n")
  cat("MSE method: ", x$MSE_method, "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat("Signif. codes: ", x$legend, "\n")
  cat("\n")
  cat("Explanatory measures:\n")
  print(x$model_select)
  cat("\n")
  cat("Residual diagnostics:\n")
  print(x$normality)
  cat("\n")
  if(is.null(x$transformation)){
    cat("Transformation: No transformation \n")
  } else {
    cat("Transformation:\n")
    print(x$transformation)
  }
}
