framework_FH <- function(combined_data, formula, vardir, domains,
                         transformation) {
  # Get sample and population data
  obs_dom <- !is.na(combined_data[[paste(lhs(formula))]])

  data <- combined_data[obs_dom == TRUE,]

  # Get response variable and model matrix from formula and data
  direct <- makeXY(formula, data)$y
  model_X <- makeXY(formula, data)$x
  vardir <- data[, vardir]
  direct_orig <- NULL
  vardir_orig <- NULL

  if (transformation == "log_crude" | transformation == "log_SM" | transformation == "log_BC2") {
    direct_orig <- direct
    vardir_orig <- vardir
    vardir <- (1 / direct)^2 * vardir
    direct <- log(direct)
  }



  if (is.null(domains)) {
    data$domains <- 1:length(direct)
    domains <- "domains"
  }

  # Number of areas
  m <- length(direct)
  M <- length(combined_data[[paste(lhs(formula))]])
  # Number of covariates
  p <- ncol(model_X)

  framework_out <- list(obs_dom = obs_dom,
                        data = data,
                        direct = direct,
                        direct_orig = direct_orig,
                        model_X = model_X,
                        vardir = vardir,
                        vardir_orig = vardir_orig,
                        domains = domains,
                        m = m,
                        M = M,
                        p = p)

  return(framework_out)
  }