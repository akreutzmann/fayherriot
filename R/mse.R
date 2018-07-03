prasad_rao <- function(framework, sigmau2, combined_data) {

  g1 <- rep(0, framework$m)
  g2 <- rep(0, framework$m)
  g3 <- rep(0, framework$m)
  mse <- rep(0, framework$m)
  # Inverse of total variance
  Vi <- 1/(sigmau2 + framework$vardir)
  # Shrinkage factor
  Bd <- framework$vardir/(sigmau2 + framework$vardir)
  # Squared inverse of total variance
  SumAD2 <- sum(Vi^2)
  # X'Vi
  XtVi <- t(Vi * framework$model_X)
  # (X'ViX)^-1
  Q <- solve(XtVi %*% framework$model_X)

  # 2 divided by squared inverse of total variance
  VarA <- 2/SumAD2

  for (d in 1:framework$m) {
    # Variance due to random effects: vardir * gamma
    g1[d] <- framework$vardir[d] * (1 - Bd[d])
    # Covariate for single domain
    xd <- matrix(framework$model_X[d, ], nrow = 1, ncol = framework$p)
    # Variance due to the estimation of beta
    g2[d] <- (Bd[d]^2) * xd %*% Q %*% t(xd)
    # Variance due to the estimation of the variance of the random effects
    g3[d] <- (Bd[d]^2) * VarA/(sigmau2 + framework$vardir[d])
    # Prasad-Rao estimator
    mse[d] <- g1[d] + g2[d] + 2 * g3[d]
  }

  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir

  # Small area MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0


  if (!all(framework$obs_dom == TRUE)) {
    h <- rep(0, framework$M - framework$m)
    mse_out <- rep(0, framework$M - framework$m)
    # Covariates for out-of-sample domains
    pred_data_tmp <- combined_data[framework$obs_dom == FALSE,]
    pred_data_tmp <- data.frame(pred_data_tmp, helper = rnorm(1,0,1))
    lhs(formula) <- quote(helper)
    pred_data <- makeXY(formula = formula, data = pred_data_tmp)
    pred_X <- pred_data$x

    for (d_out in 1:(framework$M - framework$m)) {
      xd_out <- matrix(pred_X[d_out, ], nrow = 1, ncol = framework$p)
      h[d_out] <- xd_out %*% Q %*% t(xd_out)
      mse_out[d_out] <- sigmau2 + h[d_out]
    }

    MSE_data$MSE[framework$obs_dom == FALSE] <- mse_out
    MSE_data$ind[framework$obs_dom == FALSE] <- 1
  }


  return(MSE_data)
  }


datta_lahiri <- function(framework, sigmau2, combined_data) {


  g1 <- rep(0, framework$m)
  g2 <- rep(0, framework$m)
  g3 <- rep(0, framework$m)
  mse <- rep(0, framework$m)
  Vi <- 1/(sigmau2 + framework$vardir)
  Bd <- framework$vardir/(sigmau2 + framework$vardir)
  SumAD2 <- sum(Vi^2)
  XtVi <- t(Vi * framework$model_X)
  Q <- solve(XtVi %*% framework$model_X)

  VarA <- 2/SumAD2
  b <- (-1) * sum(diag(Q %*% (t((Vi^2) * framework$model_X) %*% framework$model_X)))/SumAD2
  for (d in 1:framework$m) {
    g1[d] <- framework$vardir[d] * (1 - Bd[d])
    xd <- matrix(framework$model_X[d, ], nrow = 1, ncol = framework$p)
    g2[d] <- (Bd[d]^2) * xd %*% Q %*% t(xd)
    g3[d] <- (Bd[d]^2) * VarA/(sigmau2 + framework$vardir[d])
    mse[d] <- g1[d] + g2[d] + 2 * g3[d] - b * (Bd[d]^2)
  }

  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir

  # Small area MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0

  if (!all(framework$obs_dom == TRUE)) {
    h <- rep(0, framework$M - framework$m)
    mse_out <- rep(0, framework$M - framework$m)
    # Covariates for out-of-sample domains
    pred_data_tmp <- combined_data[framework$obs_dom == FALSE,]
    pred_data_tmp <- data.frame(pred_data_tmp, helper = rnorm(1,0,1))
    lhs(formula) <- quote(helper)
    pred_data <- makeXY(formula = formula, data = pred_data_tmp)
    pred_X <- pred_data$x

    for (d_out in 1:(framework$M - framework$m)) {
      xd_out <- matrix(pred_X[d_out, ], nrow = 1, ncol = framework$p)
      h[d_out] <- xd_out %*% Q %*% t(xd_out)
      mse_out[d_out] <- sigmau2 + b + h[d_out]
    }

    MSE_data$MSE[framework$obs_dom == FALSE] <- mse_out
    MSE_data$ind[framework$obs_dom == FALSE] <- 1
  }

  return(MSE_data)
}


li_lahiri <- function(framework, sigmau2, combined_data, method) {


    prasad_rao <- prasad_rao(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data)
    mse <- prasad_rao$MSE[framework$obs_dom == TRUE]
    X <- framework$model_X
    psi <- matrix(c(framework$vardir), framework$m, 1)
    Y <- matrix(c(framework$direct), framework$m, 1)
    Z.area <- diag(1, framework$m)
    I <- diag(1, framework$m)
    # Shrinkage factor
    Bd <- framework$vardir/(sigmau2 + framework$vardir)
    #V is the variance covariance matrix
    V <- sigmau2 * Z.area%*%t(Z.area) + I * psi[,1]
    Vi <- solve(V)
    Xt <- t(X)
    XVi <- Xt%*%Vi
    Q <- solve(XVi%*%X)
    P <- Vi - (Vi%*%X%*%Q%*%XVi)
    b.s <- Q%*%XVi%*%Y

    if (method == "AMPL") {
      Bias <- (sum(diag(P - Vi)) + (2/sigmau2)) / sum(diag(Vi^2))
      for (d in 1:framework$m) {
        # Adjusted mse
        mse[d] <- mse[d] - (Bd[d]^2) * Bias
      }
    } else if (method == "AMRL") {
      Bias <- (2/sigmau2) / sum(diag(Vi^2))
      for (d in 1:framework$m) {
        # Adjusted mse
        mse[d] <- mse[d] - (Bd[d]^2) * Bias
      }
    }

    MSE_data <- prasad_rao
    MSE_data$MSE[framework$obs_dom == TRUE] <- mse
    MSE_data$ind[framework$obs_dom == TRUE] <- 0


    if (!all(framework$obs_dom == TRUE)) {
      MSE_data$MSE[framework$obs_dom == FALSE] <- NA
      MSE_data$ind[framework$obs_dom == FALSE] <- 1

    cat("Please note that only for in-sample-domains a correction following
        Li and Lahiri (2010) is implemented. For the out-of-sample domains,
        no estimate for the MSE is returned. For the reference see help(FH_AK).")
    }

    return(MSE_data)
}


yoshimori_lahiri <- function(framework, sigmau2, combined_data, method) {


  prasad_rao <- prasad_rao(framework = framework, sigmau2 = sigmau2,
                           combined_data = combined_data)
  mse <- prasad_rao$MSE[framework$obs_dom == TRUE]
  X <- framework$model_X
  psi <- matrix(c(framework$vardir), framework$m, 1)
  Y <- matrix(c(framework$direct), framework$m, 1)
  Z.area <- diag(1, framework$m)
  I <- diag(1, framework$m)
  # Shrinkage factor
  Bd <- framework$vardir/(sigmau2 + framework$vardir)
  #V is the variance covariance matrix
  V <- sigmau2 * Z.area%*%t(Z.area) + I * psi[,1]
  Vi <- solve(V)
  Xt <- t(X)
  XVi <- Xt%*%Vi
  Q <- solve(XVi%*%X)
  P <- Vi - (Vi%*%X%*%Q%*%XVi)
  b.s <- Q%*%XVi%*%Y

  if (method == "AMPL_YL") {
    Bias <- (sum(diag(P - Vi))) / sum(diag(Vi^2))
    for (d in 1:framework$m) {
      # Adjusted mse
      mse[d] <- mse[d] - (Bd[d]^2) * Bias
    }
  } else if (method == "AMRL_YL") {
      mse <- mse
  }

  MSE_data <- prasad_rao
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0


  if (!all(framework$obs_dom == TRUE)) {
    MSE_data$MSE[framework$obs_dom == FALSE] <- NA
    MSE_data$ind[framework$obs_dom == FALSE] <- 1

    cat("Please note that only for in-sample-domains a correction following
        Yoshimori and Lahiri (2014) is implemented. For the out-of-sample domains,
        no estimate for the MSE is returned. For the reference see help(FH_AK).")
  }

  return(MSE_data)
}


slud_maiti <- function(framework, sigmau2, eblup) {

  # MSE estimation
  nu <- framework$model_X%*%eblup$coefficients$coefficients
  gamma <- eblup$gamma
  tau <- sigmau2 + framework$vardir
  # Variance of beta
  X <- framework$model_X
  psi <- matrix(c(framework$vardir), framework$m, 1)
  Y <- matrix(c(framework$direct), framework$m, 1)
  Z.area <- diag(1, framework$m)
  I <- diag(1, framework$m)
  #V is the variance covariance matrix
  V <- sigmau2 * Z.area%*%t(Z.area) + I * psi[,1]
  Vi <- solve(V)
  Xt <- t(X)
  XVi <- Xt%*%Vi
  Q <- solve(XVi%*%X)
  Var.beta <- Q
  # Variance of sigmau2
  # Identity matrix mxm
  D <- diag(1, framework$m)
  Deriv1 <- solve((sigmau2 * D) + diag(c(framework$vardir), framework$m))
  ### Inverse of fisher information matrix. That is var. sigma2u
  Var.sigma <- ((1/2) * sum(diag(Deriv1%*%Deriv1)))^(-1)

  tmp <- NULL
  for (i in 1:framework$m) {
    tmp[i] <- (t(framework$model_X)[,i]%*%Q%*%framework$model_X[i,])/(tau[i]^2)
  }

  mse <- NULL
  for (j in 1:framework$m) {
    mse[j] <- exp(2 * (nu[j,1] + sigmau2)) * (1 - exp(-gamma[j] * framework$vardir[j])) +
      ((framework$vardir[j]^2)/tau[j]^2) * exp(2 * nu[j,1] + sigmau2 * (1 + gamma[j])) * (t(framework$model_X)[,j]%*%Q%*%framework$model_X[j,]) +
      Var.sigma * ((framework$vardir[j]^2)/tau[j]^2) * exp(2 * nu[j,1] + sigmau2 * (1 + gamma[j])) *
      ((1/4) * (1 + 3 * gamma[j])^2 + (1/tau[j])) - exp(2 * (nu[j,1] + sigmau2)) *
      (2 * (1 - exp(-gamma[j] * framework$vardir[j])) * (t(framework$model_X)[,j]%*%Q%*%framework$model_X[j,]) - Var.sigma *
         (2 + (((framework$vardir[j]^2)/tau[j]^2) - 2) * exp(-gamma[j] * framework$vardir[j])) * sum(tmp) +
         Var.sigma * (2 + (((2 * (framework$vardir[j]^2))/(tau[j]^2)) - 2) * exp(-gamma[j] * framework$vardir[j]) -
                        ((framework$vardir[j]^2)/(tau[j]^3)) * exp(-gamma[j] * framework$vardir[j]) * (1 + ((framework$vardir[j]^2)/(2 * tau[j])))))
  }


  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0

  if (!all(framework$obs_dom == TRUE)) {
    MSE_data$MSE[framework$obs_dom == FALSE] <- NA
    MSE_data$ind[framework$obs_dom == FALSE] <- 1

    cat("Please note that a MSE is only returned for in-sample domains.
        For more information see help(FH_AK).")
  }

  return(MSE_data)

  }



analytical_mse <- function(framework, sigmau2, combined_data,
                           method) {

    if (method == "sae_reml" | method == "nicola_reml") {
      MSE_data <- prasad_rao(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data)
      MSE_method <- "Prasad_Rao"
    } else if (method == "AMRL" | method == "AMPL") {
      MSE_data <- li_lahiri(framework = framework, sigmau2, combined_data,
                            method = method)
      MSE_method <- paste0(method, "_corrected")
    } else if (method == "AMRL_YL" | method == "AMPL_YL") {
      MSE_data <- yoshimori_lahiri(framework = framework, sigmau2, combined_data,
                            method = method)
      MSE_method <- paste0(method, "_corrected")
    } else if (method == "ML") {
      MSE_data <- datta_lahiri(framework = framework, sigmau2, combined_data)
      MSE_method <- "Datta_Lahiri"
    }

  mse_out <- list(MSE_data = MSE_data,
                  MSE_method = MSE_method)

  return(mse_out)
  }
