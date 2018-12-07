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
    lhs(framework$formula) <- quote(helper)
    pred_data <- makeXY(formula = framework$formula, data = pred_data_tmp)
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
    lhs(framework$formula) <- quote(helper)
    pred_data <- makeXY(formula = framework$formula, data = pred_data_tmp)
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

    if (method == "ampl") {
      Bias <- (sum(diag(P - Vi)) + (2/sigmau2)) / sum(diag(Vi^2))
      for (d in 1:framework$m) {
        # Adjusted mse
        mse[d] <- mse[d] - (Bd[d]^2) * Bias
      }
    } else if (method == "amrl") {
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

  if (method == "ampl_yl") {
    Bias <- (sum(diag(P - Vi))) / sum(diag(Vi^2))
    for (d in 1:framework$m) {
      # Adjusted mse
      mse[d] <- mse[d] - (Bd[d]^2) * Bias
    }
  } else if (method == "amrl_yl") {
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


slud_maiti <- function(framework, sigmau2, eblup, combined_data) {

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

    if (method == "sae_reml" | method == "reml") {
      MSE_data <- prasad_rao(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data)
      MSE_method <- "Prasad_Rao"
    } else if (method == "amrl" | method == "ampl") {
      MSE_data <- li_lahiri(framework = framework, sigmau2, combined_data,
                            method = method)
      MSE_method <- paste0(method, "_corrected")
    } else if (method == "amrl_yl" | method == "ampl_yl") {
      MSE_data <- yoshimori_lahiri(framework = framework, sigmau2, combined_data,
                            method = method)
      MSE_method <- paste0(method, "_corrected")
    } else if (method == "ml") {
      MSE_data <- datta_lahiri(framework = framework, sigmau2, combined_data)
      MSE_method <- "Datta_Lahiri"
    }

  mse_out <- list(MSE_data = MSE_data,
                  MSE_method = MSE_method)

  return(mse_out)
  }


boot_arcsin <- function(M, m, sigmau2, vardir, combined_data, framework,
                        eblup, B = 20, method = method,
                        precision = precision, maxiter = maxiter,
                        interval = interval, alpha = alpha) {


  M <- framework$M
  m <- framework$m
  vardir <- framework$vardir
  eff_smpsize <- framework$eff_smpsize
  x <- framework$model_X

  Li <- rep(NA, M)
  Ui <- rep(NA, M)

  ### Bootstrap
    ti <- matrix(NA, M, B)
    boots_est <- matrix(NA, M, B)
    boots_par <- matrix(NA, M, B)
    in_sample <- framework$obs_dom == TRUE
    out_sample <- framework$obs_dom == FALSE

    for (b in 1:B){

      set.seed(b)

      v_boot <- rnorm(M, 0, sqrt(sigmau2))
      e_boot <- rnorm(m, 0, sqrt(vardir))

      # Get covariates for all domains
      pred_data_tmp <- combined_data
      pred_data_tmp <- data.frame(pred_data_tmp, helper = rnorm(1,0,1))
      lhs(framework$formula) <- quote(helper)
      pred_data <- makeXY(formula = framework$formula, data = pred_data_tmp)

      pred_X <- pred_data$x


      Xbeta_boot <- pred_X %*% eblup$coefficients$coefficients


      ## Theta under transformation
      theta <- Xbeta_boot[, 1] + v_boot

      ## Truncation
      true_value_boot <- Xbeta_boot + v_boot
      true_value_boot[true_value_boot < 0] <- 0
      true_value_boot[true_value_boot > (pi / 2)] <- (pi / 2)

      ## Back-transformation
      true_value_boot <- (sin(true_value_boot))^2
      boots_par[,b] <- true_value_boot
      boots_par[,b] <- Xbeta_boot + v_boot

      ystar <- Xbeta_boot[in_sample] + v_boot[in_sample] + e_boot

      ## Estimation of beta_boot
      framework2 <- framework
      framework2$direct <- ystar
      sigmau2_boot <- wrapper_estsigmau2(framework = framework2, method = method,
                                    precision = precision, maxiter = maxiter,
                                    interval = interval)


      ## Computation of the coefficients'estimator (Bstim)
      D <- diag(1, m)
      V <- sigmau2_boot*D%*%t(D) + diag(as.numeric(vardir))
      Vi <- solve(V)
      Q <- solve(t(x)%*%Vi%*%x)
      Beta.hat_boot <- Q%*%t(x)%*%Vi%*%ystar

      ## Computation of the EBLUP
      res <- ystar - c(x%*%Beta.hat_boot)
      Sigma.u <- sigmau2_boot*D
      u.hat <- Sigma.u%*%t(D)%*%Vi%*%res

      ## Small area mean
      est_mean_boot <- x%*%Beta.hat_boot+D%*%u.hat
      Bi <- as.numeric(vardir)/(sigmau2_boot + as.numeric(vardir))
      ti[in_sample, b] <- (theta[in_sample] - est_mean_boot)/sqrt(as.numeric(vardir)*(1 - Bi))

      ## Synthetic prediction for out-of-sample
      pred_out_boot <- pred_X%*%Beta.hat_boot
      ti[out_sample, b]<-(theta[out_sample] - pred_out_boot[out_sample])/sqrt(sigmau2_boot)

      print(b)
    } # End of bootstrap runs

    qi <- matrix(NA, M, 2)

    for (i in 1:M) {
      qi[i,1] <- quantile(ti[i,], prob = alpha/2, na.rm = T)
      qi[i,2] <- quantile(ti[i,], prob = (1 - alpha/2), na.rm = T)
    }

    Li <- matrix(NA, M, 1)
    Ui <- matrix(NA, M, 1)

    Di <- rep(NA, M)
    Di[in_sample] <- vardir
    Di[out_sample] <- (1 / (4 * mean(combined_data[in_sample, eff_smpsize])))
    Bi.tot <- as.numeric(Di) / (sigmau2 + as.numeric(Di))


    Li <- (eblup$EBLUP_data$EBLUP + qi * sqrt(Di * (1 - Bi.tot)))[,1]
    Ui <- (eblup$EBLUP_data$EBLUP + qi * sqrt(Di * (1 - Bi.tot)))[,2]

    ### Truncation
    Li[Li < 0] <- 0
    Ui[Ui > (pi / 2)] <- (pi / 2)

    ### Back-transformation
    Li <- (sin(Li))^2
    Ui <- (sin(Ui))^2

    conf_int <- data.frame(Li = Li, Ui = Ui)
}

jiang_jackknife <- function(framework, combined_data, sigmau2, eblup) {

  m <- framework$M
  jack_sigmau2 <- vector(length = m)
  diff_jack_eblups <- data.frame(row.names = 1:m)
  diff_jack_g1 <- data.frame(row.names = 1:m)

  g1 <- rep(0, framework$m)
  mse <- rep(0, framework$M)
  # Inverse of total variance
  Vi <- 1/(sigmau2 + framework$vardir)
  # Shrinkage factor
  Bd <- framework$vardir/(sigmau2 + framework$vardir)


  for (d in 1:framework$m) {
    # Variance due to random effects: vardir * gamma
    g1[d] <- framework$vardir[d] * (1 - Bd[d])
  }

  for (domain in 1:m) {
    print(domain)

    data_tmp <- combined_data[-domain,]

    # Framework with temporary data
    framework_tmp <- framework_FH(combined_data = data_tmp, fixed = fixed,
                              vardir = vardir, domains = domains,
                              transformation = transformation,
                              eff_smpsize = eff_smpsize)
    # Estimate sigma u
    sigmau2_tmp <- wrapper_estsigmau2(framework = framework_tmp, method = method,
                                  precision = precision, maxiter = maxiter,
                                  interval = interval)
    jack_sigmau2[domain] <- sigmau2_tmp

    Vi_tmp <- 1/(sigmau2_tmp + framework_tmp$vardir)
    # Shrinkage factor
    Bd_tmp <- framework_tmp$vardir/(sigmau2_tmp + framework_tmp$vardir)

    g1_tmp <- rep(0, framework_tmp$m)
    for (d_tmp in 1:framework_tmp$m) {
      g1_tmp[d_tmp] <- framework_tmp$vardir[d_tmp] * (1 - Bd_tmp[d_tmp])
    }

    # G1
    diff_jack_g1[, paste0(domain)] <- g1_tmp - g1

    # Standard EBLUP
    eblup_tmp <- eblup_FH(framework = framework, sigmau2 = sigmau2_tmp,
                      combined_data = combined_data)
    diff_jack_eblups[, paste0(domain)] <- eblup_tmp$EBLUP_data$EBLUP - eblup$EBLUP_data$EBLUP
  }




}
