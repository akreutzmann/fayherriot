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
      MSE_method <- "Li_Lahiri"
    } else if (method == "ampl_yl") {
      MSE_data <- yoshimori_lahiri(framework = framework, sigmau2, combined_data,
                            method = method)
      MSE_method <- "Yoshimori_Lahiri"
    } else if (method == "amrl_yl") {
      MSE_data <- yoshimori_lahiri(framework = framework, sigmau2, combined_data,
                                   method = method)
      MSE_method <- "Prasad_Rao"
    } else if (method == "ml") {
      MSE_data <- datta_lahiri(framework = framework, sigmau2, combined_data)
      MSE_method <- "Datta_Lahiri"
    }

  mse_out <- list(MSE_data = MSE_data,
                  MSE_method = MSE_method)

  return(mse_out)
  }


boot_arcsin <- function(sigmau2, vardir, combined_data, framework,
                        eblup, B, method,
                        interval, alpha) {


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

      #set.seed(b)

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

jiang_jackknife <- function(framework, combined_data, sigmau2, eblup, transformation,
                            vardir, method, interval) {


  # this MSE estimator can leed to negative values
  m <- framework$m
  jack_sigmau2 <- vector(length = m)
  diff_jack_eblups <- data.frame(row.names = 1:m)
  diff_jack_g1 <- data.frame(row.names = 1:m)

  g1 <- rep(0, framework$m)
  jack_mse <- rep(0, framework$m)
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


    data_insample <- combined_data[framework$obs_dom,]
    data_tmp <- data_insample[-domain,]

    # Framework with temporary data
    framework_tmp <- framework_FH(combined_data = data_tmp, fixed = framework$formula,
                              vardir = vardir, domains = framework$domains,
                              transformation = transformation,
                              eff_smpsize = framework$eff_smpsize)
    # Estimate sigma u
    sigmau2_tmp <- wrapper_estsigmau2(framework = framework_tmp, method = method,
                                  interval = interval)
    jack_sigmau2[domain] <- sigmau2_tmp

    Vi_tmp <- 1/(sigmau2_tmp + framework$vardir)
    # Shrinkage factor
    Bd_tmp <- framework$vardir/(sigmau2_tmp + framework$vardir)

    g1_tmp <- rep(0, framework$m)
    for (d_tmp in 1:framework$m) {
      g1_tmp[d_tmp] <- framework$vardir[d_tmp] * (1 - Bd_tmp[d_tmp])
    }

    # G1
    diff_jack_g1[, paste0(domain)] <- g1_tmp - g1

    # Standard EBLUP
    framework_insample <- framework_FH(combined_data = data_insample, fixed = framework$formula,
                                       vardir = vardir, domains = framework$domains,
                                       transformation = transformation,
                                       eff_smpsize = framework$eff_smpsize)
    eblup_tmp <- eblup_FH(framework = framework_insample, sigmau2 = sigmau2_tmp,
                      combined_data = data_insample)
    diff_jack_eblups[, paste0(domain)] <- eblup_tmp$EBLUP_data$EBLUP - eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]
  }

  jack_mse <- g1 - ((m - 1)/m) * rowSums(diff_jack_g1) + ((m - 1)/m) * rowSums(diff_jack_eblups^2)


  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir

  # Jackknife MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- jack_mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0


  if (!all(framework$obs_dom == TRUE)) {
    MSE_data$MSE[framework$obs_dom == FALSE] <- NA
    MSE_data$ind[framework$obs_dom == FALSE] <- 1

    cat("Please note that the jackknife MSE is only available for in-sample
        domains.")
  }

  mse_out <- list(MSE_data = MSE_data,
                  MSE_method = "jackknife")

  return(mse_out)
}

chen_weighted_jackknife <- function(framework, combined_data, sigmau2, eblup, transformation,
                                    vardir, method, interval) {
  
  # implementiert nach:
  # Chen S., Lahiri P. (2002), A Weighted Jackknife MSPE Estimator in Small-Area Estimation,
  # "Proceeding of the Section on Survey Research Methods", American Statistical Association,
  # pp. 473-477.
  
  # hier wurden erstmal die weights aus d-1/d festgelegt, es gibt auch andere Optionen 
  
  m <- framework$m
  jack_sigmau2 <- vector(length = m)
  diff_jack_eblups <- data.frame(row.names = 1:m)
  diff_jack_g1 <- data.frame(row.names = 1:m)
  diff_jack_g2 <- data.frame(row.names = 1:m)
  
  g1 <- rep(0, framework$m)
  jack_mse <- rep(0, framework$m)
  jack_mse_weighted <- rep(0, framework$m)
  # Inverse of total variance
  Vi <- 1/(sigmau2 + framework$vardir)
  # Shrinkage factor
  Bd <- framework$vardir/(sigmau2 + framework$vardir)
  
  
  for (d in 1:framework$m) {
    # Variance due to random effects: vardir * gamma
    g1[d] <- framework$vardir[d] * (1 - Bd[d])
  }
  
  
  Nenner <- rep(0, framework$m)
  for (d in 1:framework$m) {
    Nenner[d] <- (framework$model_X[d,] %*% framework$model_X[d,]) / 
      (framework$vardir[d] + sigmau2)
  }
  
  
  g2 <- rep(0, framework$m)
  for (d in 1:framework$m) {
    # Variance due to beta estimation
    g2[d] <- (Bd[d])^2 * framework$model_X[d,] %*% framework$model_X[d,] * (sum(Nenner))^(-1)
  }
  
  for (domain in 1:m) {
    print(domain)
    
    
    data_insample <- combined_data[framework$obs_dom,]
    data_tmp <- data_insample[-domain,]
    
    # Framework with temporary data
    framework_tmp <- framework_FH(combined_data = data_tmp, fixed = framework$formula,
                                  vardir = vardir, domains = framework$domains,
                                  transformation = transformation,
                                  eff_smpsize = framework$eff_smpsize)
    # Estimate sigma u
    sigmau2_tmp <- wrapper_estsigmau2(framework = framework_tmp, method = method,
                                      interval = interval)
    jack_sigmau2[domain] <- sigmau2_tmp
    
    Vi_tmp <- 1/(sigmau2_tmp + framework$vardir)
    # Shrinkage factor
    Bd_tmp <- framework$vardir/(sigmau2_tmp + framework$vardir)
    
    g1_tmp <- rep(0, framework$m)
    for (d_tmp in 1:framework$m) {
      g1_tmp[d_tmp] <- framework$vardir[d_tmp] * (1 - Bd_tmp[d_tmp])
    }
    
    Nenner_tmp <- rep(0, framework$m)
    for (d_tmp in 1:framework$m) {
      Nenner_tmp[d_tmp] <- (framework$model_X[d_tmp,] %*% framework$model_X[d_tmp,]) / 
        (framework$vardir[d_tmp] + sigmau2_tmp)
    }
    
    g2_tmp <- rep(0, framework$m)
    for (d_tmp in 1:framework$m) {
      g2_tmp[d_tmp] <- (Bd_tmp[d_tmp])^2 * framework$model_X[d_tmp,] %*% framework$model_X[d_tmp,] * (sum(Nenner_tmp))^(-1)
    }
    
    #G1
    diff_jack_g1[, paste0(domain)] <- g1_tmp - g1
    # negative Werte koennen erhalten werden, wenn diff_jack_g1 sehr groß ist
    
    #G1
    diff_jack_g2[, paste0(domain)] <- g1_tmp + g2_tmp - (g1 + g2)
    
    # Standard EBLUP
    framework_insample <- framework_FH(combined_data = data_insample, fixed = framework$formula,
                                       vardir = vardir, domains = framework$domains,
                                       transformation = transformation,
                                       eff_smpsize = framework$eff_smpsize)
    eblup_tmp <- eblup_FH(framework = framework_insample, sigmau2 = sigmau2_tmp,
                          combined_data = data_insample)
    diff_jack_eblups[, paste0(domain)] <- eblup_tmp$EBLUP_data$EBLUP - eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]
  }
  
  w_u  <- c()
  v_wj <- c()
  
  for(i in 1:nrow(framework$model_X)){
    w_u[i] <- 1 - t(framework$model_X[i,]) %*% solve(t(framework$model_X) %*% framework$model_X) %*% framework$model_X[i,]
    v_wj[i] <- w_u[i] * (jack_sigmau2[i] - sigmau2) * (jack_sigmau2[i] - sigmau2)
  }
  
  jack_mse <- g1 - ((m - 1)/m) * rowSums(diff_jack_g1) + ((m - 1)/m) * rowSums(diff_jack_eblups^2)
  jack_mse_weighted <- g1 + g2 - ((m - 1)/m) * rowSums(diff_jack_g2) + ((m - 1)/m) * rowSums(diff_jack_eblups^2)
  
  jack_mse <- g1 - w_u * rowSums(diff_jack_g1) + w_u * rowSums(diff_jack_eblups^2)
  jack_mse_weighted <- g1 + g2 - w_u * rowSums(diff_jack_g2) + w_u * rowSums(diff_jack_eblups^2)
  
  
  # bias correction for negative values:
  
  neg_values <- which(jack_mse_weighted < 0)
  
  jack_mse_weighted_neg <- c()
  
  if (length(neg_values) >0 ) {
    
    Sig_d<- (sigmau2 + framework$vardir)* diag(m)
    Sig_d_inv <- solve(Sig_d)
    L_d  <- framework$vardir/((sigmau2 + framework$vardir)^2) * diag(m)
    b_wj <- w_u * (jack_sigmau2 - sigmau2)
    
    # bias fuer REML ist NULL (nur fuer REML durchfuehrbar)
    app_bias_correction <- 
      sum(b_wj) * (framework$vardir^2 / (framework$vardir + sigmau2)^2) - #ueberprueft
      diag(L_d %*% Sig_d %*% t(L_d) * sum(v_wj))
    
    test <- diag(L_d %*% (framework$direct - framework$model_X %*% Beta.hat) %*% 
                   t(framework$direct - framework$model_X %*% Beta.hat) %*% t(L_d) * sum(v_wj))
    
    jack_mse_weighted_neg <- g1 + g2 + w_u * rowSums(diff_jack_eblups^2) - app_bias_correction 
  }
  
  for(i in 1:m){
    if(i %in% neg_values){
      jack_mse_weighted[i] <- jack_mse_weighted_neg[i]
    }
  }
  
  jack_mse_weighted 
  
  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir
  
  # Jackknife MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- jack_mse_weighted
  MSE_data$ind[framework$obs_dom == TRUE] <- 0
  
  
  if (!all(framework$obs_dom == TRUE)) {
    MSE_data$MSE[framework$obs_dom == FALSE] <- NA
    MSE_data$ind[framework$obs_dom == FALSE] <- 1
    
    cat("Please note that the jackknife MSE is only available for in-sample
        domains.")
  }
  
  mse_out <- list(MSE_data = MSE_data,
                  MSE_method = "jackknife weighted")
  
  return(mse_out)
}

boot_sugasawa <- function(sigmau2, vardir, combined_data, framework,
                        eblup, B, method, transformation = transformation,
                        interval, alpha) {


  M <- framework$M
  m <- framework$m
  vardir <- framework$vardir
  eff_smpsize <- framework$eff_smpsize
  x <- framework$model_X

  mse <- rep(NA, m)

  ### Bootstrap
  g1_star <- matrix(NA, m, B)
  g2_star <- matrix(NA, m, B)
  in_sample <- framework$obs_dom == TRUE
  out_sample <- framework$obs_dom == FALSE

  for (b in 1:B){

    #set.seed(b)
    # Generate bootstrap samples on transformed scale

    v_boot <- rnorm(m, 0, sqrt(sigmau2))
    e_boot <- rnorm(m, 0, sqrt(vardir))


    # Bootstrap sample on transformed scale
    h_star <- x %*% eblup$coefficients$coefficients + v_boot + e_boot

    # Truncation
    h_star[h_star < 0] <- 0
    h_star[h_star > (pi / 2)] <- (pi / 2)

    y_star <- (sin(h_star))^2

    combined_data_boot <- combined_data
    combined_data_boot[in_sample, paste(framework$formula[2])] <- y_star


    framework_boot <- framework_FH(combined_data = combined_data_boot,
                                   fixed = framework$formula,
                                   vardir = framework$vardir_orig, domains = framework$domains,
                                   transformation = transformation,
                                   eff_smpsize = framework$eff_smpsize)

    # Estimate sigma u -----------------------------------------------------------
    sigmau2_boot <- wrapper_estsigmau2(framework = framework_boot, method = method,
                                  interval = interval)
    # sigmau2 star is estimated

    # Standard EBLUP -------------------------------------------------------------
    eblup_boot <- eblup_FH(framework = framework_boot, sigmau2 = sigmau2_boot,
                      combined_data = combined_data_boot)
    # beta star is estimated


    # Estimation of g1 on transformed scale for bootstrap estimates
    g1_boot <- rep(0, framework$m)
    # Inverse of total variance
    Vi_boot <- 1/(sigmau2_boot + framework$vardir)
    # Shrinkage factor
    Bd_boot <- framework$vardir/(sigmau2_boot + framework$vardir)
    # Squared inverse of total variance
    #SumAD2_boot <- sum(Vi_boot^2)
    # X'Vi
    #XtVi_boot <- t(Vi_boot * framework$model_X)
    # (X'ViX)^-1
    #Q_boot <- solve(XtVi_boot %*% framework$model_X)

    # 2 divided by squared inverse of total variance
    #VarA_boot <- 2/SumAD2_boot

    for (d in 1:framework$m) {
      # Variance due to random effects: vardir * gamma
      g1_boot[d] <- framework$vardir[d] * (1 - Bd_boot[d])
      # Covariate for single domain
      #xd_boot <- matrix(framework$model_X[d, ], nrow = 1, ncol = framework$p)
      # Variance due to the estimation of beta
    }

    # Truncation
    #mu_star_temp <- x %*% eblup$coefficient$coefficients + v_boot

    #mu_star_temp[mu_star_temp < 0] <- 0
    #mu_star_temp[mu_star_temp > (pi / 2)] <- (pi / 2)

    #mu_star <- (sin(mu_star_temp))^2

    mu_star_int <- NULL
    for (i in 1:framework$m) {
      mu_dri_int <- x %*% eblup$coefficient$coefficients
      # Get value of first domain
      mu_dri_int <- mu_dri_int[i]


      Var_dri_int <- sigmau2
      #Var_dri_boot <- as.numeric(Var_dri_boot[i])

      integrand_int <- function(x, mean, sd){sin(x)^2 * dnorm(x, mean = mu_dri_int,
                                                          sd = sqrt(Var_dri_int))}
      mu_star_int <- c(mu_star_int, integrate(integrand_int,
                                          lower = 0,
                                          upper = pi/2)$value)
    }


    # Get the back_transformed estimate
    int_value <- NULL
    for (i in 1:framework$m) {
      mu_dri <- eblup$EBLUP_data$EBLUP[in_sample]
      # Get value of first domain
      mu_dri <- mu_dri[i]


      Var_dri <- sigmau2 * (1 - eblup$gamma)
      Var_dri <- as.numeric(Var_dri[i])

      integrand <- function(x, mean, sd){sin(x)^2 * dnorm(x, mean = mu_dri,
                                                          sd = sqrt(Var_dri))}
      int_value <- c(int_value, integrate(integrand,
                                          lower = 0,
                                          upper = pi/2)$value)
    }


    g2_boot <- (mu_star_int - int_value)^2

    g1_star[, b] <- g1_boot
    g2_star[, b] <- g2_boot
  }

  g1 <- rep(0, framework$m)
  # Inverse of total variance
  Vi <- 1/(sigmau2 + framework$vardir)
  # Shrinkage factor
  Bd <- framework$vardir/(sigmau2 + framework$vardir)

  for (d in 1:framework$m) {
    # Variance due to random effects: vardir * gamma
    g1[d] <- framework$vardir[d] * (1 - Bd[d])
  }

  g1_bc <- 2 * g1 - as.numeric(rowMeans(g1_star))
  g2_star_exp <- as.numeric(rowMeans(g2_star))

  mse <- g1_bc + g2_star_exp

  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir

  # Small area MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0
  MSE_data$MSE[framework$obs_dom == FALSE] <- NA
  MSE_data$ind[framework$obs_dom == FALSE] <- 1

  MSE_data$G1[framework$obs_dom == TRUE] <- g1_bc
  MSE_data$G1[framework$obs_dom == FALSE] <- NA

  MSE_data$G2[framework$obs_dom == TRUE] <- g2_star_exp
  MSE_data$G2[framework$obs_dom == FALSE] <- NA

  return(MSE_data)
}


boot_sugasawa2 <- function(sigmau2, vardir, combined_data, framework,
                          eblup, B, MC, method, transformation = transformation,
                          interval, alpha) {

  M <- framework$M
  m <- framework$m
  vardir <- framework$vardir
  direct <- framework$direct
  eff_smpsize <- framework$eff_smpsize
  X <- framework$model_X
  in_sample <- framework$obs_dom == TRUE
  out_sample <- framework$obs_dom == FALSE

  # Arcsin transformation
  ArcsinTrafo = function(y){
    y[y < 0] <- 0
    y[y > 1] <- 1
    asin(sqrt(y))
  }

  # Inverse of arcsin transformation
  invArcsin = function(y){
    y[y < 0] <- 0
    y[y > pi/2] <- pi/2
    sin(y)^2
  }


  generate.PTFH = function(m, sigmau2, vardir, eblup){
    # Number of domains
    # m = dim(X)[1]
    # Number of covariates
    # p = dim(X)[2]
    # Estimated lambda
    # la = para[1]
    # Estimated sigmau2
    # A = para[2]
    # A = sigmau2
    # Estimated coefficients
    #beta = para[3:(2+p)]
    beta = eblup$coefficients$coefficients

    # "True" values: synthetic part plus N(0, A)
    true = X%*%beta + rnorm(m, 0, sqrt(sigmau2))
    true = as.vector(true)
    # "Observed" values: "true" values plus N(0, D)
    obs = true + rnorm(m, 0, sqrt(vardir))
    # Before values are returned, these are back-transformed by the
    # inverse transfsormation
    return(list(invArcsin(true), invArcsin(obs)))
  }

  g1 = function(m, X, vardir, sigmau2, eblup, MC) {
    # Estimated lambda
    #la=para[1]
    # Estimated sigmau2
    # A=para[2]
    # Estimated coefficients
    beta <- eblup$coefficients$coefficients
    # Synthetic parts
    mu <- as.vector(X%*%beta)

    # G1 estimation following equation 2.5
    a = sigmau2 / (sigmau2 + vardir)
    c1 = sqrt((1 + a) / 2)
    c2 = sqrt((1 - a) / 2)

    int = c()
    for(i in 1:m) {
      # MC specifies the number of MC iterations
      # Large samples of z1 and z2 are
      z1 = rnorm(MC, 0, sqrt(sigmau2))
      z2 = rnorm(MC, 0, sqrt(sigmau2))
      rn = invArcsin(mu[i] + z1)^2 - invArcsin(mu[i] + c1[i] * z1 + c2[i] * z2) *
        invArcsin(mu[i] + c1[i] * z1 - c2[i] * z2)
      int[i] = mean(rn)
    }
    return(int)
  }


  predict = function(direct, X, m, sigmau2, vardir, eblup, MC){
    # Estimated lambda
    # la = para[1]
    # Estimated sigmau2
    # A = para[2]
    # Estimated coefficients
    beta = eblup$coefficients$coefficients
    # FH estimator: (1 - gomma) * synthetic part + gamma * h(direct)
    th = vardir / (sigmau2 + vardir) * (X%*%beta) + (sigmau2 / (sigmau2 + vardir)) * ArcsinTrafo(direct)
    # Variance of ??
    si = sigmau2 * vardir / (sigmau2 + vardir)

    # Prediction of ???
    pred = c()
    for(i in 1:m){
      rn = rnorm(MC, th[i], sqrt(si)[i])
      pred[i] = mean(invArcsin(rn))
    }
    return(pred)
  }

  # Vectors for g1 and g2 part
  g1b = matrix(NA, B, m)
  g2b = matrix(NA, B, m)

  # Bootstrap
  for(b in 1:B) {
    # Get "true" and "observed" values
    dd = generate.PTFH(m = m, sigmau2 = sigmau2, vardir = vardir, eblup = eblup)
    # Back-transformed true values
    true = dd[[1]]
    # Back-transformed observed values
    boot = dd[[2]]
    # Estimate parameters with back-transformed observed values
    combined_data_boot <- combined_data
    combined_data_boot[in_sample, paste(framework$formula[2])] <- boot

    framework_boot <- framework_FH(combined_data = combined_data_boot,
                                   fixed = framework$formula,
                                   vardir = vardir, domains = framework$domains,
                                   transformation = transformation,
                                   eff_smpsize = eff_smpsize)


    # Estimate sigma u -----------------------------------------------------------
    sigmau2_boot <- wrapper_estsigmau2(framework = framework_boot, method = method,
                                  interval = interval)


    # Standard EBLUP -------------------------------------------------------------
    eblup_boot <- eblup_FH(framework = framework_boot, sigmau2 = sigmau2_boot,
                      combined_data = combined_data_boot)

    # Predict ??? with "observed" bootstrap values and estimated parameters
    # obtained from "observed" bootstrap values
    EBE <- predict(direct = framework_boot$direct, X = X, m = m, sigmau2 = sigmau2_boot,
                   vardir = vardir, eblup = eblup_boot, MC = MC)
    # Predict ??? with "observed" bootstrap values and estimated parameters
    # obtained from original data
    BE = predict(direct = framework_boot$direct, X = X, m = m, sigmau2 = sigmau2,
                 vardir = vardir, eblup = eblup, MC = MC)
    # Calculate g1 following equation 2.5
    g1b[b, ] = g1(m = m, X = X, vardir = vardir, sigmau2 = sigmau2_boot, eblup = eblup_boot,
                  MC = MC)
    # Calculate g2 following step 3
    g2b[b, ] = (EBE - BE)^2
  }
  
  

  # Calculate MSE following equation 2.7
  mse = 2 * g1(m = m, X = X, vardir = vardir, sigmau2 = sigmau2, eblup = eblup,
               MC = MC) - apply(g1b, 2, mean) + apply(g2b, 2, mean)

  g1_bc <- 2 * g1(m = m, X = X, vardir = vardir, sigmau2 = sigmau2, eblup = eblup,
                  MC = MC) - apply(g1b, 2, mean)

  g2_boot <- apply(g2b, 2, mean)



  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir

  # Small area MSE
  MSE_data$MSE[framework$obs_dom == TRUE] <- mse
  MSE_data$ind[framework$obs_dom == TRUE] <- 0
  MSE_data$MSE[framework$obs_dom == FALSE] <- NA
  MSE_data$ind[framework$obs_dom == FALSE] <- 1

  MSE_data$G1[framework$obs_dom == TRUE] <- g1_bc
  MSE_data$G1[framework$obs_dom == FALSE] <- NA

  MSE_data$G2[framework$obs_dom == TRUE] <- g2_boot
  MSE_data$G2[framework$obs_dom == FALSE] <- NA

  return(MSE_data)
}



wrapper_MSE <- function(framework, combined_data, sigmau2, vardir, eblup,
                        transformation, method, interval, MSE) {
  MSE_data <- if (MSE == "analytical") {
    analytical_mse(framework = framework, sigmau2 = sigmau2,
                   combined_data = combined_data, method = method)
  } else if (MSE == "jackknife") {
    jiang_jackknife(framework = framework, combined_data = combined_data,
                    sigmau2 = sigmau2, vardir = vardir, eblup = eblup,
                    transformation = transformation, method = method,
                    interval = interval)
  } else if (MSE == "jackknife_w"){
    chen_weighted_jackknife(framework = framework, combined_data = combined_data,
                            sigmau2 = sigmau2, eblup = eblup, vardir = vardir,
                            transformation = transformation, method = method, 
                            interval = interval)
  }
  
  return(MSE_data)
}
