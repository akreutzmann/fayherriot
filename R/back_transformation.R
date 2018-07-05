
backtransformed <- function(framework, sigmau2, eblup, transformation,
                            combined_data, method) {

  EBLUP_data <- data.frame(Domain = combined_data[[framework$domains]])
  EBLUP_data$direct <- NA
  EBLUP_data$direct[framework$obs_dom == TRUE] <- framework$direct_orig

  MSE_data <- data.frame(Domain = combined_data[[framework$domains]])
  MSE_data$Var <- NA
  MSE_data$Var[framework$obs_dom == TRUE] <- framework$vardir_orig

  if (transformation == "log_crude") {

    estim_MSE <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                               combined_data = combined_data,
                               method = method)

    EBLUP_data$EBLUP <- exp(eblup$EBLUP_data$EBLUP + 0.5 * estim_MSE$MSE_data$MSE)
    MSE_data$MSE <- exp(eblup$EBLUP_data$EBLUP)^2 * estim_MSE$MSE_data$MSE
  } else if (transformation == "log_SM") {

    estim_MSE <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data,
                             method = method)

    EBLUP_data$EBLUP[framework$obs_dom == TRUE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == TRUE] + (0.5 * sigmau2 * (1 - eblup$gamma)))
    EBLUP_data$EBLUP[framework$obs_dom == FALSE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE] + 0.5 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE])

    SM_MSE <- slud_maiti(framework = framework, sigmau2 = sigmau2,
                         eblup = eblup, combined_data = combined_data)

    MSE_data$MSE[framework$obs_dom == TRUE] <- SM_MSE$MSE[framework$obs_dom == TRUE]
    MSE_data$MSE[framework$obs_dom == FALSE] <-  exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE])^2 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE]
  } else if (transformation == "log_BC2") {

    estim_MSE <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data,
                             method = method)

    SM_MSE <- slud_maiti(framework = framework, sigmau2 = sigmau2,
                         eblup = eblup, combined_data = combined_data)


    D <- diag(1, framework$m)
    Deriv1 <- solve((sigmau2 * D) + diag(c(framework$vardir), framework$m))
    ### Inverse of fisher information matrix. That is var. sigma2u
    II <- ((1/2) * sum(diag(Deriv1%*%Deriv1)))^(-1)

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

    A <- NULL
    for (i in 1:framework$m) {
      A[i] <- ((1 - eblup$gamma[i])^2) * (matrix(framework$model_X[i,],
                                           nrow = 1)%*%Q%*%matrix(framework$model_X[i,], ncol = 1))
    }
    tau <- sigmau2 + framework$vardir
    B1 <- (((1 - eblup$gamma)/tau) * ((framework$direct - c(framework$model_X%*%eblup$coefficients$coefficients))) + (0.5 * ((1 - eblup$gamma)^2)))^2
    B2 <- (((2 - 2 * eblup$gamma)/(tau^2)) * (framework$direct - c(framework$model_X%*%eblup$coefficients$coefficients))) + (((1 - eblup$gamma)^2)/tau)
    c1 <- exp(0.5 * (A + B1 * II))
    c2 <- exp(0.5 * (A + (B1 - B2) * II))

    EBLUP_data$EBLUP[framework$obs_dom == TRUE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == TRUE] + (0.5 * sigmau2 * (1 - eblup$gamma)))/c2
    EBLUP_data$EBLUP[framework$obs_dom == FALSE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE] + 0.5 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE])

    MSE_data$MSE[framework$obs_dom == TRUE] <- SM_MSE$MSE[framework$obs_dom == TRUE] / (c2^2)
    MSE_data$MSE[framework$obs_dom == FALSE] <-  exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE])^2 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE]

  }

  EBLUP_data$ind[framework$obs_dom == TRUE] <- 0
  EBLUP_data$ind[framework$obs_dom == FALSE] <- 1

  MSE_data$ind[framework$obs_dom == TRUE] <- 0
  MSE_data$ind[framework$obs_dom == FALSE] <- 1


  back_out <- list(EBLUP_data = EBLUP_data,
                   MSE_data = MSE_data)

  return(back_out)
  }
