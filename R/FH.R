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
#' @return fitted FH model.
#' @import formula.tools
#' @importFrom stats median model.frame model.matrix model.response optimize
#' @importFrom stats pnorm rnorm
#' @export


FH_AK <- function(formula, vardir, combined_data, domains = NULL, method,
                  back_transformation, interval = c(0, 1000), precision = 0.0001,
                  maxiter = 100) {


  # Get sample and population data
  obs_dom <- !is.na(combined_data[[paste(lhs(formula))]])

  data <- combined_data[obs_dom == TRUE,]

  # Get response variable and model matrix from formula and data
  direct <- makeXY(formula, data)$y
  model_X <- makeXY(formula, data)$x
  vardir <- data[, vardir]



  if (is.null(domains)) {
    data$domains <- 1:length(direct)
    domains <- "domains"
  }

  # Number of areas
  m <- length(direct)
  M <- length(combined_data[[paste(lhs(formula))]])
  # Number of covariates
  p <- ncol(model_X)

  # Estimate sigma u
  if (method == "sae_reml") {
    sigmau2 <- saeReml(vardir = vardir, precision = precision, maxiter = maxiter,
                       X = model_X, y = direct)
  } else if (method == "nicola_reml") {
    sigmau2 <- NicolaReml(interval = interval, vardir = vardir, x = model_X,
                          direct = direct, areanumber = m)
  } else if (method == "AMRL") {
    sigmau2 <- AMRL(interval = interval, vardir = vardir, x = model_X,
                    direct = direct, areanumber = m)
  } else if (method == "AMPL") {
    sigmau2 <- AMPL(interval = interval, vardir = vardir, x = model_X,
                    direct = direct, areanumber = m)
  }

  # Estimation of the regression coefficients
  # Identity matrix mxm
  D <- diag(1, m)
  # Total variance-covariance matrix - only values on the diagonal due to
  # independence of error terms
  V <- sigmau2 * D%*%t(D) + diag(as.numeric(vardir))
  # Inverse of the total variance
  Vi <- solve(V)
  # Inverse of X'ViX
  Q <- solve(t(model_X)%*%Vi%*%model_X)
  # Beta by (X'ViX)^-1 X'Viy
  Beta.hat <- Q%*%t(model_X)%*%Vi%*%direct

  # Inference for coefficients
  std.errorbeta <- sqrt(diag(Q))
  tvalue <- Beta.hat/std.errorbeta
  pvalue <- 2 * pnorm(abs(tvalue), lower.tail = FALSE)

  Coefficients <- data.frame(Coefficients = Beta.hat,
                             Std.Error = std.errorbeta,
                             t.value = tvalue,
                             p.value = pvalue)

  # Computation of the EBLUP
  real_res <- direct - c(model_X%*%Beta.hat)
  sigmau2Diag <- sigmau2*D
  u.hat <- sigmau2Diag%*%t(D)%*%Vi%*%real_res

  # Computation of shrinkage factor
  gamma <- sigmau2 / (sigmau2 + vardir)

  # Small area mean
  if (is.null(back_transformation)) {
    EBLUP <- model_X%*%Beta.hat + D%*%u.hat
  } else if (back_transformation == "naive") {
    EBLUP <- exp(model_X%*%Beta.hat + D%*%u.hat)
  } else if (back_transformation == "SM") {
    EBLUP <- exp(model_X%*%Beta.hat + D%*%u.hat + (0.5 * sigmau2 * (1 - gamma)))
  } else if (back_transformation == "BC2") {
    Deriv1 <- solve((sigmau2 * D) + diag(c(vardir), m))
    ### Inverse of fisher information matrix. That is var. sigma2u
    II <- ((1/2) * sum(diag(Deriv1%*%Deriv1)))^(-1)

    A <- NULL
    for (i in 1:m) {
      A[i] <- ((1 - gamma[i])^2) * (matrix(model_X[i,],
                                            nrow = 1)%*%Q%*%matrix(model_X[i,], ncol = 1))
    }
    tau <- sigmau2 + vardir
    B1 <- (((1 - gamma)/tau) * ((direct - c(model_X%*%Beta.hat))) + (0.5 * ((1 - gamma)^2)))^2
    B2 <- (((2 - 2 * gamma)/(tau^2)) * (direct - c(model_X%*%Beta.hat))) + (((1 - gamma)^2)/tau)
    c1 <- exp(0.5 * (A + B1 * II))
    c2 <- exp(0.5 * (A + (B1 - B2) * II))

    EBLUP <- exp(model_X%*%Beta.hat + D%*%u.hat + (0.5 * sigmau2 * (1 - gamma)))/c2
  }



  # Criteria for model selection
  loglike <- (-0.5) * (sum(log(2 * pi * (sigmau2 + vardir)) +
                             (real_res^2)/(sigmau2 + vardir)))
  AIC <- (-2) * loglike + 2 * (p + 1)
  BIC <- (-2) * loglike + (p + 1) * log(m)

  criteria <- data.frame(loglike = loglike,
                         AIC = AIC,
                         BIC = BIC)

  # Analytical MSE
  # Preparation of matrices to save MSE components
  if (is.null(back_transformation)) {
    g1 <- rep(0, m)
    g2 <- rep(0, m)
    g3 <- rep(0, m)
    mse <- rep(0, m)
    # Inverse of total variance
    Vi <- 1/(sigmau2 + vardir)
    # Shrinkage factor
    Bd <- vardir/(sigmau2 + vardir)
    # Squared inverse of total variance
    SumAD2 <- sum(Vi^2)
    # X'Vi
    XtVi <- t(Vi * model_X)
    # (X'ViX)^-1
    Q <- solve(XtVi %*% model_X)

    # 2 divided by squared inverse of total variance
    VarA <- 2/SumAD2
    for (d in 1:m) {
      # Variance due to random effects: vardir * gamma
      g1[d] <- vardir[d] * (1 - Bd[d])
      # Covariate for single domain
      xd <- matrix(model_X[d, ], nrow = 1, ncol = p)
      # Variance due to the estimation of beta
      g2[d] <- (Bd[d]^2) * xd %*% Q %*% t(xd)
      # Variance due to the estimation of the variance of the random effects
      g3[d] <- (Bd[d]^2) * VarA/(sigmau2 + vardir[d])
      # Prasad-Rao estimator
      mse[d] <- g1[d] + g2[d] + 2 * g3[d]
    }

    if (method == "AMPL" | method == "AMRL") {
      areanumber <- m
      x <- model_X
      psi <- matrix(c(vardir), areanumber, 1)
      Y <- matrix(c(direct), areanumber, 1)
      X <- x
      Z.area <- diag(1, areanumber)
      sigma.u_log <- interval[1]
      I <- diag(1, areanumber)
      #V is the variance covariance matrix
      V <- sigma.u_log * Z.area%*%t(Z.area) + I * psi[,1]
      Vi <- solve(V)
      Xt <- t(X)
      XVi <- Xt%*%Vi
      Q <- solve(XVi%*%X)
      P <- Vi - (Vi%*%X%*%Q%*%XVi)
      b.s <- Q%*%XVi%*%Y

      if (method == "AMPL") {
        Bias <- (sum(diag(P - Vi)) + (2/sigmau2)) / sum(diag(Vi^2))
        for (d in 1:m) {
          # Adjusted mse
          mse[d] <- mse[d] - (Bd[d]^2) * Bias
        }
      } else if (method == "AMRL") {
        Bias <- (2/sigmau2) / sum(diag(Vi^2))
        for (d in 1:m) {
          # Adjusted mse
          mse[d] <- mse[d] - (Bd[d]^2) * Bias
        }
      }
    }

  } else if (back_transformation == "SM" | back_transformation == "BC2") {
    # MSE estimation
    nu <- model_X%*%Beta.hat
    tau <- sigmau2 + vardir
    # Variance of beta
    Var.beta <- Q
    # Variance of sigmau2
    Deriv1 <- solve((sigmau2 * D) + diag(c(vardir), m))
    ### Inverse of fisher information matrix. That is var. sigma2u
    Var.sigma <- ((1/2) * sum(diag(Deriv1%*%Deriv1)))^(-1)

    tmp <- NULL
    for (i in 1:m) {
      tmp[i] <- (t(model_X)[,i]%*%Q%*%model_X[i,])/(tau[i]^2)
    }

    mse <- NULL
    for (j in 1:m) {
      mse[j] <- exp(2 * (nu[j,1] + sigmau2)) * (1 - exp(-gamma[j] * vardir[j])) +
        ((vardir[j]^2)/tau[j]^2) * exp(2 * nu[j,1] + sigmau2 * (1 + gamma[j])) * (t(model_X)[,j]%*%Q%*%model_X[j,]) +
        Var.sigma * ((vardir[j]^2)/tau[j]^2) * exp(2 * nu[j,1] + sigmau2 * (1 + gamma[j])) *
        ((1/4) * (1 + 3 * gamma[j])^2 + (1/tau[j])) - exp(2 * (nu[j,1] + sigmau2)) *
        (2 * (1 - exp(-gamma[j] * vardir[j])) * (t(model_X)[,j]%*%Q%*%model_X[j,]) - Var.sigma *
           (2 + (((vardir[j]^2)/tau[j]^2) - 2) * exp(-gamma[j] * vardir[j])) * sum(tmp) +
           Var.sigma * (2 + (((2 * (vardir[j]^2))/(tau[j]^2)) - 2) * exp(-gamma[j] * vardir[j]) -
           ((vardir[j]^2)/(tau[j]^3)) * exp(-gamma[j] * vardir[j]) * (1 + ((vardir[j]^2)/(2 * tau[j])))))
    }

    if(back_transformation == "BC2") {
      mse <- mse / (c2^2)
    }
  } else if (back_transformation == "naive") {g1 <- rep(0, m)
  g2 <- rep(0, m)
  g3 <- rep(0, m)
  mse <- rep(0, m)
  # Inverse of total variance
  Vi <- 1/(sigmau2 + vardir)
  # Shrinkage factor
  Bd <- vardir/(sigmau2 + vardir)
  # Squared inverse of total variance
  SumAD2 <- sum(Vi^2)
  # X'Vi
  XtVi <- t(Vi * model_X)
  # (X'ViX)^-1
  Q <- solve(XtVi %*% model_X)

  # 2 divided by squared inverse of total variance
  VarA <- 2/SumAD2
  for (d in 1:m) {
    # Variance due to random effects: vardir * gamma
    g1[d] <- vardir[d] * (1 - Bd[d])
    # Covariate for single domain
    xd <- matrix(model_X[d, ], nrow = 1, ncol = p)
    # Variance due to the estimation of beta
    g2[d] <- (Bd[d]^2) * xd %*% Q %*% t(xd)
    # Variance due to the estimation of the variance of the random effects
    g3[d] <- (Bd[d]^2) * VarA/(sigmau2 + vardir[d])
    # Prasad-Rao estimator
    mse[d] <- g1[d] + g2[d] + 2 * g3[d]
  }

    mse <- exp(EBLUP)^2 * mse
  }


  EBLUP_data <- data.frame(Domain = data[[domains]],
                           Direct = direct,
                           EBLUP = as.numeric(EBLUP))

  MSE_data <- data.frame(Domain = data[[domains]],
                         Var = vardir,
                         MSE = mse)

  Gamma <- data.frame(Domain = data[[domains]],
                      Gamma = sigmau2 / (sigmau2 + vardir))

  # Prediction
  if (all(obs_dom == TRUE)) {
    EBLUP_pred <- NULL
  } else {
    pred_data_tmp <- combined_data[obs_dom == FALSE,]

    pred_data_tmp <- data.frame(pred_data_tmp, helper = rnorm(1,0,1))
    lhs(formula) <- quote(helper)
    pred_data <- makeXY(formula = formula, data = pred_data_tmp)

    pred_X <- pred_data$x
    pred_y <- pred_X %*% Beta.hat

    EBLUP_pred <- data.frame(Domain = combined_data[[domains]])
    EBLUP_pred$Pred_FH[obs_dom == TRUE] <- EBLUP
    EBLUP_pred$Pred_FH[obs_dom == FALSE] <- pred_y
    EBLUP_pred$Ind[obs_dom == TRUE] <- 0
    EBLUP_pred$Ind[obs_dom == FALSE] <- 1
  }



  # Return
  out <- list(ind = EBLUP_data,
              MSE = MSE_data,
              ind_pred = EBLUP_pred,
              Coefficients = Coefficients,
              Sigmau2 = sigmau2,
              random_effects = u.hat,
              real_residuals = real_res,
              gamma = Gamma,
              model_select = criteria)

  class(out) <- "FH_AK"

  return(out)

}
