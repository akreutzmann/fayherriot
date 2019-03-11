
backtransformed <- function(framework, sigmau2, eblup, transformation,
                            combined_data, method, vardir,
                            interval, B, alpha, MSE) {

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
    MSE_method <- estim_MSE$MSE_method

  } else if (transformation == "log_naive") {

    estim_MSE <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                                combined_data = combined_data,
                                method = method)

    EBLUP_data$EBLUP <- exp(eblup$EBLUP_data$EBLUP)
    MSE_data$MSE <- exp(eblup$EBLUP_data$EBLUP)^2 * estim_MSE$MSE_data$MSE
    MSE_method <- estim_MSE$MSE_method

  } else if (transformation == "log_SM") {

    estim_MSE <- analytical_mse(framework = framework, sigmau2 = sigmau2,
                             combined_data = combined_data,
                             method = method)

    int_value <- NULL
    for (i in 1:framework$m) {

      mu_dri <- eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]
      # Get value of first domain
      #mu_dri <- as.numeric(mu_dri[i, 1])
      mu_dri <- mu_dri[i]


      Var_dri <- sigmau2 * (1 - eblup$gamma)
      Var_dri <- as.numeric(Var_dri[i])

      integrand <- function(x, mean, sd){exp(x) * dnorm(x, mean = mu_dri,
                                                        sd = sqrt(Var_dri))}

      upper_bound <- min(mean(framework$direct) + 10 * sd(framework$direct),
                         mu_dri + 100 * sqrt(Var_dri))
      lower_bound <- max(mean(framework$direct) - 10 * sd(framework$direct),
                         mu_dri - 100 * sqrt(Var_dri))


      int_value <- c(int_value, integrate(integrand,
                                          lower = lower_bound,
                                          upper = upper_bound)$value)
    }

    EBLUP_data$EBLUP[framework$obs_dom == TRUE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == TRUE] + (0.5 * sigmau2 * (1 - eblup$gamma)))
    #EBLUP_data$EBLUP[framework$obs_dom == FALSE] <- exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE] + 0.5 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE])
    EBLUP_data$EBLUP[framework$obs_dom == FALSE] <- NA

    SM_MSE <- slud_maiti(framework = framework, sigmau2 = sigmau2,
                         eblup = eblup, combined_data = combined_data)

    MSE_data$MSE[framework$obs_dom == TRUE] <- SM_MSE$MSE[framework$obs_dom == TRUE]
    # MSE_data$MSE[framework$obs_dom == FALSE] <-  exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE])^2 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE]
    MSE_data$MSE[framework$obs_dom == FALSE] <- NA
    MSE_method <- "Slud_Maiti"
    } else if (transformation == "arcsin" & MSE == "boot") {


      EBLUP_data$EBLUP <- eblup$EBLUP_data$EBLUP

      EBLUP_data$EBLUP[EBLUP_data$EBLUP < 0] <- 0
      EBLUP_data$EBLUP[EBLUP_data$EBLUP > (pi / 2)] <- (pi / 2)

      int_value <- NULL
      for (i in 1:framework$m) {

        mu_dri <- eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]
        # Get value of first domain
        mu_dri <- mu_dri[i]


        Var_dri <- sigmau2 * (1 - eblup$gamma)
        Var_dri <- as.numeric(Var_dri[i])

        integrand <- function(x, mean, sd){sin(x)^2 * dnorm(x, mean = mu_dri,
                                                            sd = sqrt(Var_dri))}

        upper_bound <- min(mean(framework$direct) + 10 * sd(framework$direct),
                           mu_dri + 100 * sqrt(Var_dri))
        lower_bound <- max(mean(framework$direct) - 10 * sd(framework$direct),
                           mu_dri - 100 * sqrt(Var_dri))

        int_value <- c(int_value, integrate(integrand,
                                            lower = 0,
                                            upper = pi/2)$value)
      }

      EBLUP_data$EBLUP_corr[eblup$EBLUP_data$ind == 0] <- int_value


      EBLUP_data$EBLUP <- (sin(EBLUP_data$EBLUP))^2
      conf_int <- boot_arcsin(sigmau2 = sigmau2, combined_data = combined_data,
                              framework = framework, eblup = eblup,
                              method = method, interval = interval,
                              B = B, alpha = alpha)
      MSE_boot <- boot_sugasawa2(sigmau2 = sigmau2, combined_data = combined_data,
                                framework = framework, transformation = transformation,
                                eblup = eblup, B = B, MC = 1000, method = method,
                                interval = interval, alpha = alpha)
      MSE_data$MSE <- MSE_boot$MSE
      MSE_data$G1 <- MSE_boot$G1
      MSE_data$G2 <- MSE_boot$G2
      MSE_data$Li <- conf_int$Li
      MSE_data$Ui <- conf_int$Ui
      MSE_method <- "boot"
    } else if (transformation == "arcsin" & MSE != "boot") {

      EBLUP_data$EBLUP <- eblup$EBLUP_data$EBLUP
      EBLUP_data$EBLUP[EBLUP_data$EBLUP < 0] <- 0
      EBLUP_data$EBLUP[EBLUP_data$EBLUP > (pi / 2)] <- (pi / 2)

      int_value <- NULL
      for (i in 1:framework$m) {

        mu_dri <- eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]
        # Get value of first domain
        mu_dri <- mu_dri[i]


        Var_dri <- sigmau2 * (1 - eblup$gamma)
        Var_dri <- as.numeric(Var_dri[i])

        integrand <- function(x, mean, sd){sin(x)^2 * dnorm(x, mean = mu_dri,
                                                          sd = sqrt(Var_dri))}

        upper_bound <- min(mean(framework$direct) + 10 * sd(framework$direct),
                           mu_dri + 100 * sqrt(Var_dri))
        lower_bound <- max(mean(framework$direct) - 10 * sd(framework$direct),
                           mu_dri - 100 * sqrt(Var_dri))

        int_value <- c(int_value, integrate(integrand,
                                            lower = 0,
                                            upper = pi/2)$value)
      }

      EBLUP_data$EBLUP_corr[eblup$EBLUP_data$ind == 0] <- int_value
      EBLUP_data$EBLUP <- (sin(EBLUP_data$EBLUP))^2


      jack_mse <- wrapper_MSE(framework = framework, combined_data = combined_data,
                              sigmau2 = sigmau2, vardir = vardir, eblup = eblup,
                              transformation = transformation, method = method,
                              interval = interval, MSE = MSE)


      back_jack_mse <- 2 * sin(eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]) * cos(eblup$EBLUP_data$EBLUP[eblup$EBLUP_data$ind == 0]) * jack_mse$MSE_data$MSE[jack_mse$MSE_data$ind == 0]

      MSE_data$MSE[framework$obs_dom == TRUE] <- back_jack_mse
      #MSE_data$MSE[framework$obs_dom == FALSE] <-  exp(eblup$EBLUP_data$EBLUP[framework$obs_dom == FALSE])^2 * estim_MSE$MSE_data$MSE[framework$obs_dom == FALSE]
      MSE_data$MSE[framework$obs_dom == FALSE] <- NA
      MSE_method <- jack_mse$MSE_method
    }

  EBLUP_data$ind[framework$obs_dom == TRUE] <- 0
  EBLUP_data$ind[framework$obs_dom == FALSE] <- 1

  MSE_data$ind[framework$obs_dom == TRUE] <- 0
  MSE_data$ind[framework$obs_dom == FALSE] <- 1


  back_out <- list(EBLUP_data = EBLUP_data,
                   MSE_data = MSE_data, MSE_method = MSE_method)

  return(back_out)
  }
