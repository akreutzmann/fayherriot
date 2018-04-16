# Estimation functions for sigmau2


saeReml <- function(vardir, precision, maxiter, X, y) {
  Aest.REML <- 0
  Aest.REML[1] <- median(vardir)
  k <- 0
  diff <- precision + 1
  while ((diff > precision) & (k < maxiter)) {
    k <- k + 1
    Vi <- 1/(Aest.REML[k] + vardir)
    XtVi <- t(Vi * X)
    Q <- solve(XtVi %*% X)
    P <- diag(Vi) - t(XtVi) %*% Q %*% XtVi
    Py <- P %*% y
    s <- (-0.5) * sum(diag(P)) + 0.5 * (t(Py) %*% Py)
    F <- 0.5 * sum(diag(P %*% P))
    Aest.REML[k + 1] <- Aest.REML[k] + s/F
    diff <- abs((Aest.REML[k + 1] - Aest.REML[k])/Aest.REML[k])
  }
  A.REML <- max(Aest.REML[k + 1], 0)
  return(sigmau_reml = A.REML)
}




NicolaReml <- function(interval, direct, x, vardir, areanumber) {

  A.reml <- function(interval, direct, x, vardir, areanumber) {
    psi <- matrix(c(vardir),areanumber,1)
    Y <- matrix(c(direct),areanumber,1)
    X <- x
    Z.area <- diag(1,areanumber)
    sigma.u_log <- interval[1]
    I <- diag(1, areanumber)
    #V is the variance covariance matrix
    V <- sigma.u_log * Z.area%*%t(Z.area) + I*psi[,1]
    Vi <- solve(V)
    Xt <- t(X)
    XVi <- Xt%*%Vi
    Q <- solve(XVi%*%X)
    P <- Vi - (Vi%*%X%*%Q%*%XVi)
    b.s <- Q%*%XVi%*%Y

    ee = eigen(V)
    - (areanumber/2) * log(2*pi)
    - 0.5 * sum(log(ee$value)) - (0.5) * log(det(t(X)%*%Vi%*%X)) - (0.5) * t(Y)%*%P%*%Y
  }
  ottimo <- optimize(A.reml, interval, maximum = TRUE,
                     vardir = vardir, areanumber = areanumber,
                     direct = direct, x = x)

  estsigma2u <- ottimo$maximum

  return(sigmau_reml = estsigma2u)
}


AMRL <- function(interval, direct, x, vardir, areanumber) {

  AR <- function(interval, direct, x, vardir, areanumber){
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

    ee <- eigen(V)
    log(sigma.u_log) - (areanumber/2) * log(2*pi)
    - 0.5 * sum(log(ee$value)) - (0.5) * log(det(t(X)%*%Vi%*%X)) - (0.5) * t(Y)%*%P%*%Y
  }


  ottimo <- optimize(AR, interval, maximum = TRUE,
                     vardir = vardir, areanumber = areanumber,
                     direct = direct, x = x)

  estsigma2u <- ottimo$maximum

  return(sigmau_amrl = estsigma2u)
}


AMPL <- function(interval, direct, x, vardir, areanumber) {

  AP <- function(interval, direct,x, vardir, areanumber){
    psi <- matrix(c(vardir), areanumber, 1)
    Y <- matrix(c(direct), areanumber, 1)
    X <- x
    Z.area <- diag(1, areanumber)
    sigma.u_log <- interval[1]
    I <- diag(1, areanumber)
    # V is the variance covariance matrix
    V <- sigma.u_log * Z.area%*%t(Z.area) + I * psi[,1]
    Vi <- solve(V)
    Xt <- t(X)
    XVi <- Xt%*%Vi
    Q <- solve(XVi%*%X)
    P <- Vi - (Vi%*%X%*%Q%*%XVi)
    b.s <- Q%*%XVi%*%Y

    ee = eigen(V)
    log(sigma.u_log) - (areanumber/2) * log(2*pi)
    - 0.5 * sum(log(ee$value)) - (0.5) * t(Y)%*%P%*%Y
  }

  ottimo <- optimize(AP, interval, maximum = TRUE,
                     vardir = vardir, areanumber = areanumber,
                     direct = direct, x = x)

  estsigma2u <- ottimo$maximum

  return(sigmau_ampl = estsigma2u)
}
