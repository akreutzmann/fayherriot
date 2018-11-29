#===================================================================================================
# SUBJECT:  Fay-Herriot function with arcsin transformation
# AUTHOR:   Timo Schmid, Fabian Bruckschen, Nicola Salvati and Till Zbiranski, June 06, 2017
# TITLE:    Example code for the paper "Constructing socio-demographic indicators for National
#           Statistical Institutes using mobile phone data: estimating literacy rates in Senegal"
#===================================================================================================

# Description of Input -----------------------------------------------------------------------------

# formula = formula object to specify the relationship between y and x
# vardir = approx. the variance of the direct estimator
# dataframe_sample = dataframe with the covariates x and y from the sample
# dataframe_pop_aux = dataframe with the covariates x from the population
# saind = id for the in-sample areas
# x.total_saind = id for the areas in the population
# area_count = Effective sample size
# Boot T/F = Running a bootstrap for constructing confidence intervals
# B = Number of bootstrap replications
# method = Fitting method for the Fay-Herriot; options =
#                 "REML": Restricted ML
#                 "AP": Adjusted maximum profile likelihood (Li and Lahiri 2010)
#                 "AR": Adjusted maximum residual likelihood (Li and Lahiri 2010)
# alpha = Quantile for constructing (1-alpha) confidence intervals
# bench_value = Benchmarking value for internal consistency, e.g. the national mean
# ratio = Weights for the areas for benchmarking, e.g. the population shares per area

# Defining Main Functions --------------------------------------------------------------------------

FH_arcsin <- function(formula,vardir,dataframe_sample, saind,dataframe_pop_aux,
                      x.total_saind,area_count,Boot=FALSE,B=100,method=c("REML","AP","AR"),
                      alpha=0.05,bench_value,ratio)
{

  library(MASS)
  library(nlme)
  library(formula.tools)

  benchmark_mean <- function(thatb, T, Phi, w, gamma, Omega) {
    sig <- Phi
    siginv <- diag(1/diag(sig))
    delta.hat <- siginv%*%(Phi %*% thatb + t(as.numeric(solve(t(w) %*% siginv %*% w))*
                                               (T - t(w) %*% siginv %*% Phi %*% thatb)%*%w))
    return(delta.hat)
  }

  ### REML log-likelihood function
  logl=function(delta,vardir_logl,areanumber_logl,direct_logl,x_logl){
    psi=matrix(c(vardir_logl),areanumber_logl,1)
    Y=matrix(c(direct_logl),areanumber_logl,1)
    X=x_logl
    Z.area=diag(1,areanumber_logl)
    sigma.u_log<-delta[1]
    I<-diag(1,areanumber_logl)
    #V is the variance covariance matrix
    V<-sigma.u_log*Z.area%*%t(Z.area)+I*psi[,1]
    Vi<-solve(V)
    Xt=t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    b.s<-Q%*%XVi%*%Y

    ee=eigen(V)
    -(areanumber_logl/2)*log(2*pi)
    -0.5*sum(log(ee$value))-(0.5)*log(det(t(X)%*%Vi%*%X))-(0.5)*t(Y)%*%P%*%Y
  }

  ### AR log-likelihood function
  ARlogl=function(delta,vardir_logl,areanumber_logl,direct_logl,x_logl){
    psi=matrix(c(vardir_logl),areanumber_logl,1)
    Y=matrix(c(direct_logl),areanumber_logl,1)
    X=x_logl
    Z.area=diag(1,areanumber_logl)
    sigma.u_log<-delta[1]
    I<-diag(1,areanumber_logl)
    #V is the variance covariance matrix
    V<-sigma.u_log*Z.area%*%t(Z.area)+I*psi[,1]
    Vi<-solve(V)
    Xt=t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    b.s<-Q%*%XVi%*%Y

    ee=eigen(V)
    log(sigma.u_log)-(areanumber_logl/2)*log(2*pi)
    -0.5*sum(log(ee$value))-(0.5)*log(det(t(X)%*%Vi%*%X))-(0.5)*t(Y)%*%P%*%Y
  }

  ### AP log-likelihood function
  APlogl=function(delta,vardir_logl,areanumber_logl,direct_logl,x_logl){
    psi=matrix(c(vardir_logl),areanumber_logl,1)
    Y=matrix(c(direct_logl),areanumber_logl,1)
    X=x_logl
    Z.area=diag(1,areanumber_logl)
    sigma.u_log<-delta[1]
    I<-diag(1,areanumber_logl)
    #V is the variance covariance matrix
    V<-sigma.u_log*Z.area%*%t(Z.area)+I*psi[,1]
    Vi<-solve(V)
    Xt=t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    b.s<-Q%*%XVi%*%Y

    ee=eigen(V)
    log(sigma.u_log)-(areanumber_logl/2)*log(2*pi)-0.5*sum(log(ee$value))-(0.5)*t(Y)%*%P%*%Y
  }

# Defining Main Vectors ----------------------------------------------------------------------------

  ni_sampled<-(table(saind))
  n<-length(saind)
  areanumber_sample<-length(ni_sampled)
  areanumber<-length(unique(x.total_saind))
  ni_ind<-rep(0,areanumber)

  for(i in 1:areanumber){
    ni_ind[i]<- max((as.numeric((attr(table(x.total_saind),"dimnames"))$x.total_saind)[i]==
                       as.numeric(attr(ni_sampled,"dimnames")$saind)))
  }

  ni<-rep(0,areanumber)
  ni[as.logical(ni_ind)]<-ni_sampled

  makeXY <- function(formula, data){
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    list(y = y,
         x = x)
  }

  daten<-makeXY(formula=formula,data=dataframe_sample)
  y<-daten$y
  x<-daten$x
  p<-ncol(x)

# Arcsin-Transformation ----------------------------------------------------------------------------

  y<-asin(sqrt(y))
  vardir <-  1/(4*area_count)

# Model Fitting ------------------------------------------------------------------------------------

  if (method=="REML")ottimo=suppressWarnings( optimize(logl,c(-5,100),maximum = TRUE,
                vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=y,x_logl=x))
  if (method=="AR")ottimo=suppressWarnings( optimize(ARlogl,c(-5,100),maximum = TRUE,
                vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=y,x_logl=x))
  if (method=="AP")ottimo=suppressWarnings( optimize(APlogl,c(-5,100),maximum = TRUE,
                vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=y,x_logl=x))

# Model Estimation ---------------------------------------------------------------------------------

  estsigma2u<-ottimo$maximum

  ### Computation of the coefficients'estimator (Bstim)
  D=diag(1,areanumber_sample)
  V<-estsigma2u*D%*%t(D)+diag(as.numeric(vardir))
  Vi<-solve(V)
  Q<-solve(t(x)%*%Vi%*%x)
  Beta.hat<-Q%*%t(x)%*%Vi%*%y

  ### Computation of the EBLUP
  res<-y-c(x%*%Beta.hat)
  I<-diag(1,areanumber_sample)
  Sigma.u=estsigma2u*I
  u.hat=Sigma.u%*%t(D)%*%Vi%*%res
  # Small area mean
  est_mean<-x%*%Beta.hat+D%*%u.hat

  ### Estimation shrinkage
  Di<-rep(NA,areanumber)
  Di[as.logical(ni_ind)]<-vardir
  Di[as.logical(abs(ni_ind-1))]<-(1/(4*mean(area_count)))
  Bi.tot<-as.numeric(Di)/(estsigma2u+as.numeric(Di))

  ### Out-of-sample-prediction
  dataframe_pop_aux<-data.frame(dataframe_pop_aux,helper=rnorm(1,0,1))
  lhs(formula)<-quote(helper)
  daten.pop<-makeXY(formula=formula,data=dataframe_pop_aux)
  x.total<-daten.pop$x

  ### Synthetic prediction for out-of-sample
  pred_out<-as.numeric(as.vector(Beta.hat)%*%t(x.total))
  est_mean_final<-rep(NA,areanumber)

  ### Merge in-sample and out-of-sample (out and in sample)
  est_mean_final[as.logical(ni_ind)]<-est_mean[1:areanumber_sample]
  est_mean_final[as.logical(abs(ni_ind-1))]<-pred_out[as.logical(abs(ni_ind-1))]

  ### Temporary variable necessary for computing the CI
  est_mean_final.MSE<-est_mean_final

  ### Truncation
  est_mean_final[est_mean_final<0]<-0
  est_mean_final[est_mean_final>(pi/2)]<-(pi/2)

  ### Back-transformation
  est_mean_final<-(sin(est_mean_final))^2

  ### Extend random effects
  u.hat_new<-rep(0,areanumber)
  u.hat_new[as.logical(ni_ind)]<-u.hat


# Benchmarking Point Estimation --------------------------------------------------------------------

  est_mean_bench<-benchmark_mean(as.numeric(est_mean_final),T=bench_value,
                                 Phi=diag(ratio/est_mean_final), w=ratio, gamma=0,
                                 Omega=diag(x.total_saind))


# Confidence Interval Estimation -------------------------------------------------------------------

  Li<-rep(NA,areanumber)
  Ui<-rep(NA,areanumber)

  ### Bootstrap
  if(Boot==TRUE){
  ti<-matrix(NA,areanumber,B)
  boots_est<-matrix(NA,areanumber,B)
  boots_par<-matrix(NA,areanumber,B)

  for (b in 1:B){

    set.seed(b)

    v_boot<-rnorm(areanumber,0,sqrt(estsigma2u))
    e_boot<-rnorm(areanumber_sample,0,sqrt(vardir))
    Xbeta_boot<-pred_out

    ## Theta under transformation
    theta<-Xbeta_boot+v_boot

    ## Truncation
    true_value_boot<-Xbeta_boot+v_boot
    true_value_boot[true_value_boot<0]<-0
    true_value_boot[true_value_boot>(pi/2)]<-(pi/2)

    ## Back-transformation
    true_value_boot<-(sin(true_value_boot))^2
    boots_par[,b]<-true_value_boot
    boots_par[,b]<-Xbeta_boot+v_boot

    ystar<-Xbeta_boot[as.logical(ni_ind)]+v_boot[as.logical(ni_ind)]+e_boot

    ## Estimation of beta_boot
    if (method=="REML")ottimo_boot=suppressWarnings(optimize(logl,c(-5,100),maximum = TRUE,
                  vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=ystar,x_logl=x))
    if (method=="AR")ottimo_boot=suppressWarnings(optimize(ARlogl,c(-5,100),maximum = TRUE,
                  vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=ystar,x_logl=x))
    if (method=="AP")ottimo_boot=suppressWarnings(optimize(APlogl,c(-5,100),maximum = TRUE,
                  vardir_logl=vardir,areanumber_logl=areanumber_sample,direct_logl=ystar,x_logl=x))

    estsigma2u_boot<-ottimo_boot$maximum

    ## Computation of the coefficients'estimator (Bstim)
    D=diag(1,areanumber_sample)
    V<-estsigma2u_boot*D%*%t(D)+diag(as.numeric(vardir))
    Vi<-solve(V)
    Q<-solve(t(x)%*%Vi%*%x)
    Beta.hat_boot<-Q%*%t(x)%*%Vi%*%ystar

    ## Computation of the EBLUP
    res<-ystar-c(x%*%Beta.hat_boot)
    Sigma.u=estsigma2u_boot*I
    u.hat=Sigma.u%*%t(D)%*%Vi%*%res

    ## Small area mean
    est_mean_boot<-x%*%Beta.hat_boot+D%*%u.hat
    Bi<-as.numeric(vardir)/(estsigma2u_boot+as.numeric(vardir))
    #Bi.out<-(1/(4*mean(area_count)))/((1/(4*mean(area_count)))+estsigma2u_boot)
    ti[as.logical(ni_ind),b]<-(theta[as.logical(ni_ind)]
                               -est_mean_boot)/sqrt(as.numeric(vardir)*(1-Bi))

    ## Synthetic prediction for out-of-sample
    pred_out_boot<-as.numeric(as.vector(Beta.hat_boot)%*%t(x.total))
    ti[as.logical(abs(ni_ind-1)),b]<-(theta[as.logical(abs(ni_ind-1))]
                                    -pred_out_boot[as.logical(abs(ni_ind-1))])/sqrt(estsigma2u_boot)
    #ti[as.logical(abs(ni_ind-1)),b]<-(theta[as.logical(abs(ni_ind-1))]-pred_out_boot[as.logical(abs(ni_ind-1))])/sqrt(as.numeric((1/(4*mean(area_count))))*(1-Bi.out))
    print(b)
  } # End of bootstrap runs

  qi<-matrix(NA,areanumber,2)

  for (i in 1:areanumber)
  {
   qi[i,1]<-quantile(ti[i,],prob=alpha/2,na.rm=T)
   qi[i,2]<-quantile(ti[i,],prob=(1-alpha/2),na.rm=T)
  }

  Li<-matrix(NA,areanumber,1)
  Ui<-matrix(NA,areanumber,1)

  Li<-(est_mean_final.MSE+qi*sqrt(Di*(1-Bi.tot)))[,1]
  Ui<-(est_mean_final.MSE+qi*sqrt(Di*(1-Bi.tot)))[,2]

  ### Truncation
    Li[Li<0]<-0
    Ui[Ui>(pi/2)]<-(pi/2)

  ### Back-transformation
    Li<-(sin(Li))^2
    Ui<-(sin(Ui))^2

  } # End of Confidence Interval Estimation


# Generating Output Frame --------------------------------------------------------------------------

  list(coefficients = Beta.hat,sigma2v=estsigma2u,samplesize = ni,rand.eff=u.hat_new,
       est_mean_FH=data.frame(area = x.total_saind,pred = as.numeric(est_mean_final)),
       est_mean_FH_bench=data.frame(area = x.total_saind,pred = as.numeric(est_mean_bench)),
       CI_FH = data.frame(area = x.total_saind, Low=as.numeric(Li),Up=as.numeric(Ui))
  )

}

