
QualityMeasure <- function(True.mean, Est.mean, MSE, MSETF) {

ResultsAnaModelBasedMSE<-function (True.mean, Est.mean, MSE){

  m<-dim(True.mean)[1]
  NoSim<-dim(True.mean)[2]
  True.mean<-t(True.mean)
  Est.mean<-t(Est.mean)
  MSE<-t(MSE)


  assign("m",m,pos=1)
  RB <-rep(0,m)
  Bias <-rep(0,m)
  RRMSE <-rep(0,m)

  True.MSE <-rep(0,m)
  Est.MSE<-rep(0,m)
  CV.MSE<-rep(0,m)
  BIAS.MSE<-rep(0,m)

  for (i in 1:m){
    RB[i]<-mean(((Est.mean[,i])-(True.mean[,i]))/(True.mean[,i]))
    Bias[i]<-mean(((Est.mean[,i])-(True.mean[,i])))
    RRMSE[i]<-sqrt(mean((((Est.mean[,i])-(True.mean[,i]))/(True.mean[,i]))^2))
    True.MSE[i]<-mean((Est.mean[,i]-True.mean[,i])^2)
    Est.MSE[i]<-mean(MSE[,i])



    CV.MSE[i]<-(sqrt(mean((sqrt(MSE[,i])-sqrt(True.MSE[i]))^2))/sqrt(True.MSE[i]))*100
    BIAS.MSE[i]<-mean(MSE[,i]-True.MSE[i])/True.MSE[i]
  }
  True.RootMSE=sqrt(True.MSE)
  Est.RootMSE=sqrt(Est.MSE)
  RB.estRRMSE=((Est.RootMSE-True.RootMSE)/True.RootMSE)

  CR<-rep(0,m)
  CR_Length<-rep(0,m)
  coverage<-matrix(0,NoSim,m)
  Length<-matrix(0,NoSim,m)


  for (k in 1:NoSim)
  {
    for (oo in 1:m)
    {Length[k,oo]<-2*1.96*sqrt(MSE[k,oo])
     if ((Est.mean[k,oo]-1.96*sqrt(MSE[k,oo])< True.mean[k,oo] && True.mean[k,oo] < Est.mean[k,oo]+1.96*sqrt(MSE[k,oo]) ))coverage[k,oo]<-1
     else
       coverage[k,oo]<- 0}
  }

  for (i in 1:m){
    CR[i]<-mean(coverage[,i])
    CR_Length[i]<-mean(Length[,i])
  }


  list(Bias=Bias,RB=RB,RRMSE=RRMSE,True.RMSE=True.RootMSE,Est.RMSE=Est.RootMSE,RB_RMSE=RB.estRRMSE,RRMSE_RMSE=CV.MSE,Coverage=CR,CI_Length=CR_Length)
}

ResultsAnaModelBased<-function (True.mean, Est.mean){

  m<-dim(True.mean)[1]
  NoSim<-dim(True.mean)[2]
  True.mean<-t(True.mean)
  Est.mean<-t(Est.mean)


  assign("m",m,pos=1)
  RB <-rep(0,m)
  ARB <-rep(0,m)
  RRMSE <-rep(0,m)
  Bias <- rep(0,m)
  True <- rep(0,m)
  Dir <- rep(0,m)

  True.MSE <-rep(0,m)
  Est.MSE<-rep(0,m)
  CV.MSE<-rep(0,m)
  BIAS.MSE<-rep(0,m)

  for (i in 1:m){
    Bias[i] <- mean((Est.mean[,i])-(True.mean[,i]))
    True[i] <- mean(True.mean[,i])
    Dir[i] <- mean(Est.mean[,i])
    RB[i] <- mean(((Est.mean[,i])-(True.mean[,i]))/True.mean[,i])
    ARB[i] <- mean(abs(((Est.mean[,i])-(True.mean[,i])))/abs(True.mean[,i]))
    RRMSE[i] <- sqrt(mean((((Est.mean[,i])-(True.mean[,i]))/(True.mean[,i]))^2))
    True.MSE[i] <- mean((Est.mean[,i]-True.mean[,i])^2)
  }

  True.RootMSE=sqrt(True.MSE)



  list(Bias=Bias,True=True,Dir=Dir,RB=RB,RRMSE=RRMSE,True.RMSE=True.RootMSE)
}

if (MSETF==TRUE) {
  results <- ResultsAnaModelBasedMSE(True.mean, Est.mean, MSE)
  return(results)
}

if (MSETF==FALSE) {
  results <- ResultsAnaModelBased(True.mean, Est.mean)
  return(results)
}
}

