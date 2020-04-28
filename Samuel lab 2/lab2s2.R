filename='C:/Users/samue/Documents/LIU/TDDE07/LABS/TDDE07-Bayesian-Learning/Samuel lab 2/WomenWork.txt'
Data<-read.table(filename,head=TRUE)
y=(as.vector(Data[,1]))
x=as.matrix(Data[,-1])
numberOfParameters=length(x[1,])
library(mvtnorm)

#a)
logPosteriorProbability=function(beta,y,x,mu,sigma){
  
  lengthBeta=length(beta)
  linearPrediction=x%*%beta
  
  logMaximumLikelihood=sum(y*linearPrediction-log(1+exp(linearPrediction)))
  

  logPrior= dmvnorm(beta, matrix(0,lengthBeta,1), sigma, log=TRUE);
  return(logMaximumLikelihood+logPrior)
  
}

tau=10
betaStart=as.vector(rep(0,numberOfParameters))
priorMean=as.vector(rep(0,numberOfParameters))
priorStd=tau^2*diag(numberOfParameters)


OptimResults<-optim(betaStart,logPosteriorProbability,gr=NULL,y,x,priorMean,priorStd,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

postMode <- OptimResults$par

postCov <- -solve(OptimResults$hessian)

y=qnorm(c(0.025,0.975),postMode[7],postCov[7,7])
#Yes it is important. If we change the nmber of childs from 0 to 1 of one woman when 
#the linear prediction is 0 (were the derivative of the logistic is at its maximum)
#we will get a probability change from 50% to 73 % and the absolute value of this variable
# is with 95% confidence abow zero (two sided)


#b)
library(MASS)
numberOfSamples=100
simulatedBetas=mvrnorm(numberOfSamples,postMode,postCov)
xForB=c(1, 10, 8, 10, 1, 40, 1, 1)
linearPrediction=xForB%*%t(simulatedBetas)
yPredictions=exp(linearPrediction)/(1+exp(linearPrediction))
hist(yPredictions)

#c

numberOfSamples=1000
simulatedBetasC=mvrnorm(numberOfSamples,postMode,postCov)
linearPredictionC=xForB%*%t(simulatedBetas)
yPredictions=exp(linearPrediction)/(1+exp(linearPrediction))
allPredictions=c()
for (i in yPredictions) {
  allPredictions=cbind(allPredictions,rbinom(1,10,i))  
}

hist(allPredictions)

