
# Task 1

# a)
a0=2;
b0=2;
s=5;
n=20;
f=n-s;
numberOfSamples=10000;
alpha=a0+s;
beta=b0+f;
betaSamples=rbeta(numberOfSamples,alpha,beta,ncp=0)

meanArray=rep(0,numberOfSamples);
prevMeanValue=0;

varArray=rep(0,numberOfSamples);

for(i in 1:numberOfSamples){

  meanArray[i]=(prevMeanValue*(i-1)+betaSamples[i])/i;
  prevMeanValue=meanArray[i];
  
  if (i>1){
    varArray[i]=var(betaSamples[1:i])
  }

}


convergingMean=(alpha)/(alpha+beta);
convergingStd=sqrt(alpha*beta/((alpha+beta)^2*(alpha+beta+1)))
  
plot(1:numberOfSamples,meanArray, type="l",col="red")
points(1:numberOfSamples,rep(convergingMean,numberOfSamples),type="l", col="purple")

plot(2:numberOfSamples,sqrt(varArray[2:numberOfSamples]), type="l",col="red")
points(2:numberOfSamples,rep(convergingStd,numberOfSamples-1),type="l", col="purple")


# b)

numberOver=ifelse(betaSamples>0.3,1,0)
probabilityOver=sum(numberOver)/length(numberOver)
trueProbability=1-pbeta(0.3, alpha, beta)

# c)
logOdds=log(betaSamples/(1-betaSamples))
hist(logOdds)
density(logOdds)
min(logOdds)
fit=density(logOdds)
plot(fit)


