filename='C:/Users/samue/Documents/LIU/TDDE07/LABS/TDDE07-Bayesian-Learning/Samuel lab 2/TempLinkoping.txt'
Data<-read.table(filename)
time=as.numeric(as.vector(Data[-1,1]))
temp=as.numeric(as.vector(Data[-1,2]))
u0=c(-10, 100, -100)
omega0=diag(x=0.01,3,3)
v0=4
sigma0sqr=1
numberOfDraws=10;

library(matlib)
library(mvtnorm)
x_draw = rchisq(numberOfDraws,v0)
inv_chi_draws = ((v0)*sigma0sqr)/x_draw

beta_draws=matrix(, nrow = numberOfDraws, ncol =3)
invOmega0=inv(omega0)

for (i in 1:numberOfDraws) {
  beta_draws[i,]=rmvnorm(1,u0,inv_chi_draws[i]*invOmega0)
}

ones=rep(1,length(time))

xmatrix=cbind(ones,time,time^2)
transxmatrix=t(xmatrix)
y=beta_draws%*%transxmatrix

plot(time*365,y[1,],type='l',col="red",ylim=c(-40,70))
for(i in 2:numberOfDraws){
  points(time*365,y[i,],col="red",type='l')
}

#C

omegaN=t(xmatrix)%*%xmatrix+omega0
invOmegaN=inv(omegaN)
vn=v0+length(time)

betahatt=inv((t(xmatrix)%*%xmatrix))%*%t(xmatrix)%*%temp
un=inv(t(xmatrix)%*%xmatrix+omega0)%*%(t(xmatrix)%*%xmatrix%*%betahatt+omega0%*%u0)
sigmasqrn=1/vn*(v0*sigma0sqr+t(temp)%*%temp+t(u0)%*%omega0%*%u0-t(un)%*%omegaN%*%un)

numberOfDraws2=100
sigmaSqrPostDraws=rchisq(numberOfDraws2,vn)

invChiPostDraws = (vn*sigmasqrn)/sigmaSqrPostDraws

betaDrawsPost=matrix(, nrow = numberOfDraws2, ncol =3)

for (i in 1:numberOfDraws2) {
  betaDrawsPost[i,]=rmvnorm(1,un,invChiPostDraws[i]*invOmegaN)
}

hist(invChiPostDraws,breaks=50)
hist(betaDrawsPost[,1],breaks = 50)
hist(betaDrawsPost[,2],breaks=50)
hist(betaDrawsPost[,3],breaks=50)

ypost=betaDrawsPost%*%t(xmatrix)
y_post_medians = apply(ypost,2,median,decreasing=F)


plot(time*365,y_post_medians,type='l',col="red",ylim=c(-25,50))
points(time*365, temp, col="pink")

quantile_y_post = apply(ypost,2,quantile, probs=c(0.025,0.975))
points(quantile_y_post[1,], col="green", type="l")
points(quantile_y_post[2,], col="green", type="l")

#Interval doesnt have all the points because the
# calculation is for the betas and not for the y values.
#This means that the model is within those barriers with 95 % confidence however that doesn't mean the observations are

#C

xMax=rep(1,length(numberOfDraws2))
xMaxbetaDrawsPost=-betaDrawsPost[,2]/(2*betaDrawsPost[,3])*365

hist(xMaxbetaDrawsPost)

#d
# We will have a laplace distribution for the betas with u0=0 and omega0=I(8)*lambda. Were lambda is the smoothing coefficient that decides
# how much the betas are allowed to differ from zero and u0 sets the mean of these to be zero in the beginning. This will make a few betas go
# inteo the fat tails of the laplace distribution and a few to end up at the prior mean.

#2


