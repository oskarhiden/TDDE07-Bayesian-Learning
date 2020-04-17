
observations=c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
n=length(observations)
mu=3.7
tauSqr=0
NDraws=10000;

tauSqr=sum((log(observations)-mu)^2)/length(observations)


PostDraws<-(n*tauSqr)/rchisq(NDraws,n)

h=hist(PostDraws,breaks=100, plot=FALSE)
h$counts=h$counts*50/sum(h$counts)
plot(h, col="red",ylim=c(0,10))

test = dinvchisq(seq(0.02, 1, by=0.01), n, tauSqr)
points(seq(0.02, 1, by=0.01),test,col="blue",type="l")



#Task 2

G = 2*pnorm(sqrt(inv_chi)/sqrt(2))-1
h=hist(G,breaks=100, plot=FALSE)
plot(h, col="red")


#Task 3

sort_g = sort(G)
low_tail = sort_g[length(G)*0.5]
high_tail = sort_g[length(G)*0.95]
c(low_tail, high_tail)

densityG=density(G)
plot(densityG)
library(HDInterval)
hdi(densityG, credMass=0.90)

