
#Assignment 1a
post_beta = 2;
post_alpha = 2;
s = 5
f= 20-s
alpha = post_alpha+s
beta = post_beta+f

beta_mean = (alpha)/(alpha+beta);
beta_std = sqrt((alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1)))

n = 100000
random_numbers = rbeta(n, alpha, beta)
hist(random_numbers, breaks = 50)
mean_value = rep(0, 10000)
std_value = rep(0, 10000)

for(i in 1:10000){
  mean_value[i] = mean(random_numbers[1:i])
  std_value[i] =sqrt(var(random_numbers[1:i]))
}

plot(mean_value, type="l")
points(x=1:10000, y=rep(beta_mean,10000), type = "l", col="red")

plot(std_value, type="l")
points(x=1:10000, y=rep(beta_std, 10000), type="l", col="red")


#Assignment 1b
count = ifelse (random_numbers > 0.3, 1, 0)
count = sum(count);
#compare
posterior_prob = count/n
exact_value = 1-pbeta(0.3, alpha, beta)
posterior_prob
exact_value

#Assignment 1c
log_ods = log(random_numbers/(1-random_numbers))
samp = mean(log_ods)
hist(log_ods)
fit = density(log_ods, kernel="gaussian")
plot(fit)
print(fit)
#mean = -1,.2413
#std = 1.0848

#Assignment 2a
y = c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
my = 3.7
n=10
tau = sum((log(y)-my)^2)/length(y)
tau  

library(LaplacesDemon)

x_draw = rchisq(10000,n)
inv_chi = ((n)*tau)/x_draw

h=hist(inv_chi,breaks=100, plot=FALSE)
h$counts=h$counts*50/sum(h$counts)
plot(h, col="red",ylim=c(0,10))

test = dinvchisq(seq(0.02, 1, by=0.01), n, tau)
points(seq(0.02, 1, by=0.01),test,col="blue",type="l")

#2b
G = 2*pnorm(sqrt(inv_chi)/sqrt(2))-1
hist(G, breaks=100)
#plot(density(G))

#2c
#90% by sorting.
sort_g = sort(G)
low_tail = sort_g[length(G)*0.5]
high_tail = sort_g[length(G)*0.95]
c(low_tail, high_tail)

fit_G = density(G)
plot(fit_G)

library(HDInterval)
hid_90 = hdi(fit_G, credMass=0.90)
hid_90

#3a
wind_d = c(40, 303, 326, 285, 296, 314, 20, 308, 299, 296)
wind_r = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

#?Bessel
k=seq(0.01, 7, by=0.01) #exp(lambda=1) as prior for k
my=2.39
#mises_distr =exp(k*cos(wind_r-my))/(2*pi*besselI(k,0)) 
p_y = rep(0, length(k))
for(i in 1:length(k)){
  p_y[i] = prod(exp(k[i]*cos(wind_r-my))/(2*pi*besselI(k[i],0)) )*dexp(k[i])
}
plot(k,p_y, type="l")

#3b
#find mode: moest frequent value of k. 
k[which.max(p_y)]
