

data = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 2/TempLinkoping.txt")

time = as.numeric(as.vector(data[-1,1]))
temp = as.numeric(as.vector(data[-1,2]))

# 2a)
mu_0 = c(-10, 100, -100)
omega_0 = 0.01*diag(3)
omega_0_inv = 100*diag(3)
v_0 = 4
sigma_0_2 = 1


#draws sigma 2
nr_draws = 100
x_draws = rchisq(nr_draws,v_0)
sigma2_draws = ((v_0)*sigma_0_2)/x_draws


#draws beta given simga
library(mvtnorm)
betas = matrix(, nrow=nr_draws, ncol=3)
for (i in 1:nr_draws) {
  covar = sigma2_draws[i]*omega_0_inv
  betas[i,] = rmvnorm(1, mu_0, covar)
}
betas
A = rep(1, length(time))
B = time
C = time^2
x = cbind(A, B, C) #beta order B0, B1, B2
x
t(x)
y_draws = betas%*%t(x) 
y_draws

#draw predictions
plot(time*365, y_draws[1,], type = "l", ylim = c(-60,70))

for (i in 2:nr_draws) {
  points( y_draws[i,], type="l")
}
# hottest in the summer


#2b
#posterior
library(matlib)
beta_hat = inv( t(x)%*%x ) %*% (t(x)%*%temp)
mu_n = inv( t(x)%*%x+omega_0 )%*%( t(x)%*%x%*%beta_hat + omega_0%*%mu_0)
omega_n = t(x)%*%x+omega_0
v_n = v_0 + length(time)
sigma_n_2 = (v_0*sigma_0_2 + (t(temp)%*%temp + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n))/v_n

#draw sigma2
nr_draws = 100
sigma_draws = rchisq(nr_draws,v_n) 
sigma2_draws = ((v_n)*sigma_n_2)/sigma_draws

betas_post = matrix(, nrow=nr_draws, ncol=3)
for (i in 1:nr_draws) {
  covar = sigma2_draws[i]*inv(omega_n)
  betas_post[i,] = rmvnorm(1, mu_n, covar)
}

#plot hist of beta and sigma
hist(sigma2_draws, breaks = 50)
hist(betas[,1], breaks = 50)
hist(betas[,2], breaks = 50)
hist(betas[,3], breaks = 50)

# calculate y, predicted temp for all draws
y_post = betas_post%*%t(x) 
y_post

#Calculate mean for y each day form all predictions
#y_mean = colMeans(y_post)
y_median = apply(y_post,2,median)
plot(temp, col="blue")
points(time*365, y_median, type="l")

#curves for 95% equal tail
#sorted_y_post = sort(y_psot)
sorted_y_post = apply(y_post,2,sort,decreasing=F)
points(sorted_y_post[round(nr_draws*0.025),], type="l", col="red")
points(sorted_y_post[round(nr_draws*0.975),], type="l", col="red")


quantile_y_post = apply(y_post,2,quantile, probs=c(0.025,0.975))
points(quantile_y_post[1,], col="green", type="l")
points(quantile_y_post[2,], col="green", type="l")
# The posterior probability intevall showes us were our real model would be wthin wiht 95% credability. (from beta)
# It does not show were the actual y is by 95% credability. 

#2c
x_hat = which.max(y_median)
y_median[x_hat]

#From previous data
all_hot_days = apply(y_post, 1, which.max)
all_hot_days
hist(all_hot_days, breaks=10)

#from derivation
all_hot_days2 = (-betas_post[,2]/(2*betas_post[,3]))*365
all_hot_days2

#2d
# We will have a laplace distribution for the betas with u0=0 and omega0=I(8)*lambda. Were lambda is the smoothing 
# coefficient that decides how much the betas are allowed to differ from zero and u0 sets the mean of these to be zero 
# in the beginning. This will make a few betas go into the fat tails of the laplace distribution and a few to end up 
# at 0, the prior mean.

# 3a
women = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 2/WomenWork.dat", header = TRUE)
women = as.numeric(women[-1,])

y = as.vector(women[,1])
X = as.matrix(women[,-1])
nr_param = length(X[1,])

# Prior values
mu = rep(0, nr_param)
tau = 10
sigma = tau^2*diag(nr_param)

posterior_loglike = function(beta, y, X, mu, sigma){
  nr_param = length(beta)
  x_B = X%*%beta
  
  #loglike of
  log_like = sum( x_B*y -log(1 + exp(x_B)))
  if (abs(log_like) == Inf) logLik = -20000 # constraint for logLike
  
  # evaluating the prior
  log_prior = dmvnorm(beta, matrix(0,nr_param,1), sigma, log=TRUE)
  
  # add the log prior and log-likelihood together to get log posterior
  return(log_like + log_prior)
  
}

#startvalues for beta vector
init_beta = as.vector(rep(0,nr_param))

#Maximising loglike by changing beta
optim_results<-optim(init_beta, posterior_loglike, gr=NULL,y,X,mu,sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

#Result
post_beta = optim_results$par
j_y = -optim_results$hessian # Posterior covariance matrix is -inv(Hessian)
post_cov = solve(j_y) #inverse

#nSmallChild = idex 7
b_sc = post_beta[7]
sigma2_sc = post_cov[7,7]

#draw from nrom
quantile_sc = dnorm(c(0.025,0.975), b_sc, sigma2_sc)
quantile_sc

#2b
library(mvtnorm)
#Predicting for
husband_inc = 10
years_ed = 8
years_exp = 10
age = 40
nr_small_ch = 1
nr_big_ch = 1
X_pred = c(1, husband_inc, years_ed, years_exp, (years_exp/10)^2, age, nr_small_ch, nr_big_ch )

#draw beta from posterior
nr_draws = 100
beta_draws = rmvnorm(100, post_beta, post_cov)
exp_xb = exp(X_pred%*%t(beta_draws))
pred_y_given_x = exp_xb/(1+exp_xb)

hist(pred_y_given_x)

#2c
nr_draws = 100
beta_draws = rmvnorm(nr_draws*10, post_beta, post_cov)
exp_xb = exp(X_pred%*%t(beta_draws))
pred_y_given_x = exp_xb/(1+exp_xb)

set.seed(12345)
pred = rbinom(nr_draws, 10, pred_y_given_x)
pred

# TEST
set.seed(12345)
allPredictions=c()
for (i in pred_y_given_x) {
  allPredictions=cbind(allPredictions,rbinom(1,10,i))  
}
hist(allPredictions)
# TEST
pred = rbinom(nr_draws, 10, c(0.1, 0.2))
hist(pred)
