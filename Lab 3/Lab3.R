


#1a
rainfall = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 3/rainfall.dat")
rainfall = as.vector(rainfall[,1])
x = rainfall
# mu Normaly distributed
mu_0 = 25
tau2_0 = 1500

# sigma scaled inv chi2
v_0 = 0
sigma2_0 = 1

sigma2_draw = 1
mu_draw = 1
n = length(x)
x_hat = mean(x)

iterations = 10000
samples = matrix(nrow = iterations, ncol = 2)

for (i in 1:iterations) {
# full conditional posteriors
# mu
w = (n/sigma2_draw)/(n/sigma2_draw + 1/tau2_0)
mu_n = w*x_hat + (1-w)*mu_0
tau2_n = 1/(n/sigma2_draw+1/tau2_0)

mu_draw = rnorm(1,mu_n, tau2_n)

# simga2
v_n = v_0+n
sigma2_n = (v_0*sigma2_0 + sum((x-mu_draw)^2))/v_n

chi_draw = rchisq(1, v_n)
sigma2_draw = v_n*sigma2_n/chi_draw

samples[i,1] = mu_draw
samples[i,2] = sigma2_draw

}
samples
plot(samples, type="b")
hist(samples[,1], breaks = 30)
hist(samples[,2], breaks = 30)
mu_hat = mean(samples[,1])
sigma2_hat = mean(samples[,2])
#1b
# Function that simulates from the Scaled Inv Chi 2
scaled_inv_chi2 <- function(nr_draws, df, scale){
  return((df*scale)/rchisq(nr_draws,df=df))
}

# Function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}


#1c
dens = density(x)
plot(dens)

lines(xGrid, dnorm(xGrid, mean = mu_hat, sd = sqrt(sigma2_hat)), type = "l", lwd = 2, col = "blue")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
legend("topright", box.lty = 1, legend = c("Kernel Density","Gibs normal","Mixted nomral"), col=c("black","red","blue"), lwd = 2)


# 2a
ebay = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 3/eBayNumberOfBidderData.dat", header = TRUE)

data = ebay[,-2]
model = glm(nBids~PowerSeller+VerifyID +Sealed +Minblem + MajBlem + LargNeg + LogBook + MinBidShare,family = poisson, data=data)
print(model$coefficients)

sort(abs(exp(model$coefficients)-1))
# MinBidShare have the biggest impact, and Sealed the second. PowerSeller seems to not impact. 


# 2b
library(mvtnorm)
X = as.matrix(ebay[,2:10])

y = ebay[,1]
prior_cov = t(X)%*%X
log_posterior_prob = function(beta, X, y){
  log_like=0
  for (i in 1:dim(X)[1]) {
    lambda = exp(t(X[i,])%*%beta)
    log_like = log_like +log(dpois(y[i], lambda = lambda))
  }
  
    
  prior_log_prob = log(dmvnorm(beta, mean = rep(0,nr_param), sigma = prior_cov))
  tot_log_like = dim(X)[1]*prior_log_prob + log_like
  return(tot_log_like)
}

log_faculty = function(y){
  log_fac = 0
  if (y != 0) {
    for (i in 1:y) {
      log_fac = log_fac + log(i) 
    }
  }
  
  return(log_fac) 
}

log_posterior_prob = function(beta, X, y){
  log_like=0
  for (i in 1:dim(X)[1]) {
    x_b = t(X[i,])%*%beta
    part = x_b*y[i]-exp(x_b)-log_faculty(y[i])
    log_like = log_like + part
  }
  
  
  prior_log_prob = log(dmvnorm(beta, mean = rep(0,nr_param), sigma = prior_cov))
  tot_log_like = dim(X)[1]*prior_log_prob + log_like
  return(tot_log_like)
}
#startvalues for beta vector
nr_param = dim(X)[2]
init_beta = as.vector(rep(0,nr_param))
result = optim(init_beta, log_posterior_prob, gr=NULL, X, y, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

beta_hat = result$par
j_y = -result$hessian
post_cov = solve(j_y)

beta_hat
post_cov


#2c
library(MASS)
library(LaplacesDemon)

proposal_generator = function(theta_prev, c, cov_mat){
  theta_p = mvrnorm(1, theta_prev, c*cov_mat)
  return(theta_p)
}

alpha_generator = function(proposal, theta_before, log_post_func, ...){
  
  new_prob = log_post_func(proposal,...)
  old_prob = log_post_func(theta_before, ...)
  
  alpha = exp(new_prob - old_prob)
  
  return(min(alpha, 1))
}

metropolis = function(init_theta, c, nr_iter, cov_mat,log_post_func, ...){
  
  metrop_sim = matrix(data=1, nrow = nr_iter, ncol = length(init_theta))
  metrop_sim[1,] = init_theta
  
  for (i in 2:nr_iter) {
    proposal = proposal_generator(metrop_sim[i-1,],c,cov_mat);
    
    #calculate apha
    alpha = alpha_generator(proposal, metrop_sim[i-1,],log_post_func, ...)
    #accept the proposal with proability alpha
    if( rbern(1, alpha)){
      metrop_sim[i,] = proposal
    }else{
      metrop_sim[i,] = metrop_sim[i-1,]
    }
  }
  return(metrop_sim)
}
init_theta = rep(0,length(beta_hat))
c=1
nr_iter=3000
cov_mat = post_cov

#proposal_generator(init_theta, c, cov_mat)
metrop_samp = metropolis(init_theta, c, nr_iter, cov_mat,log_posterior_prob,X,y)

hist(metrop_samp[,1])
hist(metrop_samp[,2])
hist(metrop_samp[,3])
hist(metrop_samp[,4])
hist(metrop_samp[,5])
hist(metrop_samp[,6])
hist(metrop_samp[,7])
hist(metrop_samp[,8])
hist(metrop_samp[,9])
beta_hat

#y1 = cumsum(metrop_samp[,1]) / seq_along(metrop_samp[,1]) 
y1 = cumsum(metrop_samp[,1]) / 1:length(metrop_samp[,1]) 
y2 = cumsum(metrop_samp[,2]) / 1:length(metrop_samp[,2])
y3 = cumsum(metrop_samp[,3]) / 1:length(metrop_samp[,3])

#convergence plotted with beta_hat, beta by optimization from maximum posterior liklehood. 
plot(rep(beta_hat[1],length(metrop_samp[,1])), type="l", main="beta_1", ylim = c(0,1.2))
lines(y1, type="l", col="blue")

plot(rep(beta_hat[2],length(metrop_samp[,1])), type="l", main="beta_2", ylim = c(-0.12,0))
lines(y2, type="l", col="blue")

plot(rep(beta_hat[3],length(metrop_samp[,1])), type="l", main="beta_3",ylim = c(-0.4,0.2))
lines(y3, type="l", col="blue")

# 2d
x_new = c(1,1,1,1,0,0,0,1,0.5)

x_b_pred = metrop_samp%*%x_new
lambda = exp(x_b_pred)
hist(exp(lambda))

min(lambda)


#poisson for all observed lambda to get y
y_prob_0 = rep(0,nr_iter)

for (i in 1:nr_iter) {
  
  y_prob_0[i] = dpois(0,lambda[i])
  
}  
sum(y_prob_0)/nr_iter

