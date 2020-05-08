


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




