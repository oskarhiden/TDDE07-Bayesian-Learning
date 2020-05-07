


#1a
rainfall = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 3/rainfall.dat")
rainfall = as.vector(rainfall[,1])
x = rainfall
# mu Normaly distributed
mu_0 = 1
tau2_0 = 1 

# sigma scaled inv chi2
v_0 = 0
sigma2_0 = 1

sigma2_draw = 0.001
mu_draw = 1
n = length(x)
x_hat = mean(x)

iterations = 1000
samples = matrix(nrow = iterations, ncol = 2)

for (i in 1:iterations) {

# full conditional posteriors
# mu
w = (n/sigma2_draw)/(n/sigma2_draw + 1/tau2_0)
mu_n = w*x_hat + (1-w)*mu_0
tau2_n = 1/(n/sigma2_draw+1/tau2_0)

mu_draw = rnorm(1,mu_n, sigma2_draw)
# simga2
v_n = v_0+n
sigma2_n = (v_0*sigma2_0 + sum((x-mu_draw)^2))/v_n

chi_draw = rchisq(1, v_n)
sigma2_draw = v_n*sigma2_n/chi_draw

samples[i,1] = mu_draw
samples[i,2] = sigma2_draw

}
samples
