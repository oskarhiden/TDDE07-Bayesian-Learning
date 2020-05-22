
# 1a
mu = 10
sigma2 = 2
T_1 = 200
x_1 = mu

AR_process = function(mu, phi, simga2, x_1, T_1){
  x = rep(0,T_1)
  x[1] = x_1
  
  for (i in 2:T_1) {
    error = rnorm(1,0,sqrt(sigma2))
    x[i] = mu+phi*(x[i-1]-mu)+error
  }
  return(x)
}

phi_1 = AR_process(mu, 1, simga2, x_1, T_1)
phi_05 = AR_process(mu, 0.5, simga2, x_1, T_1)
phi_0 = AR_process(mu, 0, simga2, x_1, T_1)
phi_neg_05 = AR_process(mu, -0.5, simga2, x_1, T_1)
phi_neg_1 = AR_process(mu, -1, simga2, x_1, T_1)

plot(phi_1, type ="l")
plot(phi_05, type="l")
plot(phi_0, type="l")
plot(phi_neg_05, type="l")
plot(phi_neg_1, type="l")

# Phi decides how the previous draw effects the next draw of x.Phi = 1 gives equlal weight to all previous error draws 
# while Phi= 0.5 only careses about the moast preavious in declining order, exponential decrease. For negative values of Phi the 
# previous deviation from mean is compensated in the other direction. Over mean then becomes under mean.


#1b
library(rstan)
set.seed(12345)
phi_03 = AR_process(mu, 0.3, simga2, x_1, T_1)
phi_095 = AR_process(mu, 0.95, simga2, x_1, T_1)

StanModel = '
data{
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
}
model {
  for (n in 2:N)
    y[n] ~ normal(mu + phi * y[n-1], sqrt(sigma2));
}
'

data_03 = list(N=T_1, y=phi_03)
burnin = 1000
niter = 2000
fit_03 = stan(model_code=StanModel,data=data_03,
           warmup=burnin,iter=niter,chains=4)

print(fit_03,digits_summary=3) 

data_095 = list(N=T_1, y=phi_095)
fit_095 = stan(model_code=StanModel,data=data_095,
              warmup=burnin,iter=niter,chains=4)

print(fit_095,digits_summary=3) 

#not true values for mu, but phi and simga2 were able to predict. mu was futher from the real value when we use phi=0.95.

#ii
post_draws_03 = extract(fit_03)
plot(post_draws_03$mu, type="l")
plot(post_draws_03$phi, type="l")

post_draws_095 = extract(fit_095)
plot(post_draws_095$mu, type="l")
plot(post_draws_095$phi, type="l")

#convergence
plot(cumsum(post_draws_03$mu) / 1:length(post_draws_03$mu), type="l")
plot(cumsum(post_draws_03$phi) / 1:length(post_draws_03$phi), type="l")

plot(cumsum(post_draws_095$mu) / 1:length(post_draws_095$mu), type="l")
plot(cumsum(post_draws_095$phi) / 1:length(post_draws_095$phi), type="l")

#  All varibles seem to convergate but some thorwards the wrong value.

#joint posterior draws
plot(post_draws_03$phi,post_draws_03$mu, type="p")
plot(post_draws_095$phi,post_draws_095$mu, type="p")

# They seem to have a quite strong correlation.

#1c

campy = read.table("/Users/oskarhiden/Git/TDDE07 Bayesian Learning/Lab 4/campy.dat", header = TRUE)

campy = as.list(campy)

StanModel = '
data{
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
  real x[N];
}
model {
  mu ~ normal(2,50); // Normal with mean 2(from log10), st.dev. 50
  sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 10
  phi ~ uniform(-1, 1);
  for (n in 2:N){
    x[n] ~ normal(mu + phi * x[n-1], sqrt(sigma2));
    y[n] ~ poisson(exp(x[n]));
  }
}
'
T_2 = 140
campy = campy$c
data = list(N=T_2, y=campy)
burnin = 1000
niter = 2000
fit_campy = stan(model_code=StanModel,data=data,
              warmup=burnin,iter=niter,chains=4)

print(fit_campy, digits_summary = 3)

post_summary = summary(fit_campy)
mean = post_summary$summary[,1]
low = post_summary$summary[,4]
high = post_summary$summary[,8]

plot(campy)
points(exp(mean[4:(length(mean)-1)]), col="red", type="l")
points(exp(low[4:(length(mean)-1)]), col="blue", type="l")
points(exp(high[4:(length(mean)-1)]), col="blue", type="l")


# 1d

StanModel_inform = '
data{
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
  real x[N];
}
model {
  mu ~ normal(2,50); // Normal with mean 2(from log10), st.dev. 50
  sigma2 ~ scaled_inv_chi_square(300,0.1); // Scaled-inv-chi2 with nu 1, sigma 10
  phi ~ uniform(-1, 1);
  for (n in 2:N){
    x[n] ~ normal(mu + phi * x[n-1], sqrt(sigma2));
    y[n] ~ poisson(exp(x[n]));
  }
}
'

fit_campy = stan(model_code=StanModel_inform,data=data,
                 warmup=burnin,iter=niter,chains=4)

print(fit_campy, digits_summary = 3)

post_summary = summary(fit_campy)
mean = post_summary$summary[,1]
low = post_summary$summary[,4]
high = post_summary$summary[,8]

plot(campy)
points(exp(mean[4:(length(mean)-1)]), col="red", type="l")
points(exp(low[4:(length(mean)-1)]), col="blue", type="l")
points(exp(high[4:(length(mean)-1)]), col="blue", type="l")

