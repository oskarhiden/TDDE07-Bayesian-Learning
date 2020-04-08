
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

