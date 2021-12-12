#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#simulate data
n = 100
p = 1
b_var <- 1
# b <- rexp(p, rate = 1 / sqrt(b_var)) * (1-rbinom(p,1,0.5)*2)
b <- 5

x <- runif(n,0, 100)
sigma2 <- 0.01
e <- rnorm(n, sd = sqrt(sigma2))
a = 0
y <- a + x* b + e

#plot
plot(x,y)

#fit bayesian model
d <- list(n = n,
          x = x,
          y = y,
          prior_var_b = 0.1)

stan_program <- "
data {
  int<lower=0> n;
  vector[n] x;
  vector[n] y;
  real<lower=0> prior_var_b;
}
parameters {
  real b;
  real<lower=0> sigma2;
  real a;
}
model {
  a ~ normal(0,10);
  sigma2 ~ exponential(0.1);
  b ~ normal(0, prior_var_b);
  
  y ~ normal(a + x * b, sqrt(sigma2));
  
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
samps <- data.frame(as_draws_df(out$draws()))
pairs(samps[,2:4])
out

mean(samps$a)
a

mean(samps$sigma2)
sigma2

