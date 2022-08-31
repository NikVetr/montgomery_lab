#libraries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(mvtnorm)

#functions
invlogit <- function(x) exp(x)/(1+exp(x))
prop_greater_than_0 <- function(x) mean(x>0)
subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

#simulation parameters
p <- 5000
n <- 2000

#simulate data
theta <- rnorm(p)
theta_diff <- rep(0,p)
theta_diff <- rnorm(p, sd = 0.5)
theta_1 <- theta + theta_diff / 2
theta_2 <- theta - theta_diff / 2
prob_1 <- invlogit(theta_1)
prob_2 <- invlogit(theta_2)
x1 <- rbinom(p, n, prob_1)
x2 <- rbinom(p, n, prob_2)

#stan model
d <- list(n = n,
          p = p,
          x1 = x1,
          x2 = x2)

# base_1 = paste0("separate_group_thetas_estimated")
# stan_program_1 <- '
# data {
#     int<lower=1> n;
#     int<lower=1> p;
#     int<lower=0,upper=n> x1[p];
#     int<lower=0,upper=n> x2[p];
# }
# parameters {
#     vector[p] theta_1_raw;
#     vector[p] theta_2_raw;
#     
#     real theta_1_mu;
#     real theta_2_mu;
#     real<lower=0> theta_1_sd;
#     real<lower=0> theta_2_sd;
# }
# transformed parameters {
#     vector[p] theta_1 = theta_1_raw * theta_1_sd + theta_1_mu;
#     vector[p] theta_2 = theta_2_raw * theta_2_sd + theta_2_mu;
# }
# model {
#     //hyperpriors
#     theta_1_mu ~ std_normal();
#     theta_2_mu ~ std_normal();
#     theta_1_sd ~ std_normal();
#     theta_2_sd ~ std_normal();
#     
#     //priors
#     theta_1_raw ~ std_normal();
#     theta_2_raw ~ std_normal();
#     
#     //likelihood
#     x1 ~ binomial_logit(n, theta_1);
#     x2 ~ binomial_logit(n, theta_2);
# 
# }
# generated quantities {
#     vector[p] theta_diff = theta_1 - theta_2;
# }
# '
# 
# base_2 = paste0("difference_between_groups_estimated")
# stan_program_2 <- '
# data {
#     int<lower=1> n;
#     int<lower=1> p;
#     int<lower=0,upper=n> x1[p];
#     int<lower=0,upper=n> x2[p];
# }
# parameters {
#     vector[p] theta_mu_raw;
#     vector[p] theta_diff_raw;
#     
#     real theta_mu_mu;
#     real theta_diff_mu;
#     real<lower=0> theta_mu_sd;
#     real<lower=0> theta_diff_sd;
# }
# transformed parameters {
#     vector[p] theta_mu = theta_mu_raw * theta_mu_sd + theta_mu_mu;
#     vector[p] theta_diff = theta_diff_raw * theta_diff_sd + theta_diff_mu;
# 
#     vector[p] theta_1 = theta_mu + theta_diff / 2;
#     vector[p] theta_2 = theta_mu - theta_diff / 2;
# }
# model {
#     //hyperpriors
#     theta_mu_mu ~ std_normal();
#     theta_diff_mu ~ std_normal();
#     theta_mu_sd ~ std_normal();
#     theta_diff_sd ~ std_normal();
#     
#     //priors
#     theta_mu_raw ~ std_normal();
#     theta_diff_raw ~ std_normal();
#     
#     //likelihood
#     x1 ~ binomial_logit(n, theta_1);
#     x2 ~ binomial_logit(n, theta_2);
# }
# '

base_1 = paste0("separate_group_thetas_estimated")
stan_program_1 <- '
data {
    int<lower=1> n;
    int<lower=1> p;
    int<lower=0,upper=n> x1[p];
    int<lower=0,upper=n> x2[p];
}
parameters {
    vector[p] theta_1;
    vector[p] theta_2;
    
    real theta_1_mu;
    real theta_2_mu;
    real<lower=0> theta_1_sd;
    real<lower=0> theta_2_sd;
}
model {
    //hyperpriors
    theta_1_mu ~ std_normal();
    theta_2_mu ~ std_normal();
    theta_1_sd ~ std_normal();
    theta_2_sd ~ std_normal();
    
    //priors
    theta_1 ~ normal(theta_1_mu, theta_1_sd);
    theta_2 ~ normal(theta_2_mu, theta_2_sd);
    
    //likelihood
    x1 ~ binomial_logit(n, theta_1);
    x2 ~ binomial_logit(n, theta_2);

}
generated quantities {
    vector[p] theta_diff = theta_1 - theta_2;
}
'

base_2 = paste0("difference_between_groups_estimated")
stan_program_2 <- '
data {
    int<lower=1> n;
    int<lower=1> p;
    int<lower=0,upper=n> x1[p];
    int<lower=0,upper=n> x2[p];
}
parameters {
    vector[p] theta_mu;
    vector[p] theta_diff;
    
    real theta_mu_mu;
    real theta_diff_mu;
    real<lower=0> theta_mu_sd;
    real<lower=0> theta_diff_sd;
}
transformed parameters {
    vector[p] theta_1 = theta_mu + theta_diff / 2;
    vector[p] theta_2 = theta_mu - theta_diff / 2;
}
model {
    //hyperpriors
    theta_mu_mu ~ std_normal();
    theta_diff_mu ~ std_normal();
    theta_mu_sd ~ std_normal();
    theta_diff_sd ~ std_normal();
    
    //priors
    theta_mu ~ normal(theta_mu_mu, theta_mu_sd);
    theta_diff ~ normal(theta_diff_mu, theta_diff_sd);
    
    //likelihood
    x1 ~ binomial_logit(n, theta_1);
    x2 ~ binomial_logit(n, theta_2);
}
'



if(!exists("curr_stan_program_1") || stan_program_1 != curr_stan_program_1){
  curr_stan_program_1 <- stan_program_1
  f_1 <- write_stan_file(stan_program_1)
}
mod_1 <- cmdstan_model(f_1)

if(!exists("curr_stan_program_2") || stan_program_2 != curr_stan_program_2){
  curr_stan_program_2 <- stan_program_2
  f_2 <- write_stan_file(stan_program_2)
}
mod_2 <- cmdstan_model(f_2)

#write model
write_stan_file(stan_program_1, dir = "~/Desktop/", basename = base_1)
write_stan_json(d_1, paste0("~/Desktop/", paste0(base_1,".json")))
write_stan_file(stan_program_2, dir = "~/Desktop/", basename = base_2)
write_stan_json(d_2, paste0("~/Desktop/", paste0(base_2,".json")))

#fit model
out_1 <- mod_1$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, 
                      adapt_delta = 0.95, refresh = 50, init = 0.1, max_treedepth = 15, thin = 2)
check_mcmc_diagnostics <- F
if(check_mcmc_diagnostics){
  summ_1 <- out_1$summary()
  summ_1[order(summ_1$ess_bulk),]
  summ_1[order(summ_1$rhat, decreasing = T),]
}

out_2 <- mod_2$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, 
                      adapt_delta = 0.95, refresh = 50, init = 0.1, max_treedepth = 15, thin = 2)
check_mcmc_diagnostics <- F
if(check_mcmc_diagnostics){
  summ_2 <- out_2$summary()
  summ_2[order(summ_2$ess_bulk),]
  summ_2[order(summ_2$rhat, decreasing = T),]
}

#extract samples
samps_1 <- data.frame(as_draws_df(out_1$draws()))
samps_2 <- data.frame(as_draws_df(out_2$draws()))

#examine output
par(mfrow = c(2,4))
plot(1,1,xlim = c(0,1), ylim = c(0,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
text(x = 0.5, y = 0.5 , labels = "model where we infer each group's theta\nand find the difference between them", xpd = NA)
theta_diffs_1 <- subset_samps("theta_diff", c("mu", "sd"), samps_1)
hist(sapply(1:p, function(i) mean(theta_diffs_1[,i] > theta_diff[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Difference between Theta_1 and Theta_2")
theta_1s_1 <- subset_samps("theta_1", c("mu", "sd"), samps_1)
hist(sapply(1:p, function(i) mean(theta_1s_1[,i] > theta_1[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Theta_1")
theta_2s_1 <- subset_samps("theta_2", c("mu", "sd"), samps_1)
hist(sapply(1:p, function(i) mean(theta_2s_1[,i] > theta_2[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Theta_2")

plot(1,1,xlim = c(0,1), ylim = c(0,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
text(x = 0.5, y = 0.5 , labels = "model where we parameterize the difference\nbetween theta_1 and theta_2 directly", xpd = NA)
theta_diffs_2 <- subset_samps("theta_diff", c("mu", "sd"), samps_2)
hist(sapply(1:p, function(i) mean(theta_diffs_2[,i] > theta_diff[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Difference between Theta_1 and Theta_2")
theta_1s_2 <- subset_samps("theta_1", c("mu", "sd"), samps_2)
hist(sapply(1:p, function(i) mean(theta_1s_2[,i] > theta_1[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Theta_1")
theta_2s_2 <- subset_samps("theta_2", c("mu", "sd"), samps_2)
hist(sapply(1:p, function(i) mean(theta_2s_2[,i] > theta_2[i])), breaks = 0:20/20,
     xlab = "quantile of the true value of the parameter\nin the marginal posterior distribution", probability = T,
     main = "Theta_2")
