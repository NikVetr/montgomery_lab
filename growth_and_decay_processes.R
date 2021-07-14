library(rethinking); set_ulam_cmdstan(F)
library(cmdstanr)
library(rstan)

#### specify simulation model ####

set.seed(777)
P <- function(K, r, t, L0){ #logistic function
  return(K / (1 + ((K-L0)/L0)*exp(-r*t) ))
}
t0 = 0; t1 = 64 #start and end time of the growth process
t <- t0:(t1*4)/4
r1 = 0.25 #growth rate
K1 <- 2500 #asymptote of growth process
L01 <- 1 #starting state of growth process
rates1 <- P(K1, r1, t, L01)
r2 = -.5 #decay rate
K2 <- 1000
L02 <- 1
L02 <- K2 - L02 #starting state of overall process
rates2 <- P(K2, r2, t, L02)
rates <- rates1 + rates2
plot(t, rates, type = "l", ylim = range(c(rates, rates1, rates2)), lwd = 2)
lines(t, rates1, col = "green")
text(x = t1, y = tail(rates1, 1) - 15, labels = "growth process", col = "green", pos = 2, srt = 1)
lines(t, rates2, col = "red")
text(x = t1, y = tail(rates2, 1) + 15, labels = "decay process", col = "red", pos = 2, srt = - 5)
lines(t, rates, col = "black")
text(x = t1, y = tail(rates, 1) + 15, labels = "overall process", col = "black", pos = 2, srt = - 5)

#### simulate data ####

n_obs <- 128 #number of observations to draw
prop_to_sample <- 1
# sampleTimes <- sort(runif(round(n_obs * prop_to_sample), range(t)[1], range(t)[2] * prop_to_sample)) #uniformly distributed sampling times
sampleTimes <- cumsum(extraDistr::rdirichlet(n = 1, alpha = rep(5, round(n_obs * prop_to_sample))) * diff(range(t)) * prop_to_sample) #dirichlet distributed sampling times
# sampleTimes <- c(1/11, 2/11, 4/11, 8/11) * 64
sampleExpRates <- P(K1, r1, sampleTimes, L01) + P(K2, r2, sampleTimes, L02) #sampled expected rates
theta <- 3 #scale parameter of gamma distribution
sampleRates <- rgamma(round(n_obs * prop_to_sample), shape = sampleExpRates / theta, scale = theta) #sampled rates at each timepoint
sampleObs <- rpois(round(n_obs * prop_to_sample), sampleRates) #poisson-distributed counts at each timepoint
d <- list(time = sampleTimes, count = sampleObs) #put data into list
plot(t, rates, type = "l", ylim = c(0,5000))
points(sampleTimes, sampleObs)

#### specify stan model ####

double_logistic_model_sc <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK1;
    real<lower=0> r1;
    real<lower=0, upper = 1> L01_prop;
    real logK2;
    real<lower=0> neg_r2;
    real<lower=0, upper = 1> L02_prop;
    real<lower=0> theta;
}
transformed parameters{
    real<lower = 0> K1 = exp(logK1);
    real<lower = 0> K2 = exp(logK2);
    real<lower = 0> L01 = K1 * L01_prop;
    real<lower = 0> L02 = K2 * L02_prop;
    real<upper = 0> r2 = -neg_r2;
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta ~ exponential( 1 );
    
    //growth model
    r1 ~ exponential( 1 );
    L01_prop ~ beta(1, 100);
    logK1 ~ normal( log(10) , 3 * log(10) );
    
    //decay model
    neg_r2 ~ exponential( 1 );
    L02_prop ~ beta(100, 1);
    logK2 ~ normal( log(10) , 3 * log(10) );
    
    //mathemagically combine them
    for ( i in 1:128 ) {
        lambda[i] = K1 / (1 + ((K1 - L01) / L01) * exp(-r1 * time[i])) + 
                    K2 / (1 + ((K2 - L02) / L02) * exp(-r2 * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    
    //evaluate likelihood
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"

double_logistic_model_sc_flatpriors <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK1;
    real<lower=0> r1;
    real<lower=0, upper = 1> L01_prop;
    real logK2;
    real<lower=0> neg_r2;
    real<lower=0, upper = 1> L02_prop;
    real<lower=0> theta;
}
transformed parameters{
    real<lower = 0> K1 = exp(logK1);
    real<lower = 0> K2 = exp(logK2);
    real<lower = 0> L01 = K1 * L01_prop;
    real<lower = 0> L02 = K2 * L02_prop;
    real<upper = 0> r2 = -neg_r2;
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta ~ exponential( 0.01 );
    
    //growth model
    r1 ~ exponential( 0.01 );
    L01_prop ~ beta(1, 1);
    logK1 ~ normal( log(10) , 50 * log(10) );
    
    //decay model
    neg_r2 ~ exponential( 0.01 );
    L02_prop ~ beta(1, 1);
    logK2 ~ normal( log(10) , 50 * log(10) );
    
    //mathemagically combine them
    for ( i in 1:128 ) {
        lambda[i] = K1 / (1 + ((K1 - L01) / L01) * exp(-r1 * time[i])) + 
                    K2 / (1 + ((K2 - L02) / L02) * exp(-r2 * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    
    //evaluate likelihood
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"


double_logistic_model_noncentered_sc <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK1_std;
    real<lower=0> r1_std;
    real<lower=0, upper = 1> L01_prop;
    real logK2_std;
    real<lower=0> neg_r2_std;
    real<lower=0, upper = 1> L02_prop;
    real<lower=0> theta_std;
}
transformed parameters{
    real<lower = 0> K1 = exp(log(10) + logK1_std * 3 * log(10));
    real<lower = 0> K2 = exp(log(10) + logK2_std * 3 * log(10));
    real<lower = 0> L01 = K1 * L01_prop;
    real<lower = 0> L02 = K2 * L02_prop;
    real<lower = 0> r1 = r1_std / 1;
    real<upper = 0> r2 = -neg_r2_std / 1;
    real<lower=0> theta = theta_std / 0.5;
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta_std ~ exponential( 1 );
    
    //growth model
    r1_std ~ exponential( 1 );
    L01_prop ~ beta(1, 100);
    logK1_std ~ std_normal();
    
    //decay model
    neg_r2_std ~ exponential( 1 );
    L02_prop ~ beta(100, 1);
    logK2_std ~ std_normal();
    
    //mathemagically combine them
    for ( i in 1:128 ) {
        lambda[i] = K1 / (1 + ((K1 - L01) / L01) * exp(-r1 * time[i])) + 
                    K2 / (1 + ((K2 - L02) / L02) * exp(-r2 * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    
    //evaluate likelihood
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"

double_logistic_model_noncentered_sc_flatpriors <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK1_std;
    real<lower=0> r1_std;
    real<lower=0, upper = 1> L01_prop;
    real logK2_std;
    real<lower=0> neg_r2_std;
    real<lower=0, upper = 1> L02_prop;
    real<lower=0> theta_std;
}
transformed parameters{
    real<lower = 0> K1 = exp(log(10) + logK1_std * 10 * log(10));
    real<lower = 0> K2 = exp(log(10) + logK2_std * 10 * log(10));
    real<lower = 0> L01 = K1 * L01_prop;
    real<lower = 0> L02 = K2 * L02_prop;
    real<lower = 0> r1 = r1_std / 0.01;
    real<upper = 0> r2 = -neg_r2_std / 0.01;
    real<lower=0> theta = theta_std / 0.01;
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta_std ~ exponential( 1 );
    
    //growth model
    r1_std ~ exponential( 1 );
    L01_prop ~ beta(1, 1);
    logK1_std ~ std_normal();
    
    //decay model
    neg_r2_std ~ exponential( 1 );
    L02_prop ~ beta(1, 1);
    logK2_std ~ std_normal();
    
    //mathemagically combine them
    for ( i in 1:128 ) {
        lambda[i] = K1 / (1 + ((K1 - L01) / L01) * exp(-r1 * time[i])) + 
                    K2 / (1 + ((K2 - L02) / L02) * exp(-r2 * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    
    //evaluate likelihood
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"

#### fit stan models #### 

double_logistic_model_sc <- gsub(x = double_logistic_model_sc, pattern = "128", replacement = round(n_obs * prop_to_sample))
double_logistic_model_noncentered_sc <- gsub(x = double_logistic_model_noncentered_sc, pattern = "128", replacement = round(n_obs * prop_to_sample))
n_chains <- 4
initf_nonc <- function(chain_id = 1) {list(L01_prop = 0.1, r1_std = 0.01, theta_std = 1, logK1_std = 1, 
                                      L02_prop = 0.9, neg_r2_std = 0.01, logK1_std = 1, alpha = chain_id)}
initf_cent <- function(chain_id = 1) {list(L01_prop = 0.1, r1 = 1, theta1 = 1, logK1 = 5, 
                                      L02_prop = 0.9, r2 = 1, theta2 = 1, logK2 = 5, alpha = chain_id)}

init_ll <- lapply(1:n_chains, function(id) initf_cent(chain_id = id))
double_logistic_model_sc_flatpriors <- gsub(x = double_logistic_model_sc_flatpriors, pattern = "128", replacement = n_obs)
double_logistic_model_noncentered_sc_flatpriors <- gsub(x = double_logistic_model_noncentered_sc_flatpriors, pattern = "128", replacement = n_obs)

fit_regpriors <- stan(model_code = double_logistic_model_sc, data = d, chains = n_chains, warmup = 3E4, 
                      iter = 4E4, cores = 4, control = list(adapt_delta = 0.995), thin = 3)
fit_flatpriors <- stan(model_code = double_logistic_model_sc_flatpriors, data = d, chains = n_chains, warmup = 1.5E4, 
                       iter = 2E4, cores = 4, init = init_ll, control = list(adapt_delta = 0.999, max_treedepth = 11), thin = 3)


#### visualize fitted model results ####

#par is c(K1, L01_prop, r1, K2, L02_prop, r2, theta)
gampois_logisticgrowth_lincomb_nopenlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  # nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(r1, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  # nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(-r2, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  # nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
varlppm_colrange <- viridis::magma(100, end = 0.8)
# for(fit in c(fit_regpriors, fit_flatpriors)){
for(fit in c(fit_regpriors)){
  
  logism <- extract.samples(fit)
  thin <- 5
  logism_expected_values <- sapply(1:(length(logism$K1) / thin), function(iter) 
    P(K = logism$K1[iter*thin], r = logism$r1[iter*thin], t = t, L0 = logism$L01[iter*thin]) +
    P(K = logism$K2[iter*thin], r = logism$r2[iter*thin], t = t, L0 = logism$L02[iter*thin]))
  logism_89HPDI <- apply(logism_expected_values, 1, HPDI)
  
  logism_sample_rates <- sapply(1:(length(logism$logK1) / thin), function(iter) 
    rgamma(n = length(t), shape = logism_expected_values[,iter] / logism$theta[iter*thin], scale = logism$theta[iter*thin]))
  logism_89HPDI_rates <- apply(logism_sample_rates, 1, HPDI)
  
  logism_sample_obs <- do.call(cbind, lapply(1:3, function(x) sapply(1:(length(logism$logK1) / thin), function(iter) 
    rpois(n = length(t), lambda = logism_sample_rates[,iter]))))
  logism_89HPDI_obs <- apply(logism_sample_obs, 1, HPDI)
  
  
  plot(t, rates, type = "l", ylim = c(0,max(c(logism_89HPDI_obs[2,],rates))), lwd = 2, col = "darkred")
  # plot(t, rates, type = "l", ylim = c(0,max(rates)), lwd = 2, col = "darkred")
  
  col_lambdaHPDI <- "#00416c"
  col_sampHPDI <- "#c24e00"
  
  lines(t, logism_89HPDI[1,], col = col_lambdaHPDI)
  lines(t, logism_89HPDI[2,], col = col_lambdaHPDI)
  polygon(y = c(logism_89HPDI[2,], rev(logism_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)
  
  # lines(t, logism_89HPDI_rates[1,], col = col_lambdaHPDI)
  # lines(t, logism_89HPDI_rates[2,], col = col_lambdaHPDI)
  # polygon(y = c(logism_89HPDI_rates[2,], rev(logism_89HPDI_rates[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5))
  
  lines(t, logism_89HPDI_obs[1,], col = col_sampHPDI)
  lines(t, logism_89HPDI_obs[2,], col = col_sampHPDI)
  polygon(y = c(logism_89HPDI_obs[2,], rev(logism_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
  polygon(y = c(logism_89HPDI[1,], rev(logism_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
  
  # color points by std.var(lppm)
  thin <- 100
  vlppms <- sapply(1:round(n_obs*prop_to_sample), function(obs) var(sapply(1:(length(logism$K1) / thin), function(iter) 
    gampois_logisticgrowth_lincomb_nopenlik(data = as.data.frame(d)[obs,], 
                                            par = c(logism$K1[iter*thin], logism$L01_prop[iter*thin], logism$r1[iter*thin], logism$K2[iter*thin], 
                                                    logism$L02_prop[iter*thin], logism$r2[iter*thin], logism$theta[iter*thin])))))
  vlppms <- vlppms / max(vlppms)
  points(sampleTimes, sampleObs, pch = 19, col = grDevices::adjustcolor(varlppm_colrange[ceiling(vlppms*100)], alpha = 0.9))
  
}


#### now use optimx to find the MLE #####

data <- as.data.frame(d)
#par is c(K1, L01_prop, r1, K2, L02_prop, r2, theta)
invlogit <- function(x){return(1 / (1+exp(-x)))}
gampois_logisticgrowth_lincomb <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}

gampois_logisticgrowth_lincomb_nopenlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  # nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(r1, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  # nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(-r2, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  # nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}

gampois_logisticgrowth_lincomb_reparam <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
#par is c(K1, L01_prop, r1, K2, L02_prop, r2, theta)
gampois_logisticgrowth_lincomb(data, par = c(2500, 0.001 / 2.5, 0.25, 1000, 0.999, -0.5, 8))
gampois_logisticgrowth_lincomb_reparam(data, par = c(0, 0, 0, 0, 0, 0, 0))
# est_all <- optimx::optimx(par = c(2000, 0.1, 1, 1000, 0.9, -1, 1), fn = gampois_logisticgrowth_lincomb, data = data, 
#                           method = "nlmimb", lower = c(1E-6, 1E-6, 1E-6, 1E-6, 1E-6, -Inf, 1E-6), 
#                           upper = c(Inf, 1-1E-6, Inf, Inf, 1-1E-6, -1E-6, Inf), control = list(maxit = 1E5, trace = 2))
# est_all <- optimx::optimx(par = c(1000, 0.1, 1, 1000, 0.9, -1, 1), fn = gampois_logisticgrowth_lincomb, data = data, 
#                           method = c('Nelder-Mead','nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rvmmin'), control = list(maxit = 1E5, trace = 2))
est_all <- optimx::optimx(par = c(5,-5,-2,5,5,-2,2), fn = gampois_logisticgrowth_lincomb_reparam, data = data, 
                          method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                          control = list(maxit = 1E4, trace = 0))
est <- est_all[which.min(est_all$value),]

# est <- optim(par = c(1000, 0.01, 1, 1), fn = gampois_logisticgrowth, data = data,
#              method = "Nelder-Mead", control = list(trace = 1))

est <- est_transf <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
gampois_logisticgrowth_lincomb_reparam(data, par = unlist(est))
gampois_logisticgrowth_lincomb(data, par = c(2500, 0.001 / 2.5, 0.25, 1000, 0.999, -0.5, 8))
names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
est$K1 <- exp(est$K1)
est$L01_prop <- invlogit(est$L01_prop)
est$r1 <- exp(est$r1)
est$K2 <- exp(est$K2)
est$L02_prop <- invlogit(est$L02_prop)
est$r2 <- -exp(est$r2)
est$theta <- exp(est$theta)
est$L01 <- est$K1 * est$L01_prop
est$L02 <- est$K2 * est$L02_prop
gampois_logisticgrowth_lincomb(data, par = unname(unlist(est[1:7])))
est_L <- P(K = est$K1, r = est$r1, t = t, L0 = est$L01) + P(K = est$K2, r = est$r2, t = t, L0 = est$L02)

plot(t, rates, type = "l", ylim = c(0,5000), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
lines(t, est_L, col = "#009ECE", lwd = 2)
lines(t, rates, col = "red", lwd = 2)
legend(x = "topleft", col=c("red", "#009ECE"), legend = c("true model", "estimated model"), lwd = 2)


#### find LOO-DeltaMass in MLE framework ####

gampois_logisticgrowth_lincomb_reparam_penlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_reparam_nopenlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_nopenlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  # nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(r1, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  # nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(-r2, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  # nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_penlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
plls <- sapply(1:n_obs, function(obs) gampois_logisticgrowth_lincomb_nopenlik(data = as.data.frame(d)[obs,], 
                                      par = c(est$K1, est$L01_prop, est$r1, est$K2, est$L02_prop, est$r2, est$theta)))
loo_plls <- rep(NA, n_obs)
est_transf_alldata <- unlist(est_transf)
for(obs in 1:n_obs){
  cat(paste0("(", obs, ") " ))
  
  lik <- "penlik"
  #fit MLE
  est_all <- optimx::optimx(par = est_transf_alldata, 
                            fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                            data = as.data.frame(d)[-obs,], method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                            control = list(maxit = 1E4, trace = 0))
  est <- est_all[which.min(est_all$value),]
  est <- est_transf <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
  
  if(lik == "penlik"){
    ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(as.data.frame(d)[-obs,], par = unlist(est))
    ll_true <- gampois_logisticgrowth_lincomb_penlik(as.data.frame(d)[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
  } else if(lik == "nopenlik"){
    ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(as.data.frame(d)[-obs,], par = unlist(est))
    ll_true <- gampois_logisticgrowth_lincomb_nopenlik(as.data.frame(d)[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
  }
  
  if(ll_est > ll_true){
    retry_max <- 1000
    retry_i <- 0
    while(retry_i < retry_max & ll_est > ll_true){
      par_init <- rnorm(n = length(est_transf_alldata), mean = est_transf_alldata, sd = 0.5 * retry_i / 100)
      est_all <- optimx::optimx(par = par_init, 
                                fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                data = as.data.frame(d)[-obs,], method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                control = list(maxit = 5E4, trace = 0))
      est <- est_all[which.min(est_all$value),]
      est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
      #compare ll of MLE to true params
      if(lik == "penlik"){
        ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(as.data.frame(d)[-obs,], par = unlist(est))
        ll_true <- gampois_logisticgrowth_lincomb_penlik(as.data.frame(d)[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      } else if(lik == "nopenlik"){
        ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(as.data.frame(d)[-obs,], par = unlist(est))
        ll_true <- gampois_logisticgrowth_lincomb_nopenlik(as.data.frame(d)[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      }
      retry_i <- retry_i + 1
    }
  }
  
  names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
  est$K1 <- exp(est$K1)
  est$L01_prop <- invlogit(est$L01_prop)
  est$r1 <- exp(est$r1)
  est$K2 <- exp(est$K2)
  est$L02_prop <- invlogit(est$L02_prop)
  est$r2 <- -exp(est$r2)
  est$theta <- exp(est$theta)
  est$L01 <- est$K1 * est$L01_prop
  est$L02 <- est$K2 * est$L02_prop
  
  loo_plls[obs] <- gampois_logisticgrowth_lincomb_nopenlik(data = as.data.frame(d)[obs,], par = c(est$K1, est$L01_prop, est$r1, 
                                                                                 est$K2, est$L02_prop, est$r2, est$theta))
  
}

delta_plss <- loo_plls - plls
delta_plss <- delta_plss / max(delta_plss)

dplss_colrange <- viridis::magma(100, end = 0.8)
plot(t, rates, type = "l", ylim = c(0,2695), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs, col = grDevices::adjustcolor(dplss_colrange[ceiling(delta_plss*100)], alpha = 0.9), pch = 19)
lines(t, est_L, col = "#009ECE", lwd = 2)
lines(t, rates, col = "red", lwd = 2)
legend(x = "topleft", col=c("red", "#009ECE"), legend = c("true model", "estimated model"), lwd = 2)


#### do nonparametric bootstrapping ####

data <- as.data.frame(d)
invlogit <- function(x){return(1 / (1+exp(-x)))}
gampois_logisticgrowth_lincomb_reparam_penlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_reparam_nopenlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_penlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_nopenlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  # nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(r1, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  # nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(-r2, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  # nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}

library(foreach)
library(parallel)
library(doParallel)
do_bootstrap <- T
n_bootstrap_replicates <- 2E3
if(do_bootstrap){
  
  if(!exists("cl")){
    cl <- makeCluster(8, outfile="")
    registerDoParallel(cl)
  }
  
  getDoParWorkers()
  
  foreach(i=76:1, .packages = c("optimx")) %dopar% {
    
    for(lik in c("nopenlik", "penlik")){
      
      fileout <- paste0("~/Documents/logistic_processes/growth_and_decay/mle_bootstrap_", lik, "_", i, "_obs")
      
      #get data subset
      dsub <- d
      dsub$time <- dsub$time[1:i]
      dsub$count <- dsub$count[1:i]
      data <- as.data.frame(dsub)
      
      #resample data w/ replacement
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) sample(x = 1:i, size = i, replace = T))
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) data[data_bootstrap[[bsr]],])
      
      #initialize output object
      fit <- list(MLE = NA, BSREPS = matrix(NA, nrow = n_bootstrap_replicates, ncol= 9), 
                  ll_true = NA, ll_est = NA, ll_true_bsreps = rep(NA, n_bootstrap_replicates), 
                  ll_est_bsreps = rep(NA, n_bootstrap_replicates), plls = rep(NA, n_obs), loo_plls = rep(NA, n_obs))
      
      #fit MLE
      est_all <- optimx::optimx(par = c(5,-5,-2,5,5,-2,2), 
                                fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                data = data, method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                control = list(maxit = 1E4, trace = 0))
      est <- est_all[which.min(est_all$value),]
      est <- est_transf_alldata <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
      est_transf_alldata <- unlist(est_transf_alldata)
      
      #compare ll of MLE to true params
      if(lik == "penlik"){
        fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data, par = unlist(est))
        fit$ll_true <- gampois_logisticgrowth_lincomb_penlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      } else if(lik == "nopenlik"){
        fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data, par = unlist(est))
        fit$ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      }
      
      if(fit$ll_est > fit$ll_true){
        retry_max <- 1000
        retry_i <- 0
        while(retry_i < retry_max & fit$ll_est > fit$ll_true){
          par_init <- rnorm(n = length(c(5,-5,-2,5,5,-2,2)), mean = c(5,-5,-2,5,5,-2,2), sd = 0.5 * retry_i / 100)
          est_all <- optimx::optimx(par = par_init, 
                                    fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                    data = data, method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                    control = list(maxit = 5E4, trace = 0))
          est <- est_all[which.min(est_all$value),]
          est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
          #compare ll of MLE to true params
          if(lik == "penlik"){
            fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data, par = unlist(est))
            fit$ll_true <- gampois_logisticgrowth_lincomb_penlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
          } else if(lik == "nopenlik"){
            fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data, par = unlist(est))
            fit$ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
          }
          retry_i <- retry_i + 1
        }
      }
      
      names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
      est$K1 <- exp(est$K1)
      est$L01_prop <- invlogit(est$L01_prop)
      est$r1 <- exp(est$r1)
      est$K2 <- exp(est$K2)
      est$L02_prop <- invlogit(est$L02_prop)
      est$r2 <- -exp(est$r2)
      est$theta <- exp(est$theta)
      est$L01 <- est$K1 * est$L01_prop
      est$L02 <- est$K2 * est$L02_prop
      fit$MLE <- unlist(est)
      
      #find MLE for all bootstrap replicates
      for(j in 1:n_bootstrap_replicates){
        cat(paste0("(", i, ", ", j, ") " ))
        est_all <- optimx::optimx(par = est_transf_alldata, fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), data = data_bootstrap[[j]], 
                                  method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                  control = list(maxit = 1E4, trace = 0))
        est <- est_all[which.min(est_all$value),]
        est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))

        if(lik == "penlik"){
          fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_penlik(data_bootstrap[[j]], par = unlist(est))
          fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_penlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        } else if(lik == "nopenlik"){
          fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data_bootstrap[[j]], par = unlist(est))
          fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_nopenlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        }
        
        if(fit$ll_est_bsreps[j] > fit$ll_true_bsreps[j]){
          retry_max <- 100
          retry_i <- 0
          while(retry_i < retry_max & fit$ll_est_bsreps[j] > fit$ll_true_bsreps[j]){
            par_init <- rnorm(n = length(est_transf_alldata), mean = est_transf_alldata, sd = 0.5 * retry_i / 10)
            
            est_all <- optimx::optimx(par = par_init, 
                                      fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), data = data_bootstrap[[j]], 
                                      method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                      control = list(maxit = 1E4, trace = 0))
            est <- est_all[which.min(est_all$value),]
            est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
            
            #compare ll of MLE to true params
            if(lik == "penlik"){
              fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_penlik(data_bootstrap[[j]], par = unlist(est))
              fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_penlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            } else if(lik == "nopenlik"){
              fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data_bootstrap[[j]], par = unlist(est))
              fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_nopenlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            }
            retry_i <- retry_i + 1
          }
        }
        
        names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
        est$K1 <- exp(est$K1)
        est$L01_prop <- invlogit(est$L01_prop)
        est$r1 <- exp(est$r1)
        est$K2 <- exp(est$K2)
        est$L02_prop <- invlogit(est$L02_prop)
        est$r2 <- -exp(est$r2)
        est$theta <- exp(est$theta)
        est$L01 <- est$K1 * est$L01_prop
        est$L02 <- est$K2 * est$L02_prop
        fit$BSREPS[j,] <- unlist(est)
        
      }
      
      #find pointwise ll for data, leaving one out and not
      fit$plls <- sapply(1:i, function(obs) gampois_logisticgrowth_lincomb_reparam_nopenlik(data = data[obs,], par = est_transf_alldata))
      for(obs in 1:i){
        cat(paste0("(", obs, ") " ))
        #fit MLE
        est_all <- optimx::optimx(par = est_transf_alldata, 
                                  fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                  data = data[-obs,], method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                  control = list(maxit = 1E4, trace = 0))
        est <- est_all[which.min(est_all$value),]
        est <- est_transf <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
        
        if(lik == "penlik"){
          ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data[-obs,], par = unlist(est))
          ll_true <- gampois_logisticgrowth_lincomb_penlik(data[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        } else if(lik == "nopenlik"){
          ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data[-obs,], par = unlist(est))
          ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        }
        
        if(ll_est > ll_true){
          retry_max <- 1000
          retry_i <- 0
          while(retry_i < retry_max & ll_est > ll_true){
            par_init <- rnorm(n = length(est_transf_alldata), mean = est_transf_alldata, sd = 0.5 * retry_i / 100)
            est_all <- optimx::optimx(par = par_init, 
                                      fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                      data = as.data.frame(d)[-obs,], method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                      control = list(maxit = 5E4, trace = 0))
            est <- est_all[which.min(est_all$value),]
            est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
            #compare ll of MLE to true params
            if(lik == "penlik"){
              ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data[-obs,], par = unlist(est))
              ll_true <- gampois_logisticgrowth_lincomb_penlik(data[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            } else if(lik == "nopenlik"){
              ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data[-obs,], par = unlist(est))
              ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data[-obs,], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            }
            retry_i <- retry_i + 1
          }
        }
        
        names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
        est$K1 <- exp(est$K1)
        est$L01_prop <- invlogit(est$L01_prop)
        est$r1 <- exp(est$r1)
        est$K2 <- exp(est$K2)
        est$L02_prop <- invlogit(est$L02_prop)
        est$r2 <- -exp(est$r2)
        est$theta <- exp(est$theta)
        est$L01 <- est$K1 * est$L01_prop
        est$L02 <- est$K2 * est$L02_prop
        
        fit$loo_plls[obs] <- gampois_logisticgrowth_lincomb_nopenlik(data = data[obs,], par = c(est$K1, est$L01_prop, est$r1, 
                                                                                                        est$K2, est$L02_prop, est$r2, est$theta))
        
      }
      
      
      #save bootstrap object
      save(fit, file = fileout)
      rm(fit)
      
    }
    
    
  }
  
  stopCluster(cl)
  rm(cl)
  
}

#### do parametric bootstrapping ####

data <- as.data.frame(d)
invlogit <- function(x){return(1 / (1+exp(-x)))}
gampois_logisticgrowth_lincomb_reparam_penlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_reparam_nopenlik <- function(data, par) {
  
  #undo params w/ transformation
  K1 <- exp(par[1])
  L01_prop <- invlogit(par[2])
  r1 <- exp(par[3])
  K2 <- exp(par[4])
  L02_prop <- invlogit(par[5])
  r2 <- -exp(par[6])
  theta <- exp(par[7])
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  # if(is.na(nloglik)){
  #   cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
  #   cat(paste0("\nr1 = ", r1, "\n"))
  #   cat(paste0("\nK1 = ", K1, "\n"))
  #   cat(paste0("\nL01 = ", L01, "\n"))
  #   cat(paste0("\nr2 = ", r2, "\n"))
  #   cat(paste0("\nK2 = ", K2, "\n"))
  #   cat(paste0("\nL02 = ", L02, "\n"))
  #   cat(paste0("\ntheta = ", theta, "\n"))
  #   return(Inf)
  # }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_penlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(r1, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(-r2, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_lincomb_nopenlik <- function(data, par) {
  
  #undo params
  K1 <- par[1]
  L01_prop <- par[2]
  r1 <- par[3]
  K2 <- par[4]
  L02_prop <- par[5]
  r2 <- par[6]
  theta <- par[7]
  
  #impose bounds on unbounded algos muahaha
  if(L01_prop < 0 | L01_prop > 1 | K1 < 0 | theta < 0 | L02_prop < 0 | L02_prop > 1 | K2 < 0 | r1 < 0 | r2 > 0){return(Inf)}
  
  # model
  L01 = K1 * L01_prop
  L02 = K2 * L02_prop
  lambda1 = K1 /(1 + ((K1 - L01)/L01) * exp(-r1 * data$time))
  lambda2 = K2 /(1 + ((K2 - L02)/L02) * exp(-r2 * data$time))
  lambda <- lambda1 + lambda2
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood = ", nloglik, "\n"))
    cat(paste0("\nr1 = ", r1, "\n"))
    cat(paste0("\nK1 = ", K1, "\n"))
    cat(paste0("\nL01 = ", L01, "\n"))
    cat(paste0("\nr2 = ", r2, "\n"))
    cat(paste0("\nK2 = ", K2, "\n"))
    cat(paste0("\nL02 = ", L02, "\n"))
    cat(paste0("\ntheta = ", theta, "\n"))
  }
  
  #penalized likelihood?
  # nloglik <- nloglik + -dnorm(log10(K1), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(r1, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L01_prop, shape1 = 1, shape2 = 100, log = T)
  # nloglik <- nloglik + -dnorm(log10(K2), mean = 1, sd = 3, log = T)
  # nloglik <- nloglik + -dexp(-r2, 2, log = T)
  # nloglik <- nloglik + -dbeta(x = L02_prop, shape1 = 100, shape2 = 1, log = T)
  # nloglik <- nloglik + -dexp(theta, 1, log = T)
  
  return(nloglik)
  
}

library(foreach)
library(parallel)
library(doParallel)
do_bootstrap <- T
n_bootstrap_replicates <- 2000
if(do_bootstrap){
  
  if(!exists("cl")){
    cl <- makeCluster(8, outfile="")
    registerDoParallel(cl)
  }
  
  getDoParWorkers()
  
  foreach(i=76:1, .packages = c("optimx")) %dopar% {
    
    for(lik in c("nopenlik", "penlik")){
      
      fileout <- paste0("~/Documents/logistic_processes/growth_and_decay/mle_parametric_bootstrap_", lik, "_", i, "_obs")
      
      #get data subset
      dsub <- d
      dsub$time <- dsub$time[1:i]
      dsub$count <- dsub$count[1:i]
      data <- as.data.frame(dsub)
      
      #initialize output object
      fit <- list(MLE = NA, BSREPS = matrix(NA, nrow = n_bootstrap_replicates, ncol= 9), 
                  ll_true = NA, ll_est = NA, ll_true_bsreps = rep(NA, n_bootstrap_replicates), ll_est_bsreps = rep(NA, n_bootstrap_replicates))
      
      #fit MLE
      est_all <- optimx::optimx(par = c(5,-5,-2,5,5,-2,2), 
                                fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                data = data, method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                control = list(maxit = 5E4, trace = 0))
      est <- est_all[which.min(est_all$value),]
      est <- est_transf_alldata <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
      est_transf_alldata <- unlist(est_transf_alldata)
      
      #compare ll of MLE to true params
      if(lik == "penlik"){
        fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data, par = unlist(est))
        fit$ll_true <- gampois_logisticgrowth_lincomb_penlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      } else if(lik == "nopenlik"){
        fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data, par = unlist(est))
        fit$ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
      }
      
      if(fit$ll_est > fit$ll_true){
        retry_max <- 1000
        retry_i <- 0
        while(retry_i < retry_max & fit$ll_est > fit$ll_true){
          par_init <- rnorm(n = length(c(5,-5,-2,5,5,-2,2)), mean = c(5,-5,-2,5,5,-2,2), sd = 0.5 * retry_i / 100)
          est_all <- optimx::optimx(par = par_init, 
                                    fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), 
                                    data = data, method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                    control = list(maxit = 5E4, trace = 0))
          est <- est_all[which.min(est_all$value),]
          est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
          #compare ll of MLE to true params
          if(lik == "penlik"){
            fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_penlik(data, par = unlist(est))
            fit$ll_true <- gampois_logisticgrowth_lincomb_penlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
          } else if(lik == "nopenlik"){
            fit$ll_est <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data, par = unlist(est))
            fit$ll_true <- gampois_logisticgrowth_lincomb_nopenlik(data, par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
          }
          retry_i <- retry_i + 1
        }
      }
      
      names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
      est$K1 <- exp(est$K1)
      est$L01_prop <- invlogit(est$L01_prop)
      est$r1 <- exp(est$r1)
      est$K2 <- exp(est$K2)
      est$L02_prop <- invlogit(est$L02_prop)
      est$r2 <- -exp(est$r2)
      est$theta <- exp(est$theta)
      est$L01 <- est$K1 * est$L01_prop
      est$L02 <- est$K2 * est$L02_prop
      fit$MLE <- unlist(est)
      
      #resimulate data from fitted model
      sampleExpRates_bootstrap <- P(est$K1, est$r1, data$time, est$L01) + P(est$K2, est$r2, data$time, est$L02)
      sampleRates_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) rgamma(i, shape = sampleExpRates_bootstrap / est$theta , scale = est$theta))
      sampleObs_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) rpois(i, sampleRates_bootstrap[[bsr]]))
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) data.frame(time = data$time, count = sampleObs_bootstrap[[bsr]]))
      
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) sample(x = 1:i, size = i, replace = T))
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) data[data_bootstrap[[bsr]],])
      
      #find MLE for all bootstrap replicates
      for(j in 1:n_bootstrap_replicates){
        cat(paste0("(", i, ", ", j, ") " ))
        est_all <- optimx::optimx(par = est_transf_alldata, 
                                  fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), data = data_bootstrap[[j]], 
                                  method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                  control = list(maxit = 1E4, trace = 0))
        est <- est_all[which.min(est_all$value),]
        est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))
        
        if(lik == "penlik"){
          fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_penlik(data_bootstrap[[j]], par = unlist(est))
          fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_penlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        } else if(lik == "nopenlik"){
          fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data_bootstrap[[j]], par = unlist(est))
          fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_nopenlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
        }
        
        if(fit$ll_est_bsreps[j] > fit$ll_true_bsreps[j]){
          retry_max <- 100
          retry_i <- 0
          while(retry_i < retry_max & fit$ll_est_bsreps[j] > fit$ll_true_bsreps[j]){
            par_init <- rnorm(n = length(est_transf_alldata), mean = est_transf_alldata, sd = 0.5 * retry_i / 10)
            
            est_all <- optimx::optimx(par = par_init, 
                                      fn = ifelse(lik == "penlik", gampois_logisticgrowth_lincomb_reparam_penlik, gampois_logisticgrowth_lincomb_reparam_nopenlik), data = data_bootstrap[[j]], 
                                      method = c('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'nlm', 'nlminb'), 
                                      control = list(maxit = 1E4, trace = 0))
            est <- est_all[which.min(est_all$value),]
            est <- as.list(as.numeric(est[names(est)[grep(names(est), pattern = "p")]]))

            #compare ll of MLE to true params
            if(lik == "penlik"){
              fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_penlik(data_bootstrap[[j]], par = unlist(est))
              fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_penlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            } else if(lik == "nopenlik"){
              fit$ll_est_bsreps[j] <- gampois_logisticgrowth_lincomb_reparam_nopenlik(data_bootstrap[[j]], par = unlist(est))
              fit$ll_true_bsreps[j] <- gampois_logisticgrowth_lincomb_nopenlik(data_bootstrap[[j]], par = c(K1, L01 / K1, r1, K2, L02 / K2, r2, theta))
            }
            retry_i <- retry_i + 1
          }
        }
        
        names(est) <- c("K1", "L01_prop", "r1", "K2", "L02_prop", "r2", "theta")
        est$K1 <- exp(est$K1)
        est$L01_prop <- invlogit(est$L01_prop)
        est$r1 <- exp(est$r1)
        est$K2 <- exp(est$K2)
        est$L02_prop <- invlogit(est$L02_prop)
        est$r2 <- -exp(est$r2)
        est$theta <- exp(est$theta)
        est$L01 <- est$K1 * est$L01_prop
        est$L02 <- est$K2 * est$L02_prop
        fit$BSREPS[j,] <- unlist(est)
        
      }
      
      #save bootstrap object
      save(fit, file = fileout)
      rm(fit)
      
    }
    
    
  }
  
  stopCluster(cl)
  rm(cl)
  
}


#### fit all replicate models from 2 obs to 128 ####

#first the regularizing priors
do_regpriors <- T
use_noncentered_model = T
n_chains <- 4
pre_multiplier <- 8
initf_nonc <- function(chain_id = 1) {list(L01_prop = 0.1, r1_std = 0.01, theta_std = 1, logK1_std = 1, 
                                           L02_prop = 0.9, neg_r2_std = 0.01, logK1_std = 1, alpha = chain_id)}
initf_cent <- function(chain_id = 1) {list(L01_prop = 0.1, r1 = 1, theta1 = 1, logK1 = 5, 
                                           L02_prop = 0.9, r2 = 1, theta2 = 1, logK2 = 5, alpha = chain_id)}
if(use_noncentered_model){init_ll <- lapply(1:n_chains, function(id) initf_nonc(chain_id = id))}else{init_ll <- lapply(1:n_chains, function(id) initf_cent(chain_id = id))}

if(do_regpriors){
  
  for(i in 10:5){
    
    fileout <- paste0("~/Documents/logistic_processes/growth_and_decay/regpriors_", i, "_obs")
    print(fileout)
    
    #check if previous fit exists and load if so, to avoid recompiling
    if(file.exists(fileout)){
      load(fileout)
    }
    
    #get data subset
    dsub <- d
    dsub$time <- dsub$time[1:i]
    dsub$count <- dsub$count[1:i]
    
    #get model code
    double_logistic_model_noncentered_sc_sub <- gsub(x = double_logistic_model_noncentered_sc, pattern = n_obs, replacement = i)
    double_logistic_model_sc_sub <- gsub(x = double_logistic_model_sc, pattern = n_obs, replacement = i)
    
    #specify fit requirements
    n_acceptable_divergent_transitions <- 1
    min_req_ess <- 2E3
    max_acceptable_rhat <- 1.01
    min_ess <- 0
    n_divergent_transitions <- 100
    max_rhat <- 10
    n_warmup <- 0.75E4 * pre_multiplier
    n_iter <- 1E4 * pre_multiplier
    adapt_delta = 0.9995
    max_td = 15
    thin = pre_multiplier
    while(min_ess < min_req_ess | n_divergent_transitions > n_acceptable_divergent_transitions | max_rhat > max_acceptable_rhat){
      
      cat(paste0(i, " "))
      
      #fit model
      if(!exists("fit")){
        if(use_noncentered_model){
          print("using noncentered model")
          fit <- stan(model_code = double_logistic_model_noncentered_sc_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
        } else {
          print("using centered model")
          fit <- stan(model_code = double_logistic_model_sc_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                      iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
        }
      } else{
        fit <- stan(fit = fit, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
      }
      
      #do automated diagnostics
      sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
      n_divergent_transitions <- max(sapply(sampler_params, function(x) max(x[, "divergent__"])))
      max_rhat <- max(summary(fit)$summary[, "Rhat"])
      min_ess <- min(summary(fit)$summary[, "n_eff"])
      
      # #do more intensive diagnostics in case of just one or two failed chains for cases where 
      # if(i < 30 & n_chains > 4){
      #   s <- extract.samples(fit)
      #   samp_per_chain <- length(s$logK1_std) / n_chains
      #   sep_chains <- sapply(1:n_chains, function(chain) lapply(s, function(param))
      #   
      # }
      
      #make necessary adjustments
      if(min_ess < min_req_ess | max_rhat > max_acceptable_rhat){
        n_warmup <- n_warmup * 2
        n_iter <- n_iter * 2
        thin <- thin * 2
        if(n_iter > 200000){min_req_ess <- 1000; max_acceptable_rhat <- 1.01}
        print(paste0("ess = ", min_ess))
        print(paste0("rhat = ", max_rhat))
      } 
      if(n_divergent_transitions > n_acceptable_divergent_transitions){
        adapt_delta <- mean(c(adapt_delta, 1))
        if(adapt_delta > 0.999){n_acceptable_divergent_transitions <- n_acceptable_divergent_transitions +1}
        print(paste0("# divergent transitions: ", n_divergent_transitions))
      } 
      
    }
    
    #save stanfit object
    save(fit, file = fileout)
    rm(fit)
    
  }
  
}


do_flatpriors <- T
use_noncentered_model = T
n_chains <- 4
initf_nonc <- function(chain_id = 1) {list(L01_prop = 0.1, r1_std = 0.01, theta_std = 1, logK1_std = 1, 
                                           L02_prop = 0.9, neg_r2_std = 0.01, logK1_std = 1, alpha = chain_id)}
initf_cent <- function(chain_id = 1) {list(L01_prop = 0.1, r1 = 1, theta1 = 1, logK1 = 5, 
                                           L02_prop = 0.9, r2 = 1, theta2 = 1, logK2 = 5, alpha = chain_id)}
if(use_noncentered_model){init_ll <- lapply(1:n_chains, function(id) initf_nonc(chain_id = id))}else{init_ll <- lapply(1:n_chains, function(id) initf_cent(chain_id = id))}

if(do_flatpriors){
  
  for(i in 128:2){
    
    pre_multiplier <- ceiling( (129-i) / 10)
    
    fileout <- paste0("~/Documents/logistic_processes/growth_and_decay/flatpriors_", i, "_obs")
    print(fileout)
    
    #check if previous fit exists and load if so, to avoid recompiling
    # if(file.exists(fileout)){
    #   load(fileout)
    # }
    
    #get data subset
    dsub <- d
    dsub$time <- dsub$time[1:i]
    dsub$count <- dsub$count[1:i]
    
    #get model code
    double_logistic_model_noncentered_sc_flatpriors_sub <- gsub(x = double_logistic_model_noncentered_sc_flatpriors, pattern = n_obs, replacement = i)
    double_logistic_model_sc_flatpriors_sub <- gsub(x = double_logistic_model_sc_flatpriors, pattern = n_obs, replacement = i)
    
    #specify fit requirements
    n_acceptable_divergent_transitions <- 1
    min_req_ess <- 2E3
    max_acceptable_rhat <- 1.01
    min_ess <- 0
    n_divergent_transitions <- 100
    max_rhat <- 10
    n_warmup <- 0.75E4 * pre_multiplier
    n_iter <- 1E4 * pre_multiplier
    adapt_delta = 0.9 + min(c(0.099, pre_multiplier / 100 - 0.01))
    max_td = 15
    thin = pre_multiplier
    while(min_ess < min_req_ess | n_divergent_transitions > n_acceptable_divergent_transitions | max_rhat > max_acceptable_rhat){
      
      cat(paste0(i, " "))
      
      #fit model
      if(!exists("fit")){
        if(use_noncentered_model){
          print("using noncentered model")
          fit <- stan(model_code = double_logistic_model_noncentered_sc_flatpriors_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                      iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
        } else {
          print("using centered model")
          fit <- stan(model_code = double_logistic_model_sc_flatpriors_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                      iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
        }
      } else{
        fit <- stan(fit = fit, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = n_chains, thin = thin, control = list(adapt_delta = adapt_delta, max_treedepth = max_td), init = init_ll)
      }
      
      #do automated diagnostics
      sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
      n_divergent_transitions <- max(sapply(sampler_params, function(x) max(x[, "divergent__"])))
      max_rhat <- max(summary(fit)$summary[, "Rhat"])
      min_ess <- min(summary(fit)$summary[, "n_eff"])
      
      #do more intensive diagnostics in case of just one or two failed chains for cases where
      if(i < 100 & n_chains > 4){
        all_chains <- extract.samples(fit)
        samp_per_chain <- length(all_chains$logK1) / n_chains
        sep_chains <- sapply(1:n_chains, function(chain) lapply(all_chains, function(param) param[(1+(chain-1)*samp_per_chain):(chain*samp_per_chain)]))
      }
      
      #make necessary adjustments
      if(min_ess < min_req_ess | max_rhat > max_acceptable_rhat){
        n_warmup <- n_warmup * 4
        n_iter <- n_iter * 4
        thin <- thin * 4
        if(n_iter > 200000){min_req_ess <- 1000; max_acceptable_rhat <- 1.01}
        print(paste0("ess = ", min_ess))
        print(paste0("rhat = ", max_rhat))
      } 
      if(n_divergent_transitions > n_acceptable_divergent_transitions){
        adapt_delta <- mean(c(adapt_delta, 1))
        if(adapt_delta > 0.999){n_acceptable_divergent_transitions <- n_acceptable_divergent_transitions +1}
        print(paste0("# divergent transitions: ", n_divergent_transitions))
      } 
      
    }
    
    #save stanfit object
    save(fit, file = fileout)
    rm(fit)
    
  }
  
}


#### show samples from the prior ####

library(rethinking)
P <- function(K, r, t, L0){ #logistic function
  return(K / (1 + ((K-L0)/L0)*exp(-r*t) ))
}

#regularizing prior
t0 = 0; t1 = 64 #start and end time of the growth process
t <- t0:(t1*4)/4
nsamp <- 10000
K1s <- 10^rnorm(nsamp, mean = 1, sd = 3)
r1s <- rexp(nsamp, 1)
L01_props <- rbeta(nsamp, shape1 = 1, shape2 = 100)
L01s <- K1s * L01_props
K2s <- 10^rnorm(nsamp, mean = 1, sd = 3)
r2s <- -rexp(nsamp, 1)
L02_props <- rbeta(nsamp, shape1 = 100, shape2 = 1)
L02s <- K2s * L02_props
thetas <- rexp(nsamp, 1)
ratess <- t(sapply(1:nsamp, function(x) P(K = K1s[x], r = r1s[x], t = t, L0 = L01s[x]) + P(K = K2s[x], r = r2s[x], t = t, L0 = L02s[x])))
prior_89HPDI <- apply(ratess, 2, HPDI)
prior_sample_rates <- sapply(1:nsamp, function(iter) 
  rgamma(n = length(t), shape = ratess[iter,] / thetas[iter], scale = thetas[iter]))
prior_sample_obs <- sapply(1:nsamp, function(iter) rpois(n = length(t), lambda = prior_sample_rates[,iter]))
prior_89HPDI_obs <- apply(prior_sample_obs, 1, HPDI)

plot(t, rates, type = "l", ylim = c(0,max(prior_89HPDI_obs)), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
nlines_toplot <- 1E3
for(i in 1:nlines_toplot){lines(t, ratess[i * trunc(nsamp/nlines_toplot),], col = grDevices::adjustcolor("black", 0.35), lwd = 1)}

col_lambdaHPDI <- "#00416c"
col_sampHPDI <- "#c24e00"
lines(t, prior_89HPDI[1,], col = col_lambdaHPDI)
lines(t, prior_89HPDI[2,], col = col_lambdaHPDI)
polygon(y = c(prior_89HPDI[2,], rev(prior_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)
lines(t, prior_89HPDI_obs[1,], col = col_sampHPDI)
lines(t, prior_89HPDI_obs[2,], col = col_sampHPDI)
polygon(y = c(prior_89HPDI_obs[2,], rev(prior_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
polygon(y = c(prior_89HPDI[1,], rev(prior_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)

#flat prior
t0 = 0; t1 = 64 #start and end time of the growth process
t <- t0:(t1*4)/4
nsamp <- 10000
K1s <- 10^rnorm(nsamp, mean = 1, sd = 50)
r1s <- rexp(nsamp, 0.01)
L01_props <- rbeta(nsamp, shape1 = 1, shape2 = 1)
L01s <- K1s * L01_props
K2s <- 10^rnorm(nsamp, mean = 1, sd = 50)
r2s <- -rexp(nsamp, 0.01)
L02_props <- rbeta(nsamp, shape1 = 1, shape2 = 1)
L02s <- K2s * L02_props
thetas <- rexp(nsamp, 0.01)
ratess <- t(sapply(1:nsamp, function(x) P(K = K1s[x], r = r1s[x], t = t, L0 = L01s[x]) + P(K = K2s[x], r = r2s[x], t = t, L0 = L02s[x])))
prior_89HPDI <- apply(ratess, 2, HPDI)
prior_sample_rates <- sapply(1:nsamp, function(iter) 
  rgamma(n = length(t), shape = ratess[iter,] / thetas[iter], scale = thetas[iter]))
prior_sample_obs <- sapply(1:nsamp, function(iter) rpois(n = length(t), lambda = prior_sample_rates[,iter]))
prior_89HPDI_obs <- apply(prior_sample_obs, 1, HPDI)

plot(t, rates, type = "l", ylim = c(0,max(prior_89HPDI_obs)), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
nlines_toplot <- 1E3
for(i in 1:nlines_toplot){lines(t, ratess[i * trunc(nsamp/nlines_toplot),], col = grDevices::adjustcolor("black", 0.35), lwd = 1)}

col_lambdaHPDI <- "#00416c"
col_sampHPDI <- "#c24e00"
lines(t, prior_89HPDI[1,], col = col_lambdaHPDI)
lines(t, prior_89HPDI[2,], col = col_lambdaHPDI)
polygon(y = c(prior_89HPDI[2,], rev(prior_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)
lines(t, prior_89HPDI_obs[1,], col = col_sampHPDI)
lines(t, prior_89HPDI_obs[2,], col = col_sampHPDI)
polygon(y = c(prior_89HPDI_obs[2,], rev(prior_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
polygon(y = c(prior_89HPDI[1,], rev(prior_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)


#### fit single obs model ####


fileout <- paste0("~/Documents/logistic_processes/growth_and_decay/regpriors_1_obs")
#get data subset
dsub <- as.data.frame(d)[1,]
double_logistic_model_noncentered_sc_flatpriors_singleobs <- "
data{
    int count;
    real time;
}
parameters{
    real logK1_std;
    real<lower=0> r1_std;
    real<lower=0, upper = 1> L01_prop;
    real logK2_std;
    real<lower=0> neg_r2_std;
    real<lower=0, upper = 1> L02_prop;
    real<lower=0> theta_std;
}
transformed parameters{
    real<lower = 0> K1 = exp(log(10) + logK1_std * 3 * log(10));
    real<lower = 0> K2 = exp(log(10) + logK2_std * 3 * log(10));
    real<lower = 0> L01 = K1 * L01_prop;
    real<lower = 0> L02 = K2 * L02_prop;
    real<lower = 0> r1 = r1_std / 1;
    real<upper = 0> r2 = -neg_r2_std / 1;
    real<lower=0> theta = theta_std / 0.5;
}
model{
    real lambda;
    real alpha;
    real beta;
    theta_std ~ exponential( 1 );
    
    //growth model
    r1_std ~ exponential( 1 );
    L01_prop ~ beta(1, 100);
    logK1_std ~ std_normal();
    
    //decay model
    neg_r2_std ~ exponential( 1 );
    L02_prop ~ beta(100, 1);
    logK2_std ~ std_normal();
    
    //mathemagically combine them
    lambda = K1 / (1 + ((K1 - L01) / L01) * exp(-r1 * time)) + 
                K2 / (1 + ((K2 - L02) / L02) * exp(-r2 * time));
    alpha = lambda / theta;
    beta = 1 / theta;
  
    
    //evaluate likelihood
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"
initf_nonc <- function(chain_id = 1) {list(L01_prop = 0.1, r1_std = 0.01, theta_std = 1, logK1_std = 1, 
                                           L02_prop = 0.9, neg_r2_std = 0.01, logK1_std = 1, alpha = chain_id)}
init_ll <- lapply(1:n_chains, function(id) initf_nonc(chain_id = id))
fit <- stan(model_code = double_logistic_model_noncentered_sc_flatpriors_singleobs, data = dsub, chains = 4, warmup = 3E4, 
            iter = 4E4, cores = 4, thin = 1, control = list(adapt_delta = 0.975, max_treedepth = 15), init = init_ll)
save(fit, file = fileout)
rm(fit)