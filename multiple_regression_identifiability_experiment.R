#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)
library(TESS)
library(adephylo)

#overall graphical params
par(mfrow = c(2,3))


#specify a few functions
rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

#simulate data

#tree-like covariance
speciation <- function(t){
  return(5 / (1 + (20)*exp(-0.75*t) ))
}
rbf <- function(x1, x2, params = list(s = 1, l = 1)){
  params$s^2*exp(-(x2-x1)^2/(2*params$l^2))
}
brownian <- function(x1, x2, params = list(s=1)){return(params$s^2*min(x1, x2))}
tree <- tess.sim.taxa(n = 1,
                      nTaxa = n,
                      max = 1E4,
                      lambda = speciation,
                      mu = 0)[[1]]
vcov_tree <- vcv.phylo(tree)
# dist_tree <- as.matrix(distTips(tree))
# (sapply(1:n, function(ri) sapply(1:n, function(ci) brownian(dist_tree[ri], dist_tree[ci]))))
vcov_tree <- vcov_tree / vcov_tree[1,1] * sigma2

#population params
n = 200
p = 50
b_var <- 1
# b <- rexp(p, rate = 1 / sqrt(b_var)) * (1-rbinom(p,1,0.5)*2)
b <- rnorm(p, sd = sqrt(b_var)) 
b[sample(1:p, size = p / 10)] <- rnorm(p / 10, mean = 0, sd = sqrt(b_var * 100))

x <- matrix(rnorm(n*p),n,p)
corr_x <- 0.5
corrmat_x <- diag(n)*(1-corr_x)+corr_x
# corrmat_x <- cov2cor(vcov_tree)
x <- t(chol(corrmat_x)) %*% x
sigma2 <- 25
e <- rnorm(n, sd = 1)
# e <- t(vcov_tree) %*% e
hist(e)
a = rnorm(1, sd = 10)
y <- a + x %*% b + e

#fit ols model
ols_bs_mr_cf <- lm(y ~ 1 + x)$coefficients[-1]

#marginally

#check understanding
cov_x <- cov(x) #true sample cov(x)
cov_x <- vcov_tree #true sample cov(x)
eig <- eigen(cov(x))
v <- eig$vectors
l <- eig$values
mean_center <- function(x) sapply(1:ncol(x), function(ci) x[,ci] - mean(x[,ci]))

#compute PCs
PCs <- mean_center(x) %*% v
n_PCs_to_include <- 3
partial_rank_cov_x <- v[,1:n_PCs_to_include] %*% diag(l[1:n_PCs_to_include]) %*% t(v[,1:n_PCs_to_include])

ols_bs <- sapply(1:p, function(i) lm(y ~ 1 + x[,i] + PCs[,1:n_PCs_to_include])$coefficients[2])
ols_bs_noPCs <- sapply(1:p, function(i) lm(y ~ 1 + x[,i])$coefficients[2])
cov_x_y <- ols_bs * apply(x, 2, var)
(cov2cor(cov_x) - cov2cor(partial_rank_cov_x))
cov_x_to_use <- diag(sqrt(diag(cov_x))) %*% (cov2cor(cov_x) - cov2cor(partial_rank_cov_x) + diag(p)) %*% diag(sqrt(diag(cov_x)))
cov_x_to_use <- cov_x

ols_bs_mr <- c(pracma::pinv(cov_x_to_use) %*% cov_x_y)
# plot(b, ols_bs, ylim = range(c(ols_bs, ols_bs_mr)));abline(0,1)
# plot(b, ols_bs_mr, ylim = range(c(ols_bs, ols_bs_mr)));abline(0,1)
plot(ols_bs_mr_cf, ols_bs_mr, ylim = range(c(ols_bs, ols_bs_mr)));abline(0,1)

ols_bs_noPCs_true_mr <- lm(y ~ x)$coefficients[-1]
ols_bs_noPCs_mr <- c(pracma::pinv(cov_x) %*% (ols_bs_noPCs * apply(x, 2, var)))
plot(b, ols_bs); abline(0,1)
plot(b, ols_bs_noPCs); abline(0,1)
plot(b, ols_bs_noPCs_mr); abline(0,1)
plot(b, ols_bs_noPCs_true_mr); abline(0,1)


# cor(b, ols_bs)
# cor(b, ols_bs_mr)
# 
# cov_x <- cov(x)
# cov_xy <- cbind(rbind(cov_x, cov_x_y), c(cov_x_y, var(y)))
# r2 <- t(cov2cor(cov_xy)[p+1,-(p+1)]) %*% pracma::pinv(cov2cor(cov_x)) %*% t(t(cov2cor(cov_xy)[p+1,-(p+1)]))
# adj_r2 <- 1-((1- r2) * (n-1) / (n-p-1))
# cov_coef <- pracma::pinv(cov_x) / c(((n-1) / (var(y) * (1-adj_r2))))
# coef_SE <- sqrt(diag(cov_coef))
# plot(coef_SE, summary(ols_bs_mr)$coefficients[-1,2]); abline(0,1,col=2)
# 
# #fit bayesian model
# d <- list(n = n,
#          p = p,
#          x = x,
#          y = c(y))
# 
# stan_program <- "
# data {
#   int<lower=0> n;
#   int<lower=0> p;
#   matrix[n, p] x;
#   vector[n] y;
# }
# parameters {
#   vector[p] b;
#   real<lower=0> b_var;
#   real<lower=0> sigma2;
#   real a;
# }
# model {
#   a ~ normal(0,10);
#   sigma2 ~ exponential(0.1);
#   b_var ~ exponential(1);
#   // b ~ double_exponential(0, sqrt(b_var));
#   b ~ normal(0, sqrt(1000));
#   
#   y ~ normal(a + x * b, sqrt(sigma2));
#   
# }
# "
# 
# if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
#   curr_stan_program <- stan_program
#   f <- write_stan_file(stan_program)  
# }
# mod <- cmdstan_model(f)
# 
# #fit model
# out <- mod$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
# samps <- data.frame(as_draws_df(out$draws()))
# out
# 
# mean(samps$a)
# a
# 
# mean(samps$sigma2)
# sigma2
# 
# mean(samps$b_var)
# b_var
# 
# par(mfrow = c(3,1))
# plot(b, apply(samps[,grep("b\\.", colnames(samps))],2,mean),
#      main = round(cor(b, apply(samps[,grep("b\\.", colnames(samps))],2,mean)), 3)) 
# abline(0,1,col=2)
# plot(b, ols_bs,main = round(cor(b, ols_bs), 3)) 
# abline(0,1,col=2)
# 
# plot(ols_bs, apply(samps[,grep("b\\.", colnames(samps))],2,mean), 
#      main = round(cor(apply(samps[,grep("b\\.", colnames(samps))],2,mean), ols_bs), 3)) 
# abline(0,1,col=2)
# 
# 
# plot(abs(ols_bs), abs(ols_bs - apply(samps[,grep("b\\.", colnames(samps))],2,mean)))
# 
# par(mfrow = c(1,2))
# plot(sqrt(diag(pracma::pinv(cov_x))), apply(samps[,grep("b\\.", colnames(samps))],2,sd))
# plot(cov2cor(pracma::pinv(cov_x))[upper.tri(diag(p))], cor(samps[,grep("b\\.", colnames(samps))])[upper.tri(diag(p))]);abline(0,1)
