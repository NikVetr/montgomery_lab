#load libraries

#specify functions
#rbf
rbf <- function(x1, x2, params = list(s = 1, l = 1)){
  params$s^2*exp(-(x2-x1)^2/(2*params$l^2))
}
rbf(1, 2, list(s = 1, l = 1))

#brownian
brownian <- function(x1, x2, params = list(s=1)){return(params$s^2*min(x1, x2))}
brownian(2,3, params = list(s=2))


#turn into covariance matrix
cov_mat <- function(x, params, fun){
  n <- length(x)
  matdat <- lapply(2:n, function(ri) sapply(1:(ri-1), function(ci) fun(x[ri], x[ci], params)))
  mat <-  diag(n) - diag(n)
  for(ri in 2:n){
    for(ci in 1:(ri-1)){
      mat[ri,ci] <- matdat[[ri-1]][ci]
    }
  }
  mat <- mat + t(mat)
  if(identical(fun, rbf)){
    diag(mat) <- rep(params$s^2, n)
  }
  if(identical(fun, brownian)){
    diag(mat) <- x*s^2
  }
  return(mat)
}
logit <- function(p) log(p/(1-p))
invlogit <- function(x) exp(x) / (1+exp(x))
  
#visualization params
dim <- 50
x <- seq(0, 15, length.out = dim)

# incorporate observations!
true_param_vals <- list(muy = 4, muz = 5, coef_yz = 0, error_sig2 = 0.05, error_sig2_rbf = 1.5, lengthscale_rbf = 2)
nobs <- 5
obs <- list(x = rnorm(nobs, mean = 5, sd = 3))
obs$y <- mvtnorm::rmvnorm(1, mean = rep(true_param_vals$muy, nobs), sigma = 
                            cov_mat(obs$x, params = list(s = sqrt(true_param_vals$error_sig2_rbf), 
                                                         l = true_param_vals$lengthscale_rbf), fun = rbf) + diag(nobs) * true_param_vals$error_sig2)
obs$z <- mvtnorm::rmvnorm(1, mean = rep(true_param_vals$muz, nobs) + obs$y * true_param_vals$coef_yz, sigma = 
                            cov_mat(obs$x, params = list(s = sqrt(true_param_vals$error_sig2_rbf), 
                                                         l = true_param_vals$lengthscale_rbf), fun = rbf) + diag(nobs) * true_param_vals$error_sig2)

plot(obs$y, obs$z)
summary(lm(c(obs$z) ~ 1 + c(obs$y)))
# obs <- list(x = rnorm(3, mean = 5, sd = 1), y = rnorm(3, 0, 0.5))
# obs$x <- c(2:7); obs$y <- c(5:7, 6:4)
# nobs <- length(obs$x)

#use grid approx to find posterior mean of mu, error_sig2, lengthscale
quantile_spacing_logit_scale <- 1
n_quantiles_to_use <- 7
grid_quantiles <- invlogit(seq(from = -(n_quantiles_to_use-1)/2*quantile_spacing_logit_scale,
                      to = (n_quantiles_to_use-1)/2*quantile_spacing_logit_scale,
                      length.out = n_quantiles_to_use))



#specify mean vector prior
mu_prior_norm_mu_sig <- c(0,3)
mu_prior_dist <- qnorm(p = grid_quantiles, mean = mu_prior_norm_mu_sig[1], sd = mu_prior_norm_mu_sig[2])
mu_prior_priordens <- dnorm(x = mu_prior_dist, mean = mu_prior_norm_mu_sig[1], sd = mu_prior_norm_mu_sig[2])
mu_prior_mean <- sum(mu_prior_dist * mu_prior_priordens) / sum(mu_prior_priordens)
mu_prior_predictive <- rep(mu_prior_mean, dim)

#specify error variance on obs prior
error_sig2_exprate <- 1
error_sig2_dist <- qexp(p = grid_quantiles, rate = error_sig2_exprate)
error_sig2_priordens <- dexp(x = error_sig2_dist, rate = error_sig2_exprate)
error_sig2_prior <- sum(error_sig2_dist * error_sig2_priordens) / sum(error_sig2_priordens)

#specify rbf sd prior
error_sig2_rbf_exprate <- 0.1
error_sig2_rbf_dist <- qexp(p = grid_quantiles, rate = error_sig2_rbf_exprate)
error_sig2_rbf_priordens <- dexp(x = error_sig2_rbf_dist, rate = error_sig2_rbf_exprate)
error_sig2_rbf_prior <- sum(error_sig2_rbf_dist * error_sig2_rbf_priordens) / sum(error_sig2_rbf_priordens)
error_sig_rbf_prior <- sqrt(error_sig2_rbf_prior)

#specify rbf lengthscale prior
lengthscale_rbf_exprate <- 0.2
lengthscale_rbf_dist <- qexp(p = grid_quantiles, rate = lengthscale_rbf_exprate)
lengthscale_rbf_priordens <- dexp(x = lengthscale_rbf_dist, rate = lengthscale_rbf_exprate)
lengthscale_rbf_prior <- sum(lengthscale_rbf_dist * lengthscale_rbf_priordens) / sum(lengthscale_rbf_priordens)

prior_means <- c(mu_prior_mean = mu_prior_mean, error_sig2_prior = error_sig2_prior, 
                 error_sig_rbf_prior = error_sig_rbf_prior, lengthscale_rbf_prior = lengthscale_rbf_prior)

#compute approx posterior of GP params
possible_param_vals <- expand.grid(mu_prior_dist = mu_prior_dist, 
                                   error_sig2_dist = error_sig2_dist, 
                                   error_sig2_rbf_dist = error_sig2_rbf_dist, 
                                   lengthscale_rbf_dist = lengthscale_rbf_dist)
possible_param_vals_prior_dens <- expand.grid(mu_prior_dens = mu_prior_priordens, 
                                              error_sig2_dens = error_sig2_priordens, 
                                              error_sig2_rbf_dens = error_sig2_rbf_priordens, 
                                              lengthscale_rbf_dens = lengthscale_rbf_priordens)
possible_param_vals_joint_prior_dens <- apply(possible_param_vals_prior_dens, 1, prod)

possible_param_vals$likelihood <- sapply(1:nrow(possible_param_vals), function(ri){
  cvmat_grid <- cov_mat(obs$x, params = list(s = sqrt(possible_param_vals$error_sig2_rbf_dist[ri]), 
                               l = possible_param_vals$lengthscale_rbf_dist[ri]), fun = rbf) + diag(nobs) * possible_param_vals$error_sig2_dist[ri]
  cvmat_dens <- mvtnorm::dmvnorm(x = obs$y, mean = rep(possible_param_vals$mu_prior_dist[ri], nobs), sigma = cvmat_grid)
  return(prod(cvmat_dens))
})

possible_param_vals$posterior_dens <- possible_param_vals$likelihood * possible_param_vals_joint_prior_dens

#get posterior means
posterior_means <- sapply(1:(ncol(possible_param_vals)-2), function(param){
  sum(possible_param_vals[,param] * possible_param_vals$posterior_dens) / sum(possible_param_vals$posterior_dens)
})
names(posterior_means) <- stringr::str_remove(colnames(possible_param_vals)[1:(ncol(possible_param_vals)-2)], pattern = "_dist")

error_sig2_posterior <- posterior_means["error_sig2"]
error_sig2_rbf_posterior <- posterior_means["error_sig2_rbf"]
mu_posterior_mean <- posterior_means["mu_prior"]
lengthscale_rbf_posterior <- posterior_means["lengthscale_rbf"]


#compute implied prior predictive
cvmat_prior_predictive <- cov_mat(x, params = list(s = error_sig_rbf_prior, l = lengthscale_rbf_prior), fun = rbf) + diag(dim) * error_sig2_prior
# cvmat <- cov_mat(x, params = list(s = 1), fun = brownian)


#compute posterior predictive
# cvmat_posterior <- cov_mat(x, params = list(s = error_sig2_rbf_posterior, l = lengthscale_rbf_posterior), fun = rbf) + diag(dim) * error_sig2_posterior
cvmat_posterior_incl_obs <- cov_mat(c(x, obs$x), params = list(s = sqrt(error_sig2_rbf_posterior), l = lengthscale_rbf_posterior), fun = rbf) + 
                        diag(dim+nobs) * error_sig2_posterior

cv11 <- cvmat_posterior_incl_obs[1:dim, 1:dim]
cv12 <- cvmat_posterior_incl_obs[1:dim, (dim+1):(dim+nobs)]
cv21 <- t(cv12)
cv22 <- cvmat_posterior_incl_obs[(dim+1):(dim+nobs), (dim+1):(dim+nobs)]
mu_posterior_incl_obs <- rep(mu_posterior_mean, nobs+dim)
mu_posterior_predictive <- mu_posterior_incl_obs[1:dim] + cv12 %*% solve(cv22) %*% t(obs$y - mu_posterior_incl_obs[(dim+1):(dim+nobs)])
cvmat_posterior_predictive <- cv11 - cv12 %*% solve(cv22) %*% cv21


#remove noise from gaussian process model
cvmat_prior_predictive <-  cvmat_prior_predictive - diag(dim) * error_sig2_prior
cvmat_posterior_predictive <-  cvmat_posterior_predictive - diag(dim) * error_sig2_posterior

nrep = 50
y_prior <- mvtnorm::rmvnorm(n = nrep, mean = mu_prior_predictive, sigma = cvmat_prior_predictive)
y_posterior <- mvtnorm::rmvnorm(n = nrep, mean = mu_posterior_predictive, sigma = cvmat_posterior_predictive)


#do some plotting
par(mfrow = c(2,2), mar = c(4,4,1,1))
opacity_power <- 0.5
use_cov_mats <- T
cmat_cols <- c("darkred", "darkblue")

#plot prior predictive cmat
plot(1,1,col="white", xlim = c(0,dim+1), ylim = c(0,dim+1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", frame = F)

if(use_cov_mats){
  for(ri in 1:dim){
    for(ci in 1:dim){
      base_col <- ifelse(sign(cvmat_prior_predictive[ri,ci]) == 1, cmat_cols[1], cmat_cols[2])
      rect(xleft = ri - 1/2, xright = ri + 1/2, ybottom = dim - ci + 1 - 1/2, ytop = dim - ci + 1 + 1/2, border = NA,
             col = adjustcolor(base_col, alpha.f = abs(cvmat_prior_predictive[ri,ci] / max(c(cvmat_prior_predictive, cvmat_posterior_predictive)))^opacity_power))
    }
  }
} else {
  cormat_prior_predictive <- cov2cor(cvmat_prior_predictive)
  for(ri in 1:dim){
    for(ci in 1:dim){
      base_col <- ifelse(sign(cvmat_prior_predictive[ri,ci]) == 1, cmat_cols[1], cmat_cols[2])
      rect(xleft = ri - 1/2, xright = ri + 1/2, ybottom = dim - ci + 1 - 1/2, ytop = dim - ci + 1 + 1/2, border = NA,
           col = adjustcolor(base_col, alpha.f = abs(cormat_prior_predictive[ri,ci])^opacity_power))
    }
  }
}

#plot legend for cmat
par(xpd = NA)
if(use_cov_mats){
  cmat_val_range <- max(c(cvmat_prior_predictive, cvmat_posterior_predictive))
  cmat_val_range <- c(-cmat_val_range, cmat_val_range)
} else {
  cmat_val_range <- c(-1,1)
}

ndiv <- 100
corrvals <- seq(-1,1,length.out = ndiv)
base_cols <- ifelse(sign(corrvals) == 1, cmat_cols[1], cmat_cols[2])
alphas <- abs(corrvals)^opacity_power
xl <- dim * 1.1; xr <- dim * 1.175; yb <- dim * 0.1; yt <- dim * 0.9
rect(xleft = xl, border = NA,
    xright = xr,
    ybottom = seq(yb,yt, length.out = ndiv)[-ndiv],
    ytop = seq(yb,yt, length.out = ndiv)[-1], xpd = NA,
    col = sapply(1:ndiv, function(x) adjustcolor(base_cols[x], alphas[x]))
)
ndiv_text <- 11
corrvals_text <- round(seq(cmat_val_range[1],cmat_val_range[2],length.out = ndiv_text), 2)
text(x = xr * 0.99, y = seq(yb, yt, length.out = ndiv_text), pos = 4, labels = corrvals_text)
text(x = (xl+xr)/2, y = yt, labels = ifelse(use_cov_mats, "Cov.", "Corr."), pos = 3, font = 2, cex = 1.25)

#print prior means
text(labels = latex2exp::TeX(paste0(c("$\\mu_i$", "$\\sigma^2_{error}$", "$\\sigma^2_{RBF}$", "$l_{RBF}$"), 
                                    "$ \\approx $", round(prior_means, 2), collapse = ", ")), x = dim/2, y = 0, pos = 1)

#plot posterior predictive cmat
plot(1,1,col="white", xlim = c(0,dim+1), ylim = c(0,dim+1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", frame = F)
if(use_cov_mats){
  for(ri in 1:dim){
    for(ci in 1:dim){
      base_col <- ifelse(sign(cvmat_prior_predictive[ri,ci]) == 1, "darkred", "darkblue")
      rect(xleft = ri - 1/2, xright = ri + 1/2, ybottom = dim - ci + 1 - 1/2, ytop = dim - ci + 1 + 1/2, border = NA,
           col = adjustcolor(base_col, alpha.f = abs(cvmat_posterior_predictive[ri,ci] / max(c(cvmat_prior_predictive, cvmat_posterior_predictive)))^opacity_power))
    }
  }
} else {
  cormat_posterior_predictive <- cov2cor(cvmat_posterior_predictive)
  for(ri in 1:dim){
    for(ci in 1:dim){
      base_col <- ifelse(sign(cvmat_prior_predictive[ri,ci]) == 1, "darkred", "darkblue")
      rect(xleft = ri - 1/2, xright = ri + 1/2, ybottom = dim - ci + 1 - 1/2, ytop = dim - ci + 1 + 1/2, border = NA,
           col = adjustcolor(base_col, alpha.f = abs(cormat_posterior_predictive[ri,ci])^opacity_power))
    }
  }
}

#print posterior means
text(labels = latex2exp::TeX(paste0(c("$\\mu_i$", "$\\sigma^2_{error}$", "$\\sigma^2_{RBF}$", "$l_{RBF}$"), 
                                    "$ \\approx $", round(posterior_means, 2), collapse = ", ")), x = dim/2, y = 0, pos = 1)


#print true vals
text(labels = "true values", x = -dim/6, y = -6, pos = 3, font = 4, cex = 1.5)
text(labels = latex2exp::TeX(paste0(c("$\\mu_i$", "$\\sigma^2_{error}$", "$\\sigma^2_{RBF}$", "$l_{RBF}$"), 
                                    "$ \\approx $", 
                                    round(as.numeric(true_param_vals[-match(c("muz", "coef_yz"), names(true_param_vals))]), 2), 
                                    collapse = ", ")), x = -dim/6, y = -5, pos = 1)


cols <- RColorBrewer::brewer.pal(8, "Dark2")
# cols <- c(viridisLite::inferno(10), viridisLite::magma(10), viridisLite::cividis(10))
cols <- colorRampPalette(sample(cols, size = 100, replace = T))(1000)
plot(x, y_prior[1,], type = "l", xlim = range(x), ylim = range(c(y_prior, y_posterior)), 
     col = sample(cols, 1), lwd = 2, xlab = "input", ylab = "output", cex.lab = 1.25)
text(x = mean(x), y = max(c(y_prior, y_posterior)) + diff(range(c(y_prior, y_posterior))) / 20, 
     pos = 3, xpd = NA, cex = 2, font = 2, labels = "Prior Predictive")
box("plot", lwd = 2)
for(i in 2:nrep){lines(x, y_prior[i,], col = sample(cols, 1), lwd = 2)}
plot(x, y_posterior[1,], type = "l", xlim = range(x), ylim = range(c(y_prior, y_posterior)), 
     col = sample(cols, 1), lwd = 2, xlab = "input", ylab = "output", cex.lab = 1.25)
text(x = mean(x), y = max(c(y_prior, y_posterior)) + diff(range(c(y_prior, y_posterior))) / 20, 
     pos = 3, xpd = NA, cex = 2, font = 2, labels = "Posterior Predictive")
box("plot", lwd = 2)
for(i in 2:nrep){lines(x, y_posterior[i,], col = sample(cols, 1), lwd = 2)}
points(obs, pch = 19, cex = 1.5, col = adjustcolor("black", 0.5))

print(prior_means)
print(posterior_means)
print(unlist(true_param_vals))
