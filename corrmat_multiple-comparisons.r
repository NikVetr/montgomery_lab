library(mvtnorm)

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

d <- 10
n <- 100
FT_SE <- 1 / sqrt(n-3)
zscore_thresh_2t <- abs(qnorm(0.025, 0, 1)) #for 2-tailed test

corrmat <- rlkj(K = d, eta = 1)
x <- rmvnorm(n, sigma = corrmat)
obs_corrmat <- cor(x)

#let's just keep it in symmetric matrix form for clarity's sake
atanh_corrmat <- atanh_obs_corrmat <- diag(d)

atanh_corrmat[upper.tri(atanh_corrmat)] <- atanh(corrmat[upper.tri(corrmat)])
atanh_corrmat <- atanh_corrmat + t(atanh_corrmat)

atanh_obs_corrmat[upper.tri(atanh_obs_corrmat)] <- atanh(obs_corrmat[upper.tri(obs_corrmat)])
atanh_obs_corrmat <- atanh_obs_corrmat + t(atanh_obs_corrmat)

diag(atanh_corrmat) <- diag(atanh_obs_corrmat) <- NA

#compute sample z-scores
abs_zscores <- abs(atanh_corrmat - atanh_obs_corrmat) / FT_SE

#look at false positive rate
sum(abs_zscores > zscore_thresh_2t, na.rm = T) / 2 / choose(d, 2)
#huh look at that

#let's just more explicitly confirm in a linear modeling context
y <- x*0 + matrix(rnorm(d*n), nrow = nrow(x), ncol = ncol(x))
Bs <- sapply(1:d, function(i) (solve(t(cbind(1, x[,i])) %*% cbind(1, x[,i])) %*% t(cbind(1, x[,i])) %*% t(t(y[,i])))[2])
Bs_SEs <- sqrt(sapply(1:d, function(i) var(y[,i]) / var(x[,i]) * (1-cor(x[,i], y[,i])^2) / (n-3)))
Bs_Zs <- abs(Bs / Bs_SEs)
sum(Bs_Zs > zscore_thresh_2t) / d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#let's just do one more quick confirmation replicated at low d
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#doing this cos lkj(eta=1) concentrates probability density around identity which might throw things off

low_d_replicate <- function(d, n, bonf = F){
  
  FT_SE <- 1 / sqrt(n-3)
  zscore_thresh_2t <- abs(qnorm(0.025, 0, 1)) #for 2-tailed test
  zscore_thresh_2t_pairwise_bonf <- abs(qnorm(0.025 / choose(d,2), 0, 1)) #for 2-tailed test
  zscore_thresh_2t_totalD_bonf <- abs(qnorm(0.025 / d, 0, 1)) #for 2-tailed test
  
  
  corrmat <- rlkj(K = d, eta = 1)
  x <- rmvnorm(n, sigma = corrmat)
  obs_corrmat <- cor(x)
  
  #let's just keep it in symmetric matrix form for clarity's sake
  atanh_corrmat <- atanh_obs_corrmat <- diag(d)
  
  atanh_corrmat[upper.tri(atanh_corrmat)] <- atanh(corrmat[upper.tri(corrmat)])
  atanh_corrmat <- atanh_corrmat + t(atanh_corrmat)
  
  atanh_obs_corrmat[upper.tri(atanh_obs_corrmat)] <- atanh(obs_corrmat[upper.tri(obs_corrmat)])
  atanh_obs_corrmat <- atanh_obs_corrmat + t(atanh_obs_corrmat)
  
  diag(atanh_corrmat) <- diag(atanh_obs_corrmat) <- NA
  
  #compute sample z-scores
  abs_zscores <- abs(atanh_corrmat - atanh_obs_corrmat) / FT_SE
  
  #let's just more explicitly confirm in a linear modeling context
  y <- x*0 + matrix(rnorm(d*n), nrow = nrow(x), ncol = ncol(x))
  Bs <- sapply(1:d, function(i) (solve(t(cbind(1, x[,i])) %*% cbind(1, x[,i])) %*% t(cbind(1, x[,i])) %*% t(t(y[,i])))[2])
  Bs_SEs <- sqrt(sapply(1:d, function(i) var(y[,i]) / var(x[,i]) * (1-cor(x[,i], y[,i])^2) / (n-3)))
  Bs_Zs <- abs(Bs / Bs_SEs)
  
  if(bonf){
    return(c(sum(abs_zscores > zscore_thresh_2t_pairwise_bonf, na.rm = T) / 2, sum(Bs_Zs > zscore_thresh_2t_totalD_bonf)))
  } else {
    return(c(sum(abs_zscores > zscore_thresh_2t, na.rm = T) / 2, sum(Bs_Zs > zscore_thresh_2t)))
  }

}

nrep <- 1E4
d = 5
n = 100
apply(replicate(nrep, low_d_replicate(d, n)), 1, sum) / c(choose(d, 2) * nrep, nrep * d)
apply(replicate(nrep, low_d_replicate(d, n, T) > 0), 1, sum) / nrep
#yep yep seems to hold still! okie-dokie