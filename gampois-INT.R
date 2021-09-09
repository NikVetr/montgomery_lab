d_mat <- 10
n_samp <- 100
n_rep <- 1E3
covs <- rWishart(n_rep, n_samp, diag(d_mat)) / n_samp
L <- sapply(1:dim(covs)[3], function(x) eigen(covs[,,x])$values)
hist(L[1,])

n = 1E4
gamma_rates <- rgamma(n,shape = 10, scale = 10)
poisson_counts <- rpois(n, gamma_rates)
INT_vals <- qnorm(rank(poisson_counts) / n, sd = 0.5)
plot(log2(poisson_counts), INT_vals)
abline(-6.5,2, col = 2, lty = 2, lwd = 2)
