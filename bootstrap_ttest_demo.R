n_sims <- 1000
n_obs <- 50
n_bootstrap_rep <- 1000

out <- data.frame(t(replicate(n_sims, {
  
  x1 <- rnorm(n_obs)
  x2 <- rnorm(n_obs)

  bootstrapped_means <- t(replicate(n_bootstrap_rep, {
    x1b <- sample(x1, n_obs, T)
    x2b <- sample(x2, n_obs, T)
    c(mean(x1b), mean(x2b))
  }))
  
  ttest_pvalue <- t.test(x = bootstrapped_means[,1], y = bootstrapped_means[,2])$p.value
  
  n_diff <- 5E3
  bootstrapped_diffs <- sample(bootstrapped_means[,1], n_diff, T) - sample(bootstrapped_means[,2], n_diff, T)
  diff_pvalue <- 1 - abs(0.5 - mean(bootstrapped_diffs > 0)) * 2
  
  return(c(ttest = ttest_pvalue, diff = diff_pvalue))
})))

par(mfrow = c(2,1))
hist(out$ttest, main = "bootstrapped t.test p-values")
hist(out$diff, main = "bootstrap difference p-value")
