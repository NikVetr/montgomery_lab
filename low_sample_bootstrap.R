n <- 1000

foo <- function(n, CI = 0.95){
  x <- rnorm(n)
  sum(quantile(replicate(1E3, mean(sample(x, n, T))), c((1-CI)/2,(1+CI)/2)) > 0) == 1
}

mean(replicate(5E2, foo(n)))
