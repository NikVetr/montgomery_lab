n <- 10000
p <- 3
x <- matrix(rnorm(p*n), nrow = n, ncol = p)
b <- rnorm(p)
y <- c(x %*% t(t(b)))
prop_var_x <- 0.5
x_var_sum <- sum(b^2)
y <- y / sqrt(x_var_sum) * sqrt(prop_var_x / (1-prop_var_x)) + rnorm(n, sd = 1)
y <- (y - mean(y)) / sd(y)
sdx <- apply(x, 2, sd)
x <- x %*% diag(1/sdx)

# est_b <- coefficients(lm(y ~ x))[-1]
est_b <- sapply(1:p, function(i) coefficients(lm(y ~ x[,i]))[-1])
est_ss <- sum(est_b^2)
est_ss

pchisq(est_ss, df = p / (n-3))

