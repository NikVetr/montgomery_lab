n <- 1000
x <- rnorm(n)
a <- 1
b <- 2
sigma <- 2
y <- a + b*x + rnorm(n, sd = sigma)

fit <- lm(y ~ 1 + x)
sfit <- summary(fit)

sfit$sigma^2
var(fit$residuals)




k <- 10
seq(1, n, length.out = k+1)
