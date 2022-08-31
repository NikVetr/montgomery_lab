#####

n <- 1E3
x <- rnorm(n, 0, 0.1)
y <- rbinom(n, 1, 0.5)
a <- 1
b <- 2
c <- 3
d <- 0
z <- exp(a + b * x + c * y + d * x * y + rnorm(n, sd = 0.5))

fit1 <- lm(log(z) ~ 1 + x + y + x * y)
fit2 <- lm(z ~ 1 + x + y + x * y)

par(mfrow = c(2, 1))

plot(x,
     z,
     pch = 19,
     col = adjustcolor(y + 1, 0.2),
     main = "red and black lines are not parallel,\ninteraction term is estimated to not be 0")
abline(
  a = fit2$coefficients[1],
  fit2$coefficients[2],
  lwd = 2,
  lty = 2
)
abline(
  a = sum(fit2$coefficients[c(1, 3)]),
  sum(fit2$coefficients[c(2, 4)]),
  lwd = 2,
  lty = 2,
  col = 2
)

plot(x,
     log(z),
     pch = 19,
     col = adjustcolor(y + 1, 0.2),
     main = "red and black lines are parallel,\ninteraction term is estimated to be 0")
abline(
  a = fit1$coefficients[1],
  fit1$coefficients[2],
  lwd = 2,
  lty = 2
)
abline(
  a = sum(fit1$coefficients[c(1, 3)]),
  sum(fit1$coefficients[c(2, 4)]),
  lwd = 2,
  lty = 2,
  col = 2
)


#####

foo <- function() {
  n <- 1E2
  x <- rnorm(n, 0, 0.1)
  y <- rbinom(n, 1, 0.5)
  a <- 1
  b <- 2
  c <- 3
  d <- 0
  z <- exp(a + b * x + c * y + d * x * y + rnorm(n, sd = 0.5))
  
  fit1 <- lm(log(z) ~ 1 + x + y + x * y)
  fit2 <- lm(z ~ 1 + x + y + x * y)
  
  c(fit1 = summary(fit1)$coefficients["x:y", "Pr(>|t|)"],
    fit2 = summary(fit2)$coefficients["x:y", "Pr(>|t|)"])
}

dev.off()
pvals <- t(replicate(1000, foo()))
hist(pvals, col = "blue")
hist(pvals[, 1], add = T, col = "orange")
legend(
  x = "topright",
  legend = c("identity link", "log link"),
  col = c("blue", "orange"),
  pch = 15,
  pt.cex = 2
)
