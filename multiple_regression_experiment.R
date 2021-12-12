#number of parameters, i.e. 'X's
np <- 20

#number of observations
n <- 1E3

#function to sample some correlation matrix
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

#data-generating correlation / covariance matrix
corrmat <- rlkj(np, 0.5)

#sample our 'X's
X <- mvtnorm::rmvnorm(n, rep(0, np), corrmat)

#autocorrelated alleles
X_freq <- replicate(np, rbeta(1,1,2)) #allele frequency across loci
std_normal_thresholds <- qnorm(X_freq)
Xc1 <- mvtnorm::rmvnorm(n, rep(0, np), corrmat)
Xc2 <- mvtnorm::rmvnorm(n, rep(0, np), corrmat)
X <- sapply(1:np, function(gi) as.integer(Xc1[,gi] > std_normal_thresholds[gi]) + as.integer(Xc2[,gi] > std_normal_thresholds[gi]))

#sample our 'B's
B <- rnorm(np, 0, 2)

#sample our intercept
A <- 3

#sample our 'Y's
Y <- A + X%*%t(t(B)) + rnorm(n, sd = 20)

#estimate simple regression coefficients
C_hat <- sapply(1:np, function(p) lm(Y ~ 1 + X[,p])$coefficients[2])

#retrieve covariances between 'X's and 'Y's
CovXY <- C_hat * apply(X, 2, var)

#solve for 'B_hat's using 'C_hat's
CtB_hat <- pracma::pinv(cov(X)) %*% CovXY

#estimate multple regression coefficients
B_hat <- lm(Y~1+X)$coefficients[-1]

#compare
par(mfrow = c(2,1))
plot(B_hat, CtB_hat, pch = 19, col = adjustcolor("blue", 0.5), cex = 1.25); abline(0,1,col=2,lwd=2)
points(B_hat, B, pch = 19, col = adjustcolor("green", 0.5), cex = 1.25); abline(0,1,col=2,lwd=2)

#now obtain multiple regresssion standard errors
# C_hat_SE <- sapply(1:np, function(p) summary(lm(Y ~ 1 + X[,p]))$coefficients[2,2])
B_hat_SE <- summary(lm(Y~1+X))$coefficients[-1,2]
diag(solve(cov(X))) / B_hat_SE^2
var(summary(lm(Y~1+X))$residuals)
var(Y) * (1-summary(lm(Y~1+X))$r.squared)
(n-1) / (var(Y) * (1-summary(lm(Y~1+X))$adj.r.squared))
(n-1) / summary(lm(Y~1+X))$sigma^2

covariance_matrix <- cbind(rbind(cov(X), CovXY), c(CovXY, var(Y)))
r2 <- t(cov2cor(covariance_matrix)[np+1,-(np+1)]) %*% solve(cor(X)) %*% t(t(cov2cor(covariance_matrix)[np+1,-(np+1)]))
summary(lm(Y~1+X))$r.squared
r2
adj_r2 <- 1-((1- r2) * (n-1) / (n-np-1))
summary(lm(Y~1+X))$adj.r.squared
adj_r2

B_hat_SE <- summary(lm(Y~1+X))$coefficients[-1,2]
B_hat_SE
Cov_B <- solve(cov(X)) / c(((n-1) / (var(Y) * (1-adj_r2))))
est_B_hat_SE <- sqrt(diag(Cov_B))
est_B_hat_SE
all(abs(B_hat_SE - est_B_hat_SE) < 1E-6)
plot(B_hat_SE, est_B_hat_SE, pch = 19, col = adjustcolor("blue", 0.5), cex = 1.25); abline(0,1,col=2,lwd=2)

#### let's try adding unreported covariates ####
nZ <- 5
Z <- matrix(rnorm(n*nZ, 0, 1), nrow = n)
Z <- Z# - X[,1:nZ]
B_Z <- rnorm(nZ, 0, 1)
Y_Z <- Y + Z%*%t(t(B_Z))
D_hat <- sapply(1:np, function(p) lm(Y ~ 1 + X[,p] + Z)$coefficients[2])
CovXY_Z <- D_hat * apply(X, 2, var)
CtB_hat_Z <- pracma::pinv(cov(X)) %*% CovXY_Z
B_hat_Z <- lm(Y~1+X+Z)$coefficients[-1][-((np+1):(np+nZ))]
plot(B_hat_Z, CtB_hat_Z, pch = 19, col = adjustcolor(c(rep("orange", nZ), rep("blue", np-nZ)), 0.5), cex = 1.25,
     xlab = "full multiple regression estimate", ylab = "retrieved estimate from partial regressions"); abline(0,1,col=2,lwd=2)
points(B_hat_Z, B, pch = 19, col = adjustcolor("green", 0.5), cex = 1.25); abline(0,1,col=2,lwd=2)
legend("topleft", cex = 0.75, col = adjustcolor(c("blue", "orange", "green"), 0.5), pch = 19, pt.cex = 1.25,
       legend = c("independent from unreported regressors", "covary with unreported regressors", "true values"))
legend("bottomright", cex = 0.75, col = "red", lwd = 2, legend = "1-to-1 line")

covariance_matrix_Z <- cbind(rbind(cov(X), CovXY_Z), c(CovXY_Z, var(Y)))
r2_Z <- t(cov2cor(covariance_matrix_Z)[np+1,-(np+1)]) %*% solve(cor(X)) %*% t(t(cov2cor(covariance_matrix_Z)[np+1,-(np+1)]))
summary(lm(Y~1+X+Z))$r.squared
r2_Z
adj_r2_Z <- 1-((1- r2_Z) * (n-1) / (n-np-1))
summary(lm(Y~1+X+Z))$adj.r.squared
adj_r2_Z

summary(lm(Y_Z~1+Z))$r.squared

B_hat_SE_Z <- summary(lm(Y~1+X+Z))$coefficients[-1,2][-((np+1):(np+nZ))]
B_hat_SE_Z
est_B_hat_SE_Z <- sqrt(diag(solve(cov(X))) / c(((n-1) / (var(Y) * (1-adj_r2_Z)))))
est_B_hat_SE_Z
mean(abs(B_hat_SE_Z - est_B_hat_SE_Z))
plot(B_hat_SE_Z, est_B_hat_SE_Z, pch = 19, col = adjustcolor(c(rep("orange", nZ), rep("blue", np-nZ)), 0.5), cex = 1.25,
     xlab = "full multiple regression estimate", ylab = "retrieved estimate from partial regressions"); abline(0,1,col=2,lwd=2)
legend("topleft", cex = 0.75, col = adjustcolor(c("blue", "orange"), 0.5), pch = 19, pt.cex = 1.25,
       legend = c("independent from unreported regressors", "covary with unreported regressors"))
legend("bottomright", cex = 0.75, col = "red", lwd = 2, legend = "1-to-1 line")
cor(B_hat_SE_Z, est_B_hat_SE_Z)
# plot(B_hat_SE^2, diag(solve(cov(X))))
# n / anova(lm(Y~X))[[3]][2]
# mdev <- function(x) x - mean(x)
# plot(t(apply(X, 2, mdev) ) %*% apply(X, 2, mdev) / (n-1), cov(X))
# solve(cov(X)) %*% B_hat_SE

# points(B_hat, B_that, pch = 19, col = adjustcolor("red", 0.5), cex = 1.25)


# plot(B_hat, B); abline(0,1)
# plot(correlation::cor_to_spcor(cov = covariance_matrix)[np+1,-(np+1)] * sd(Y) / apply(X, 2, sd), lm(Y~1+X)$coefficients[-1]); abline(0,1)
# plot(correlation::cor_to_pcor(cor = covariance_matrix)[np+1,-(np+1)] * sd(Y) / apply(X, 2, sd), lm(Y~1+X)$coefficients[-1]); abline(0,1)
# # plot(correlation::cor_to_spcor(cov = covariance_matrix)[-(np+1),np+1], lm(Y~1+X)$coefficients[-1])

# sum(correlation::cor_to_spcor(cov = covariance_matrix)[np+1,-(np+1)]^2)
# sum(correlation::cor_to_pcor(cor = covariance_matrix)[np+1,-(np+1)]^2) 
