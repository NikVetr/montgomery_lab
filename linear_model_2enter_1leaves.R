#load libraries
library(MASS)
library(mvtnorm)

#specify functions
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

b3 = round(rnorm(1,0,2),2)
b2 = round(rnorm(1,0,2),2)
# b2 = NA
sim <- function(b2 = NA, b3 = NA){

  n = 1E3 #n individuals at locus
  
  x_freq <- rbeta(1,1,2) #allele frequency at locus
  x1 <- rbinom(n, 2, prob = x_freq) #genotypes in pop 1
  x2 <- rbinom(n, 2, prob = x_freq) #genotypes in comparable pop2
  
  b1 <- rnorm(1, sd = rexp(1, 1)) #effect of genotype on gene expression
  
  if(is.na(b2)){b2 <- rnorm(1, sd = rexp(1, 1))} #effect of genotype on phenotype
  if(is.na(b3)){b3 <- rnorm(1, sd = rexp(1, 1))} #effect of expression on phenotype
  
  y1_obs <- x1*b1 + rnorm(n, sd = rexp(1, 1)) #expression in pop 1
  y1_unobs <- x2*b1 + rnorm(n, sd = rexp(1, 1)) #expression in pop 2
  
  y2_obs <- x2*b2 + b3*y1_unobs + rnorm(n, sd = rexp(1, .1))  #phenotype in pop 2
  
  fit1 <- lm(y1_obs ~ x1)
  fit2 <- lm(y2_obs ~ x2)
  
  c(b1 = fit1$coefficients["x1"], b4 = fit2$coefficients["x2"])

}

fits <- t(replicate(1E3, sim(b2, b3)))
plot(fits, pch = 19, col = adjustcolor(1, 0.5))
coef_fit <- lm(fits[,"b4.x2"] ~ fits[,"b1.x1"])
abline(coef_fit, col = 2, lty = 2)
abline(a = b2, b = b3, col = "green", lty = 2)
text(paste0(paste0(rep(" ", 100), collapse = ""), "true value of b3 = ", b3, ", estimated slope = ", round(coef_fit$coefficients["fits[, \"b1.x1\"]"], 3)),
     x = par("usr")[1], y = par("usr")[3], pos = 3)
text(paste0(paste0(rep(" ", 100), collapse = ""), "true value of b2 = ", ifelse(is.na(b2), "Exp. 0", b2), ", estimated intercept = ", round(coef_fit$coefficients["(Intercept)"], 3)),
     x = par("usr")[1], y = par("usr")[4], pos = 1)

#check coverage?
sim2 <- function(n){
  b3 = round(rnorm(1,0,2),2)
  fits <- t(replicate(n, sim(b3 = b3)))
  coef_fit <- rlm(fits[,"b4.x"] ~ 0 + fits[,"b1.x"], maxit = 100)
  in95CI <- abs(summary(coef_fit)$coefficients[1] - b3) < qt(0.975, df = n-1)*summary(coef_fit)$coefficients[2]
  tscore <- (summary(coef_fit)$coefficients[1] - b3) / summary(coef_fit)$coefficients[2]
  return(c(in95CI = in95CI, tscore = tscore))
}

n = 1E2
coverage_sims <- as.data.frame(t(replicate(5000, sim2(n))))
mean(coverage_sims$in95CI)
hist(coverage_sims$tscore, breaks = seq(min(coverage_sims$tscore), max(coverage_sims$tscore), length.out = 50))
hist(pt(coverage_sims$tscore, df = n-1), breaks = seq(0,1,length.out = 20))

#excess association in estimate b4 due to b3?

#explore things in a gaussian dag / covariance framework
plot(dagitty::dagitty("dag{y1 <- x -> y2; y1 -> y2}"), lwd =2 )
n = 1E3
b3 = 1.2
x <- rnorm(n)
b1 <- 0.5
e1 <- rnorm(n)
y1 <- x*b1 + e1
b2 <- 1.25
e2 <- rnorm(n)
y2 <- x*b2 + b3*y1 + e2
cov(x,y2) / b1 / var(x) - b2/b1
lm(y2~x)$coefficients["x"] / b1 - b2/b1

cov(x,y2)

lm(y1~x)
cov(y1,x)/var(x)
lm(y2~x)
cov(y2,x)/var(x)
lm(y2~y1)
cov(y2,y1)/var(y1)
b1*b2*var(x) + b3*var(y1)
cov(y1, y2)

#nicer looking DAG
par(mfrow = c(1,1))
plot(NULL, xlim = c(-2,2), ylim = c(-2,2), frame = F, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
text(x = c(-1, 0, 1), y = c(-1,1,-1), labels = c("Y₁", "X", "Y₂"), family="Helvetica", font = 2, cex = 2)
arrows(x0 = -0.05, y0 = 0.85, x1 = -0.9, y1 = -0.85, lwd = 4)
arrows(x0 = 0.05, y0 = 0.85, x1 = 0.9, y1 = -0.85, lwd = 4)
arrows(x0 = -0.9, y0 = -1, x1 = 0.9, y1 = -1, lwd = 4)
text(x = c(-0.55, 0.55, 0), y = c(0.2,0.2,-1.2), labels = c("β₁", "β₂", "β₃"), family="Helvetica", font = 2, cex = 2, col= "darkgreen" )

#### what if I extended this whole thing into more dimensions? ####
n_y <- 2
bs_wanted = round(rnorm(n = n_y,0,2),2)
corr_mat_ys <- rlkj(n_y, 1)
x_fixed <- cumsum(rnorm(n_obs, sd = sqrt(1 / n_obs))) + rev(cumsum(rnorm(n_obs, sd = sqrt(1 / n_obs))))
sim3 <- function(bs_wanted, corr_mat_ys){
  
  n_obs = 1E4
  n_y <- length(bs_wanted)
  # x <- rnorm(n_obs, sd = 1)
  # x <- cumsum(rnorm(n_obs, sd = sqrt(1 / n_obs))) + rev(cumsum(rnorm(n_obs, sd = sqrt(1 / n_obs)))) #autocorrelate the 'x's
  x <- x_fixed
  bs_fitted <- rnorm(n_y, sd = 1)
  ys_exp <- crossprod(t(x), bs_fitted)
  #sd_ys <- rexp(n_y, 0.1)
  sd_ys <- rep(5, n_y) 
  ys <- t(sapply(1:n_obs, function(obs) 
    rmvnorm(n = 1, mean = ys_exp[obs,], sigma = diag(sd_ys) %*% corr_mat_ys %*% diag(sd_ys))
  ))
  z <- x*b1 + ys %*% bs_wanted + rnorm(n_obs, sd = rexp(1, 0.01))
  
  
  fit1 <- lm(ys ~ x)
  fit2 <- lm(z ~ x)
  
  list(bys = fit1$coefficients["x",], bz = fit2$coefficients["x"])
  
}

fits <- do.call(rbind, lapply(1:100, function(x) sim3(bs_wanted, corr_mat_ys)))
bys <- do.call(rbind, fits[,"bys"])
bzs <- do.call(rbind, fits[,"bz"])

coef_fit <- rlm(bzs ~ bys, maxit = 100)
coef_fit <- lm(bzs ~ bys)
par(mar = c(4,6,4,4))
plot(bs_wanted, y = coef_fit$coefficients[-1], cex.lab = 1.5, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5,
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from z ~ $\\Sigma y_i$"))
for(i in 1:length(coef_fit$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = coef_fit$coefficients[-1][i] + 2*summary(coef_fit)$coefficients[,"Std. Error"][-1][i], 
           x1 = bs_wanted[i], y1 = coef_fit$coefficients[-1][i] - 2*summary(coef_fit)$coefficients[,"Std. Error"][-1][i])}
abline(0,1, col = 2, lty = 2, lwd = 2)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))
legend(lwd = 1, x = "bottomright", title = latex2exp::TeX("correlations in error_{ij} on 'y's"),
       legend = c(latex2exp::TeX("+ r_{ij}"), latex2exp::TeX("- r_{ij}")), lty = c(1,1), col = c("orange","blue"))

combos <- apply(do.call(rbind, strsplit(unique(apply(t(apply(expand.grid(1:n_y, 1:n_y), 1, sort)), 
                                           1, paste0, collapse = "-")), "-")), 2, as.numeric)

for(ri in 1:nrow(combos)){
  segments(x0 = bs_wanted[combos[ri,1]], y0 = coef_fit$coefficients[-1][combos[ri,1]], 
           x1 = bs_wanted[combos[ri,2]], y1 = coef_fit$coefficients[-1][combos[ri,2]],
           lwd = corr_mat_ys[combos[ri,1], combos[ri,2]]^2*5, 
           col = adjustcolor(ifelse(sign(corr_mat_ys[combos[ri,1], combos[ri,2]]) == 1, "orange", "blue"), 0.5))
}

#### ok now let us extend it in a slightly more realistic fashion ####

fastsim_MVN <- function(corr_step_multiplier = 0, n){
  x <- cumsum(rnorm(n, sd = sqrt(1 / n))) + rev(cumsum(rnorm(n, sd = sqrt(1 / n)))) / sqrt(1+1/n)
  x <- (x + rnorm(n, sd = sqrt(corr_step_multiplier / n))) / sqrt(1 + corr_step_multiplier / n)
  x
}

n_y <- 15 #number of genes (or tissues)
bs_wanted = round(rnorm(n = n_y,0,2),2) #effect of each gene expression on z, the phenotype; if high, most of genotype signal should come from here
h2_y <- rbeta(n_y,1,1)
prop_resid_var_z <- rbeta(1,10,10)
corr_mat_ys <- rlkj(n_y, 1)

n_obs = c(1E3, 2E4) #total number of individuals in samples 1 and 2
n_rep <- 10 #total number of loci influencing each gene
x_freq <- replicate(n_y, rbeta(n_rep,1,2)) #allele frequency across loci

#uncorrelated alleles
x1 <- lapply(1:n_y, function(gi) t(sapply(x_freq[,gi], function(xf) rbinom(n_obs[1], 2, prob = xf)))) #genotypes in pop 1
x2 <- lapply(1:n_y, function(gi) t(sapply(x_freq[,gi], function(xf) rbinom(n_obs[2], 2, prob = xf)))) #genotypes in comparable pop2

#autocorrelated alleles
# std_normal_thresholds <- qnorm(x_freq)
# x1 <- lapply(1:n_y, function(gi) replicate(n_obs, c(fastsim_MVN(n = n_rep) > std_normal_thresholds[,gi]) + c(fastsim_MVN(n = n_rep) > std_normal_thresholds[,gi])))
# x2 <- lapply(1:n_y, function(gi) replicate(n_obs, c(fastsim_MVN(n = n_rep) > std_normal_thresholds[,gi]) + c(fastsim_MVN(n = n_rep) > std_normal_thresholds[,gi])))

#effect of each marginal variant at each locus on gene expression
bs_expr <- lapply(1:n_y, function(gi) rnorm(n_rep, sd = 1)) #uncorrelated effects on expression
corr_loci <- 0.1
bs_corrmat <- diag(n_y) + corr_loci - diag(n_y) * corr_loci
bs_expr <- sapply(1:n_rep, function(li) rmvnorm(1, sigma = bs_corrmat)) #correlated effects on expression
bs_expr <- lapply(1:n_y, function(gi) bs_expr[gi,])
y1_exp <- lapply(1:n_y, function(gi) t(bs_expr[[gi]]) %*% x1[[gi]])
y2_exp <- lapply(1:n_y, function(gi) t(bs_expr[[gi]]) %*% x2[[gi]])

sd_ys_obs <- sapply(1:n_y, function(gi) sd(c(y1_exp[[gi]], y2_exp[[gi]]))) #can work out exactly later
sd_ys_remaining <- sqrt(sd_ys_obs^2 * (1-h2_y) / h2_y)

#expression levels in genes for pop1
y1 <- lapply(1:n_y, function(gi) y1_exp[[gi]] + rnorm(n = n_obs[1], 0, sd_ys_remaining[gi]))
resids_ys_1 <- rmvnorm(n = n_obs[1], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y1 <- lapply(1:n_y, function(gi) y1_exp[[gi]] + resids_ys_1[,gi]) #correlated residuals

#expression levels in genes for pop2
y2 <- lapply(1:n_y, function(gi) y2_exp[[gi]] + rnorm(n = n_obs[2], 0, sd_ys_remaining[gi]))
resids_ys_2 <- rmvnorm(n = n_obs[2], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y2 <- lapply(1:n_y, function(gi) y2_exp[[gi]] + resids_ys_2[,gi]) #correlated residuals

#phenotype z for each population
bs_pheno <- lapply(1:n_y, function(gi) rnorm(n_rep, mean = 10, sd = 1))

z1_exp <- sapply(1:n_obs[1], function(indiv)
  sum(sapply(1:n_y, function(gi) sum(x1[[gi]][,indiv] * bs_pheno[[gi]]))) + #direct genetic effect 
  sum(sapply(1:n_y, function(gi) y1[[gi]][,indiv] %*% bs_wanted[[gi]])) #effect of gene expression
)

z2_exp <- sapply(1:n_obs[2], function(indiv)
  sum(sapply(1:n_y, function(gi) sum(x2[[gi]][,indiv] * bs_pheno[[gi]]))) + #direct genetic effect 
    sum(sapply(1:n_y, function(gi) y2[[gi]][,indiv] %*% bs_wanted[[gi]])) #effect of gene expression
)

sd_z <- sqrt(var(c(z1_exp, z2_exp)) * prop_resid_var_z / (1-prop_resid_var_z))
z1 <- z1_exp + rnorm(n_obs[1], 0, sd_z)
z2 <- z2_exp + rnorm(n_obs[2], 0, sd_z)

#fit single regressions
# fit1 <- do.call(rbind, lapply(1:n_y, function(gi) sapply(1:n_rep, function(li) lm(c(y1[[gi]]) ~ x1[[gi]][li,])$coefficients[2])))
# fit2 <- do.call(rbind, lapply(1:n_y, function(gi) sapply(1:n_rep, function(li) lm(z2 ~ x2[[gi]][li,])$coefficients[2])))

#alternatively, fit multiple regressions
fit1 <- do.call(rbind, lapply(1:n_y, function(gi) lm(t(y1[[gi]]) ~ t(x1[[gi]]))$coefficients[-1]))
fit2 <- do.call(rbind, lapply(1:n_y, function(gi) lm(t(t(z2)) ~ t(x2[[gi]]))$coefficients[-1]))

pad_with_0s <- function(fit_mat){
  n_genes <- nrow(fit_mat)
  n_loci <- ncol(fit_mat)
  t(sapply(1:n_genes, function(gi) c(rep(0, (gi-1)*n_loci), fit_mat[gi,], rep(0, (n_genes-gi)*n_loci))))
}

fit1 <- t(pad_with_0s(fit1))
rownames(fit1) <- colnames(fit1) <- NULL
fit1_intercepts <- diag(n_rep * n_y)
# fit1 <- cbind(fit1, fit1_intercepts)
fit2 <- t(t(c(t(fit2))))

#fit coefficients model
standardize <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
par(mfrow = c(1,2), mar = c(5,7,4,4))
fit3 <- rlm(fit2 ~ fit1, na.action = na.exclude)
# plot(standardize(resid(fit3)), standardize(unlist(bs_pheno)))
# fit3 <- lm(fit2 ~ fit1)
fit3_summary <- summary(fit3)
plot(bs_wanted, y = fit3$coefficients[-1], cex.lab = 1.5, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5,
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i}$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = fit3$coefficients[-1][i] + 2*fit3_summary$coefficients[,"Std. Error"][-1][i], 
           x1 = bs_wanted[i], y1 = fit3$coefficients[-1][i] - 2*fit3_summary$coefficients[,"Std. Error"][-1][i])}
abline(0,1, col = 2, lty = 2, lwd = 2)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.02*(par("usr")[4] - par("usr")[3]), labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, fit3$coefficients[-1]), 3))), pos = 2)


combos <- apply(do.call(rbind, strsplit(unique(apply(t(apply(expand.grid(1:n_y, 1:n_y), 1, sort)), 
                                                     1, paste0, collapse = "-")), "-")), 2, as.numeric)

for(ri in 1:nrow(combos)){
  segments(x0 = bs_wanted[combos[ri,1]], y0 = fit3$coefficients[-1][combos[ri,1]], 
           x1 = bs_wanted[combos[ri,2]], y1 = fit3$coefficients[-1][combos[ri,2]],
           lwd = corr_mat_ys[combos[ri,1], combos[ri,2]]^2*5, 
           col = adjustcolor(ifelse(sign(corr_mat_ys[combos[ri,1], combos[ri,2]]) == 1, "orange", "blue"), 0.5))
}


in95CI <- abs(fit3$coefficients[-1] - bs_wanted) < qt(0.975, df = n_y-1)*fit3_summary$coefficients[,"Std. Error"][-1]
mean(in95CI)
tscore <- (fit3$coefficients[-1] - bs_wanted) / fit3_summary$coefficients[,"Std. Error"][-1]
hist(pt(tscore, df = n_y-1), breaks = seq(0,1,length.out = 20))
