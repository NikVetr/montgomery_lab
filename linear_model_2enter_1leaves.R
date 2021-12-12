#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

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
plot(bs_wanted, y = coef_fit$coefficients[-1], cex.lab = 1.25, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5,
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

n_y <- 20 #number of genes (or tissues)
bs_wanted_var <- 4
bs_wanted = round(rnorm(n = n_y,0,sqrt(bs_wanted_var)),2) #effect of each gene expression on z, the phenotype; if high, most of genotype signal should come from here
h2_y <- rbeta(n_y,1,1)
prop_resid_var_z <- rbeta(1,10,10)
corr_mat_ys <- rlkj(n_y, 1)

n_obs = c(2E3, 8E3) #total number of individuals in samples 1 and 2
n_rep <- 200 #total number of loci influencing each gene
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
bs_pheno_mean <- 10
bs_pheno <- lapply(1:n_y, function(gi) rnorm(n_rep, mean = bs_pheno_mean, sd = 1))

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
fit1 <- fit1_orig <- do.call(rbind, lapply(1:n_y, function(gi) lm(t(y1[[gi]]) ~ t(x1[[gi]]))$coefficients[-1]))
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
plot(bs_wanted, y = fit3$coefficients[-1], cex.lab = 1.25, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5,
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


# try out the Bayesian model #
d <- list(fit1 = fit1, 
          fit2 = c(fit2),
          n_y = n_y,
          n_rep = n_rep)

stan_program <- "
functions {
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }
}
data {
  int<lower=0> n_y;
  int<lower=0> n_rep;
  real fit2[n_y * n_rep];
  matrix[n_y * n_rep, n_y] fit1;
}
parameters {
  vector[n_y] beta;
  real<lower=0> beta_var;
  real<lower=0> sigma2;
  vector<lower=0>[n_y] extra_sigma2;
  real<lower=0> extra_sigma2_rate_reg;
  real alpha;
}
model {
  vector[n_y * n_rep] mu;
  vector[n_y * n_rep] sigma;

  alpha ~ normal(0,5);
  sigma2 ~ exponential(0.5);
  beta_var ~ exponential(0.5);
  beta ~ normal(0, sqrt(beta_var));
  extra_sigma2_rate_reg ~ exponential(5);
  extra_sigma2 ~ exponential(extra_sigma2_rate_reg);

  mu = alpha + fit1 * beta;
  sigma = rep_each(sqrt(sigma2 + extra_sigma2 .* square(beta)), n_rep);
  // fit2 ~ student_t(5, mu, sqrt(sigma2));
  fit2 ~ student_t(5, mu, sigma);
  
}
"

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)


out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4)
samps <- data.frame(as_draws_df(out$draws()))
bs_pheno_mean
mean(samps[,grep("alpha", colnames(samps))])
par(mfrow = c(1,1))
plot(bs_wanted, apply(samps[,grep("beta\\.", colnames(samps))],2,mean)); abline(0,1)

bs_wanted_var
mean(samps[,grep("beta_var", colnames(samps))])

#### OK, now let's have it do multi-tissues, one gene, same loci affecting everything ####

## specify simulation parameters
n_y <- 15 #number of genes (or tissues)
bs_wanted_var <- 20 #variance of tissdue-specific effect sizes
bs_wanted = round(rnorm(n = n_y,0,sqrt(bs_wanted_var)),2) #effect of each gene expression on z, the phenotype; if high, most of genotype signal should come from here
h2_y <- rbeta(n_y,15,5) #heritability of gene expression
prop_resid_var_z <- rbeta(1,10,10) #residual variance in trait, conditional on gene expression and SNPs
weight_with_perfect_correlation <- 0.5
corr_mat_ys <- rlkj(n_y, 1) * (1-weight_with_perfect_correlation) + 
  weight_with_perfect_correlation #correlation structure between SNP effects on trait
n_obs = c(800, 1.5E4) #total number of individuals in samples 1 and 2
n_rep <- 100 #total number of loci influencing each gene
x_freq <- rbeta(n_rep,1,2) * 0.99 + 0.01 #allele frequency across loci


# simulate uncorrelated alleles
# x1 <- sapply(x_freq, function(xf) rbinom(n_obs[1], 2, prob = xf)) #genotypes in pop 1
# x2 <- sapply(x_freq, function(xf) rbinom(n_obs[2], 2, prob = xf)) #genotypes in comparable pop2

# simulate autocorrelated alleles by thresholding a MVN random variable 
# rbf_d <- function(dist, params = list(s = 1, l = 1)){params$s^2*exp(-(dist)^2/(2*params$l^2))}
ou_d <- function(dist, params = list(s = 1, l = 1)){params$s^2*exp(-abs(dist)/params$l)}
# brownian <- function(x1, x2, params = list(s=1)){return(params$s^2*min(x1, x2))}
# fastsim_MVN <- function(corr_step_multiplier = 0, n){
#   x <- cumsum(rnorm(n, sd = sqrt(1 / n))) + rev(cumsum(rnorm(n, sd = sqrt(1 / n)))) / sqrt(1+1/n)
#   x <- (x + rnorm(n, sd = sqrt(corr_step_multiplier / n))) / sqrt(1 + corr_step_multiplier / n)
#   x
# }
std_normal_thresholds <- qnorm(x_freq)
x_locations <- sort(runif(n_rep, 0, 2))
x_dists <- as.matrix(dist(x_locations))
x_cov <- ou_d(x_dists, params = list(s = 1, l = 0.15)) #covariance in underlying MVN 

#genotypes in population 1
x1 <- replicate(2, rmvnorm(n = n_obs[1], mean = rep(0, n_rep), sigma = x_cov))
x1 <- sapply(1:n_rep, function(li) (x1[,li,1] > std_normal_thresholds[li]) + (x1[,li,2] > std_normal_thresholds[li]))
#genotypes in population 2
x2 <- replicate(2, rmvnorm(n = n_obs[2], mean = rep(0, n_rep), sigma = x_cov))
x2 <- sapply(1:n_rep, function(li) (x2[,li,1] > std_normal_thresholds[li]) + (x2[,li,2] > std_normal_thresholds[li]))

#effect of each marginal variant at each locus on gene expression
bs_expr <- t(sapply(1:n_rep, function(li) rmvnorm(1, sigma = corr_mat_ys))) #correlated effects on expression

#expected values for gene expression in pops 1 and 2
y1_exp <- sapply(1:n_y, function(gi) x1 %*% t(t(bs_expr[,gi])))
y2_exp <- sapply(1:n_y, function(gi) x2 %*% t(t(bs_expr[,gi])))

#residual variance
sd_ys_obs <- sapply(1:n_y, function(gi) sd(c(y1_exp[,gi], y2_exp[,gi]))) #can work out exactly later
sd_ys_remaining <- sqrt(sd_ys_obs^2 * (1-h2_y) / h2_y)

#expression levels in genes for pop1
# y1 <- sapply(1:n_y, function(gi) y1_exp[,gi] + rnorm(n = n_obs[1], 0, sd_ys_remaining[gi]))
resids_ys_1 <- rmvnorm(n = n_obs[1], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y1 <- sapply(1:n_y, function(gi) y1_exp[,gi] + resids_ys_1[,gi]) #correlated residuals

#expression levels in genes for pop2
# y2 <- sapply(1:n_y, function(gi) y2_exp[,gi] + rnorm(n = n_obs[2], 0, sd_ys_remaining[gi]))
resids_ys_2 <- rmvnorm(n = n_obs[2], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y2 <- sapply(1:n_y, function(gi) y2_exp[,gi] + resids_ys_2[,gi]) #correlated residuals

#direct SNP effects on phenotype z for each population
bs_pheno_mean <- 2
bs_pheno <- rnorm(n_rep, mean = bs_pheno_mean, sd = 1)

z1_exp <- sapply(1:n_obs[1], function(indiv)
  sum(x1[indiv,] * bs_pheno) + #direct genetic effect 
    sum(y1[indiv,] * bs_wanted) #effect of gene expression
)

z2_exp <- sapply(1:n_obs[2], function(indiv)
  sum(x2[indiv,] * bs_pheno) + #direct genetic effect 
    sum(y2[indiv,] * bs_wanted) #effect of gene expression
)

#... can also work out exactly later
sd_z <- sqrt(var(c(z1_exp, z2_exp)) * prop_resid_var_z / (1-prop_resid_var_z))
z1 <- z1_exp + rnorm(n_obs[1], 0, sd_z)
z2 <- z2_exp + rnorm(n_obs[2], 0, sd_z)

#fit single regressions
# fit1 <- do.call(rbind, lapply(1:n_y, function(gi) sapply(1:n_rep, function(li) lm(c(y1[[gi]]) ~ x1[[gi]][li,])$coefficients[2])))
# fit2 <- do.call(rbind, lapply(1:n_y, function(gi) sapply(1:n_rep, function(li) lm(z2 ~ x2[[gi]][li,])$coefficients[2])))

#alternatively, fit multiple regressions, and pretend we got them from the single regressions
fit1 <- do.call(rbind, lapply(1:n_y, function(gi) lm(y1[,gi] ~ x1)$coefficients[-1])) #fitting the OLS eQTL model
fit1_SE <- do.call(rbind, lapply(1:n_y, function(gi) summary(lm(y1[,gi] ~ x1))$coefficients[-1,2])) #not very efficient but w/e
fit2 <-  lm(z2 ~ x2)$coefficients[-1] #fitting the OLS GWAS model
fit2_SE <-  summary(lm(z2 ~ x2))$coefficients[-1,2]

#fit coefficients model
# standardize <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
par(mfrow = c(3,2), mar = c(4,7,3,3), xpd = NA)
fit3 <- rlm(t(t(fit2)) ~ t(fit1), na.action = na.exclude)
# plot(standardize(resid(fit3)), standardize(unlist(bs_pheno)))
# fit3 <- lm(t(t(fit2)) ~ t(fit1), na.action = na.exclude)
fit3_summary <- summary(fit3)
plot(bs_wanted, y = fit3$coefficients[-1], cex.lab = 1.25, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5, main = "Robust IWLS ('M-Estimator')",
     ylim = range(c(fit3$coefficients[-1] + 2*fit3_summary$coefficients[,"Std. Error"][-1], 
                    fit3$coefficients[-1] - 2*fit3_summary$coefficients[,"Std. Error"][-1])),
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i} + \\epsilon$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = fit3$coefficients[-1][i] + 2*fit3_summary$coefficients[,"Std. Error"][-1][i], 
           x1 = bs_wanted[i], y1 = fit3$coefficients[-1][i] - 2*fit3_summary$coefficients[,"Std. Error"][-1][i])}
abline(0,1, col = 2, lty = 2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.03*(par("usr")[4] - par("usr")[3]), labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, fit3$coefficients[-1]), 3))), pos = 2)

#95% CI
in95CI <- abs(fit3$coefficients[-1] - bs_wanted) < qt(0.975, df = n_y-1)*fit3_summary$coefficients[,"Std. Error"][-1]
text(x = par("usr")[2], y = par("usr")[3] + 0.12*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("Prop. in CI_{95} = ", round(mean(in95CI), 3))), pos = 2)

#additional coverage check
tscore <- (fit3$coefficients[-1] - bs_wanted) / fit3_summary$coefficients[,"Std. Error"][-1]
hist(pt(tscore, df = n_y-1), breaks = seq(0,1,length.out = 20), xlab = "Quantile of True Values in T-Distribution")

#plot raw coefficient estimates
plot(bs_expr, t(fit1), pch = 19, col = adjustcolor(1, 0.5), 
     ylim = range(c(t(fit1) + 2*t(fit1_SE), t(fit1) - 2*t(fit1_SE))),
     xlab = "true eQTL effect sizes", ylab = "estimated (OLS) eQTL effect sizes"); 
for(i in 1:length(fit1_SE)){
  segments(x0 = bs_expr[i], y0 = t(fit1)[i] + 2*t(fit1_SE)[i], 
           x1 = bs_expr[i], y1 = t(fit1)[i] - 2*t(fit1_SE)[i],
           col = adjustcolor(1, 0.5))}
abline(0,1,col=2,lty=2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))

plot(bs_pheno, fit2, pch = 19, col = adjustcolor(1, 0.5), 
     ylim = range(fit2 + 2*fit2_SE, fit2 - 2*fit2_SE),
     xlab = "true horizontal pleiotropic (i.e. direct) effect sizes", ylab = "estimated (OLS) total SNP effect sizes"); 
for(i in 1:length(bs_pheno)){
  segments(x0 = bs_pheno[i], y0 = (fit2)[i] + 2*(fit2_SE)[i], 
           x1 = bs_pheno[i], y1 = (fit2)[i] - 2*(fit2_SE)[i],
           col = adjustcolor(1, 0.5))}
abline(0,1,col=2,lty=2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))

bs_pheno_mean
fit3$coefficients[1]

# try out the Bayesian model #
d <- list(fit1 = t(fit1),
          fit1_SE = t(fit1_SE),
          fit2 = c(fit2),
          fit2_SE = c(fit2_SE),
          n_y = n_y,
          n_rep = n_rep)

stan_program <- "
data {
  int<lower=0> n_y;
  int<lower=0> n_rep;
  real fit2[n_rep];
  real fit2_SE[n_rep];
  matrix[n_rep, n_y] fit1;
  matrix[n_rep, n_y] fit1_SE;
}
parameters {
  vector[n_y] beta;
  real<lower=0> beta_var;
  real<lower=0> sigma2;
  real alpha;
  matrix[n_rep, n_y] fit1_EST;
  real fit2_EST[n_rep];
}
model {
  // intermediate variables
  vector[n_y * n_rep] mu;
  real sigma;

  // priors
  alpha ~ normal(0,5);
  sigma2 ~ exponential(0.1);
  beta_var ~ exponential(0.1);
  beta ~ normal(0, sqrt(beta_var));
  
  // upstream inferential uncertainty, since normals are symmetric
  // normal() doesn't seem to accept matrix arguments but does accept vector args
  for(ri in 1:n_rep){
    fit1[ri,] ~ normal(fit1_EST[ri,], fit1_SE[ri,]);  
  }
  fit2 ~ normal(fit2_EST, fit2_SE);

  // actual model and likelihood
  mu = alpha + fit1_EST * beta;
  sigma = sqrt(sigma2);
  // fit2_EST ~ student_t(5, mu, sigma);
  fit2_EST ~ normal(mu, sigma);
  
}
"

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 4E3, iter_warmup = 4E3, data = d, parallel_chains = 4)
samps <- data.frame(as_draws_df(out$draws()))

#print out some hyperparameter estimates
bs_pheno_mean
mean(samps[,grep("alpha", colnames(samps))])

bs_wanted_var
mean(samps[,grep("beta_var", colnames(samps))])

#visualize posterior means and sds
plot(bs_wanted, y = apply(samps[,grep("beta\\.", colnames(samps))],2,mean), cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Hierarchical Bayesian Analysis",
     ylim = range(c(apply(samps[,grep("beta\\.", colnames(samps))],2,mean) + 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd),
                    apply(samps[,grep("beta\\.", colnames(samps))],2,mean) - 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))),
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i} + \\epsilon$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = (apply(samps[,grep("beta\\.", colnames(samps))],2,mean) + 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))[i],
           x1 = bs_wanted[i], y1 = (apply(samps[,grep("beta\\.", colnames(samps))],2,mean) - 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))[i])}
abline(0,1, col = 2, lty = 2, lwd = 2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SD"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.02*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, apply(samps[,grep("beta\\.", colnames(samps))],2,mean)), 3))), pos = 2)
#check 95% coverage
in95CI_Bayes <- apply(samps[,grep("beta\\.", colnames(samps))],2,quantile, probs = c(0.025,0.975))
in95CI_Bayes <- bs_wanted > in95CI_Bayes[1,] & bs_wanted < in95CI_Bayes[2,]  
text(x = par("usr")[2], y = par("usr")[3] + 0.12*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("Prop. in CI_{95} = ", round(mean(in95CI_Bayes), 3))), pos = 2)


#compare bayesian vs OLS estimates
plot(fit3$coefficients[-1], apply(samps[,grep("beta\\.", colnames(samps))],2,mean),
     ylab = "Posterior Means", xlab = "Robust IWLS Estimates", cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Shrinkage Comparison") 
abline(0,1, lty = 2, col = 2, xpd = F)
abline(lm(apply(samps[,grep("beta\\.", colnames(samps))],2,mean) ~ fit3$coefficients[-1]), col = 4, lty = 3, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "OLS fit to Points"), lty = c(2,1), col = c(2,4))


