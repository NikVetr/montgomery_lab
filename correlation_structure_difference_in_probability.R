#libraries
library(parallel)
library(mvtnorm)

#### load functions ####
source(file = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts/deg-trait_functions.R")

#custom functions for this script
mcprint <- function(...){
  system(sprintf('echo "%s"', paste0(..., collapse="")))
}

pbivnorm_infs <- function(upper1, upper2, rhos){
  nInf <- (upper1 == -Inf | upper2 == -Inf)
  Inf1 <- upper1 == Inf
  Inf2 <- upper2 == Inf
  n2_1 <- Inf1 & !Inf2
  n1_2 <- Inf2 & !Inf1
  include <- !(nInf | Inf1 | Inf2)
  
  results <- rep(0, length(include))
  results[Inf1 & Inf2] <- 1
  results[n2_1] <- pnorm(q = upper2[n2_1])
  results[n1_2] <- pnorm(q = upper1[n1_2])
  results[include] <- pbivnorm(upper1[include], upper2[include], rhos[include])
  return(results)
}

#set high-level simulation parameters
n <- 1E4 #number of individuals
ncol <- 100 #number of columns, e.g. hobbies
nrow <- 2 #number of rows, e.g. subjects

#sample a correlation matrix for columns / hobbies 'rcol' with some block structure
nblocks <- 20 
baseline_correlation <- 0 #correlation among hobbies outside of blocks
block_correlation <- runif(nblocks, 0, 1) #correlation among hobbies inside of each of the nblocks
rcol <- diag(ncol) + baseline_correlation - diag(ncol) * baseline_correlation
block_cuts <- round(cumsum(c(0,gtools::rdirichlet(1, rep(1, nblocks-1)))) * (ncol-1))
for(i in 1:(nblocks-1)){
  dims <- (block_cuts[i+1]) - (block_cuts[i]+1) + 1
  if(dims == 0){next()}
  rcol[(block_cuts[i]+1):(block_cuts[i+1]),(block_cuts[i]+1):(block_cuts[i+1])] <- diag(dims) + block_correlation[i] - diag(dims) * block_correlation[i]
}

#or just use one that has all the correlations along a linear gradient
rcol <- 1 - as.matrix(dist(cbind(1:ncol, 1:ncol), method = "manhattan")) / (ncol * 2)
for(i in seq(1, ncol, by = 2)){ #prob some linear algebraic way to do this
  rcol[i,] <- -rcol[i,]
  rcol[,i] <- -rcol[,i]
}


simulate_hits <- function(nrowhits, rcol){
  #retrieve high-level parameters
  nrow <- length(nrowhits)
  ncol <- dim(rcol)[1]
  
  #simulate which individuals are hits in the row
  rowhits <- lapply(nrowhits, function(nh) sample(1:n, size = nh, replace = F))
  
  #simulate which individuals are hits int he columns
  liabilities <- rmvnorm(n, sigma = rcol)
  ncolhits <- ceiling(invlogit(rmvnorm(1, sigma = rcol * 2)) * n)
  thresholds <- sapply(1:ncol, function(ci) sort(liabilities[,ci])[n-ncolhits[ci]])
  colhitsmat <- sapply(1:ncol, function(ci) as.numeric(liabilities[,ci] > thresholds[ci]))
  colhits <- lapply(1:ncol, function(ci) which(colhitsmat[,ci] == 1))
  
  #aggregate to table
  hits_in_row <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], rowhits[[ri]]))))
  hits_not_in_row <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], setdiff(1:n, rowhits[[ri]])))))
  y <- rbind(hits_in_row, hits_not_in_row)
  
  return(list(y = y, colhits = colhits, rowhits = rowhits))
}

prop_teachers <- "1 / 500"
nbeta <- 1E3
bp <- c(1,1)
nrowhits <- round(n * eval(parse(text=prop_teachers))) #number of individuals in row, e.g. teachers in subject
MLE_diffs_in_prop <- do.call(rbind, mclapply(1:5E3, function(i) {
  y <- simulate_hits(nrowhits, rcol)$y
  prop_hits <- (y + 1) / (rbind(rep(nrowhits, ncol), rep(n-nrowhits, ncol)) + 2) #posterior means of each value
  diff_posterior_means <- apply(logit(prop_hits), 2, diff)
  # posterior_mean_diff <- sapply(1:ncol, function(ci) mean(logit(rbeta(nbeta, 1 + y[1,ci], 1 + nrowhits - y[1,ci])) - logit(rbeta(nbeta, 1 + y[2,ci], 1 + (n - nrowhits) - y[2,ci]))))
  posterior_mean_diff <- sapply(1:ncol, function(ci) 
    mean(qnorm(rbeta(nbeta, bp[1] + y[1,ci], bp[2] + nrowhits - y[1,ci])) - qnorm(rbeta(nbeta, bp[1] + y[2,ci], bp[2] + (n - nrowhits) - y[2,ci])))
  )
  posterior_mean_diff
}, mc.cores = 12))
cor_among_diffs <- cor(MLE_diffs_in_prop)

#### plot output ####

#hobby version
plot(cor_among_diffs[upper.tri(cor_among_diffs)], rcol[upper.tri(rcol)], xlim = c(-1,1), ylim = c(-1,1),
     xlab = paste0("pairwise hobby_i x hobby_j correlations among posterior means for differences in probability\n",
                   "of hits (probit-scale) between two schools computed across independent simulations"),
     ylab = "true correlation between hobby liabilities", 
     main = paste0("relationship when ", prop_teachers, " of ", n," individuals are teachers,\n",
                   "priors on Pr(hit) in each group were Beta(", bp[1], ", ", bp[2], ")s"), 
     pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, lty = 2, lwd = 2, col = 2)
legend(x = "topleft", legend = c("1-to-1 line", "pairwise correlation between hobbies"), lty = c(2, NA), lwd = c(2, NA), col = c(2, adjustcolor(1, 0.5)), pch = c(NA, 19))

#bio version
plot(cor_among_diffs[upper.tri(cor_among_diffs)], rcol[upper.tri(rcol)], xlim = c(-1,1), ylim = c(-1,1),
     xlab = paste0("pairwise trait_i x trait_j correlations among posterior means for differences in probability\n",
                   "of hits (probit-scale) between two tissues computed across independent simulations"),
     ylab = "correlation between trait liabilities", 
     main = paste0("relationship when ", prop_teachers, " of ", n," genes are DEGs,\n",
                   "priors on Pr(hit) in each group were Beta(", bp[1], ", ", bp[2], ")s"), 
     pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, lty = 2, lwd = 2, col = 2)
legend(x = "topleft", legend = c("1-to-1 line", "pairwise correlation between traits"), lty = c(2, NA), lwd = c(2, NA), col = c(2, adjustcolor(1, 0.5)), pch = c(NA, 19))


#estimate bivariate correlations in a single set of outcomes
y <- simulate_hits(nrowhits, rcol)$colhits

probs_quadrants <- function(mus, rij){
  probs <- matrix(0, 2, 2)
  probs[2,1] <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[1,1] <- mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[1,2] <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[2,2] <- 1 - sum(probs)
  return(probs)
}

ll_biv_probit <- function(data, par) {
  
  #undo params
  rij <- par[1]
  
  #undo data
  ct2x2 <- data$ct2x2
  mus <- data$mus
  
  #compute log likelihood
  probs <- probs_quadrants(mus, rij)
  log_probs <- log(probs)
  ll <- sum(log_probs * ct2x2)
  
  #return log likelihood
  return(-ll)
  
}

ll_biv_probit_mus <- function(data, par) {
  
  #undo params
  rij <- par[1]
  mus <- par[2:3]
  
  #undo data
  ct2x2 <- data$ct2x2
  
  #compute log likelihood
  probs <- probs_quadrants(mus, rij)
  log_probs <- log(probs)
  ll <- sum(log_probs * ct2x2)
  
  #return log likelihood
  return(-ll)
  
}

estimate_correlations <- function(y, n, optimize_jointly = F, ncores = 12, nearPDadj = T, print_progress = F) {
  
  #retrieve metadata
  ncol <- length(y)
  
  #iterate through all pairs of dimensions, optimizing bivariate probit probability
  estimated_corrs <- mclapply(1:(ncol-1), function(ci1) sapply((ci1+1):ncol, function(ci2){
  
    #print current location
    if(print_progress){
      mcprint(paste0(ci1, ", ", ci2, ", ", optim_out$p1))  
    }
    
    #snag 2x2 contingency table
    ct2x2 <- matrix(c(length(setdiff(y[[ci1]],y[[ci2]])), 
                      n - length(union(y[[ci1]],y[[ci2]])), 
                      length(intersect(y[[ci1]],y[[ci2]])), 
                      length(setdiff(y[[ci2]],y[[ci1]]))), 
                    2, 2)
    
    #can optimize jointly or marginally, knowing the location of the mus
    if(optimize_jointly){
      optim_out <- optimx::optimx(par = rep(0,3), 
                     fn = ll_biv_probit_mus, 
                     data = list(ct2x2 = ct2x2), method = "nlm", lower = c(-1,-Inf,-Inf), upper = c(1,Inf,Inf),
                     control = list(maxit = 1E3, trace = 0, dowarn = F))
    } else {
      mus <- qnorm(c(sum(ct2x2[1,]), sum(ct2x2[,2])) / sum(ct2x2))
      optim_out <- optimx::optimx(par = 0.5, 
                                  fn = ll_biv_probit, 
                                  data = list(ct2x2 = ct2x2, mus = mus), method = "nlm", lower = -1, upper = 1,
                                  control = list(maxit = 1E3, trace = 0, dowarn = F))
    }
    
    #return estimate
    optim_out$p1
    
  }), mc.cores = ncores)
  
  #reconstruct and (optionally) adjust estimated matrix
  estimated_corr_mat <- matrix(0, nrow = ncol, ncol = ncol)
  for(ri in 1:(ncol-1)){
    estimated_corr_mat[ri, (ri+1):ncol] <- estimated_corrs[[ri]]
  }
  estimated_corr_mat <- estimated_corr_mat + t(estimated_corr_mat) + diag(ncol)
  if(nearPDadj){
    estimated_corr_mat <- Matrix::nearPD(estimated_corr_mat, corr = T)$mat  
  }
  
  return(estimated_corr_mat)

}

estimated_corr_mat <- estimate_correlations(y, n)

plot(estimated_corr_mat, rcol); abline(0,1)


plot(cor(colhitsmat)[upper.tri(cor(colhitsmat))], rcol[upper.tri(rcol)], xlim = c(-1,1)); abline(0,1)
plot(cor(colhitsmat)[upper.tri(cor(colhitsmat))], cor_among_diffs[upper.tri(cor_among_diffs)], xlim = c(-1,1)); abline(0,1)


#see if the correlation among individuals for a given trait corresponds to the correlation matrix of interest
