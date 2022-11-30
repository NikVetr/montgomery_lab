#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x){ out <- exp(x) / (1 + exp(x)); ifelse(is.na(out), 1, out) }

#set simulation parameters
n_indiv <- 1000
logit_scale <- 1
n_loci <- 10
noise_var <- n_loci * 1
ngen <- 100
p.t0 <- rep(0.5, n_loci)

#sample initial biallelic indiv
x.t0 <- lapply(1:n_loci, function(li) cbind(
  rbinom(n_indiv, 1, prob = p.t0[li]), 
  rbinom(n_indiv, 1, prob = p.t0[li])
))

#find phenotype values
b <- rep(1, n_loci)
y <- matrix(NA, nrow = n_indiv, ncol = ngen)
y[,1] <- Reduce("+", lapply(x.t0, function(xi) apply(xi, 1, sum) * b)) + rnorm(n_indiv, sd = sqrt(noise_var))
x <- lapply(1:ngen, function(i) x.t0)

gen_i <- 2

for(gen_i in 2:ngen){
  if(gen_i %% 10 == 0){print(gen_i)}
  
  # sample mating pairs -- initialize to optimum, then specify a markov chain for swaps
  matches <- matrix(order(y[,gen_i-1]), ncol = 2, byrow = T)
  dy <- as.matrix(dist(y[,gen_i-1]))
  diag(dy) <- Inf
  niter <- 2E1 * n_indiv / 2
  for(i in 1:niter){
    propswap <- sample(1:(n_indiv/2), 2, replace = F)
    currdist <- dy[matches[propswap[1],1], matches[propswap[1],2]] + dy[matches[propswap[2],1], matches[propswap[2],2]]
    propdist <- dy[matches[propswap[1],1], matches[propswap[2],2]] + dy[matches[propswap[1],2], matches[propswap[2],1]]
    prob_swap <- invlogit((currdist - propdist) * logit_scale)
    swap <- sample(c(T,F), 1, prob = c(prob_swap, 1-prob_swap))
    if(swap){
      foo <- matches[propswap[1],1]
      matches[propswap[1],1] <- matches[propswap[2],1]
      matches[propswap[2],1] <- foo
    }
  }
  
  #check average distance for debugging
  mean(apply(matches, 1, function(xi) dy[xi[1], xi[2]]))
  
  #new generation
  x[[gen_i]] <- lapply(1:n_loci, function(li){
    gametes <- t(apply(x[[gen_i-1]][[li]], 1, sample, replace = T))
    kids <- do.call(rbind, apply(matches, 1, function(mi) 
                      rbind(c(gametes[mi[1], 1], gametes[mi[2], 1]), 
                            c(gametes[mi[1], 2], gametes[mi[2], 2])), 
                      simplify = F))
  })
  
  #new gen phenotypes
  y[,gen_i] <- Reduce("+", lapply(x[[gen_i]], function(xi) apply(xi, 1, sum) * b)) + rnorm(n_indiv, sd = sqrt(noise_var))
  
}

corrs_genotypes <- sapply(1:ngen, function(gen_i){
  cor(apply(x[[gen_i]][[1]], 1, sum), 
      apply(x[[gen_i]][[2]], 1, sum))
}) #or do ld score D = p(AB) - P(A)P(B)

D <- sapply(1:ngen, function(gen_i){
  D_unnorm <- (mean(x[[gen_i]][[1]][,1] + x[[gen_i]][[2]][,1] == 2) - mean(x[[gen_i]][[1]][,1]) * mean(x[[gen_i]][[2]][,1]) +
     mean(x[[gen_i]][[1]][,2] + x[[gen_i]][[2]][,2] == 2) - mean(x[[gen_i]][[1]][,2]) * mean(x[[gen_i]][[2]][,2]) +
     mean(x[[gen_i]][[1]][,1] + x[[gen_i]][[2]][,2] == 2) - mean(x[[gen_i]][[1]][,1]) * mean(x[[gen_i]][[2]][,2]) +
     mean(x[[gen_i]][[1]][,2] + x[[gen_i]][[2]][,1] == 2) - mean(x[[gen_i]][[1]][,2]) * mean(x[[gen_i]][[2]][,1])) / 4
  D_norm <- if(D_unnorm > 0){
    D_unnorm / min(c(mean(x[[gen_i]][[1]]) * (1 - mean(x[[gen_i]][[2]]))), 
        c(mean(x[[gen_i]][[2]]) * (1 - mean(x[[gen_i]][[1]]))))
  } else {
    D_unnorm / max(c(-mean(x[[gen_i]][[1]])*mean(x[[gen_i]][[2]]),
          -(1-mean(x[[gen_i]][[1]]))*(1-mean(x[[gen_i]][[2]]))))
  }
  D_norm
})

afreq <- t(sapply(1:ngen, function(gen_i){
  sapply(1:n_loci, function(li) mean(x[[gen_i]][[li]]))
}))


#### plotting ####
plot(1:ngen, corrs_genotypes, type = "l", lwd = 2, ylim = c(-0.1,1),
     xlab = "# of Generation", ylab = "Value of Statistic",
     main = latex2exp::TeX(paste0(ngen, " generations of ", n_indiv, " individuals, ", 
                                  n_loci, " independent loci, ", 
                                  "logit_scale = ", logit_scale, 
                                  ", $h^2_{t_0}$ = ", round(sum(p.t0 * (1-p.t0) * 2 * b^2) / (sum(p.t0 * (1-p.t0) * 2 * b^2) + noise_var), 2))))
lines(1:ngen, D, col = "grey", lwd = 2)
lines(1:ngen, afreq[,1], col = "red")
lines(1:ngen, afreq[,2], col = "blue")
legend("topleft", 
       legend = latex2exp::TeX(c("Genotype Correlation(G$_1$, G$_2$)", "Normalized D\' for Loci 1 & 2", 
                                 "Allele Freq at Locus #1", "Allele Freq at Locus #2")),
       col = c(1, "grey", "red", "blue"), lwd = c(2,2,1,1), lty = 1
       )
       

