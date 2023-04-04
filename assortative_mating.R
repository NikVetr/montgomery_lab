#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x){ out <- exp(x) / (1 + exp(x)); ifelse(is.na(out), 1, out) }
softmax <- function(x) exp(x) / sum(exp(x))

#set simulation parameters
n_indiv <- 200
logit_scale_matches <- 2
softmax_scale_selection <- 0
n_loci <- 2
noise_var_y <- n_loci * 1
noise_var_z <- n_loci * 1
ngen <- 100
p.t0 <- rep(0.5, n_loci)

#sample initial biallelic indiv
x.t0 <- lapply(1:n_loci, function(li) cbind(
  rbinom(n_indiv, 1, prob = p.t0[li]), 
  rbinom(n_indiv, 1, prob = p.t0[li])
))

#find phenotype values
b_y <- rep(1, n_loci)
b_y <- rnorm(n_loci); b_y <- b_y / sqrt(sum(b_y^2)) * n_loci #effects for the mating trait
b_z <- rnorm(n_loci); b_z <- b_z / sqrt(sum(b_z^2)) * n_loci #effects for the mating trait
h2 <- sum(p.t0 * (1-p.t0) * 2 * b_y^2) / (sum(p.t0 * (1-p.t0) * 2 * b_y^2) + noise_var_y)
y <- matrix(NA, nrow = n_indiv, ncol = ngen)
z <- matrix(NA, nrow = n_indiv, ncol = ngen)
y[,1] <- Reduce("+", lapply(seq_along(x.t0), function(xi) apply(x.t0[[xi]], 1, sum) * b_y[xi])) + rnorm(n_indiv, sd = sqrt(noise_var_y))
z[,1] <- Reduce("+", lapply(seq_along(x.t0), function(xi) apply(x.t0[[xi]], 1, sum) * b_z[xi])) + rnorm(n_indiv, sd = sqrt(noise_var_z))
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
    prob_swap <- invlogit((currdist - propdist) * logit_scale_matches)
    swap <- sample(c(T,F), 1, prob = c(prob_swap, 1-prob_swap))
    if(swap){
      foo <- matches[propswap[1],1]
      matches[propswap[1],1] <- matches[propswap[2],1]
      matches[propswap[2],1] <- foo
    }
  }
  
  #check average distance for debugging
  mean(apply(matches, 1, function(xi) dy[xi[1], xi[2]]))
  
  #match future children with current parents to represent selection
  child_pis <- rep(1:nrow(matches), each = 2) #0 reproductive variance
  child_pis <- sample(1:nrow(matches), size = n_indiv, replace = T) #uniform sampling variance
  child_pis <- sample(1:nrow(matches), size = n_indiv, 
                      prob = softmax(apply(matches, 1, function(xi) mean(z[xi,gen_i-1])) * 
                                       softmax_scale_selection),
                      replace = T) #variance propto parental midpoint y
  
  x[[gen_i]] <- lapply(1:n_loci, function(li){
    kids <- do.call(rbind, lapply(child_pis, function(mii){ 
      mi <- matches[mii,]
      cbind(sample(x[[gen_i-1]][[li]][mi[1],], 1, replace = T),
            sample(x[[gen_i-1]][[li]][mi[2],], 1, replace = T))
    }))
  })
  
  #new gen phenotypes
  y[,gen_i] <- Reduce("+", lapply(seq_along(x[[gen_i]]), function(xi) apply(x[[gen_i]][[xi]], 1, sum) * b_y[xi])) + rnorm(n_indiv, sd = sqrt(noise_var_y))
  z[,gen_i] <- Reduce("+", lapply(seq_along(x[[gen_i]]), function(xi) apply(x[[gen_i]][[xi]], 1, sum) * b_z[xi])) + rnorm(n_indiv, sd = sqrt(noise_var_z))
  
  #simulate mutations
  
  
}

corrs_genotypes <- sapply(1:ngen, function(gen_i){
  suppressWarnings(cor(apply(x[[gen_i]][[1]], 1, sum), 
      apply(x[[gen_i]][[2]], 1, sum)))
}) #or do ld score D = p(AB) - P(A)P(B)
corrs_genotypes[is.na(corrs_genotypes)] <- 0

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
D[is.na(D)] <- 0


afreq <- t(sapply(1:ngen, function(gen_i){
  sapply(1:n_loci, function(li) mean(x[[gen_i]][[li]]))
}))


#### plotting ####
plot(1:ngen, corrs_genotypes, type = "l", lwd = 2, ylim = c(-0.1,1),
     xlab = "# of Generation", ylab = "Value of Statistic",
     main = latex2exp::TeX(paste0(ngen, " generations of ", n_indiv, " individuals, ", 
                                  n_loci, " independent loci, ", 
                                  "logit_scale_matches = ", logit_scale_matches, 
                                  ", $h^2_{t_0}$ = ", round(h2, 2))))
lines(1:ngen, D, col = "grey", lwd = 2)
lines(1:ngen, afreq[,1], col = "red")
lines(1:ngen, afreq[,2], col = "blue")
legend("topleft", 
       legend = latex2exp::TeX(c("Genotype Correlation(G$_1$, G$_2$)", "Normalized D\' for Loci 1 & 2", 
                                 "Allele Freq at Locus #1", "Allele Freq at Locus #2")),
       col = c(1, "grey", "red", "blue"), lwd = c(2,2,1,1), lty = 1
       )
       

