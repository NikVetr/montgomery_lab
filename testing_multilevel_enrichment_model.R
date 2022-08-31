#libraries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(mvtnorm)

#functions
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")

#high level spec
n <- 1E4
ncol <- 100
nrow <- 20
ncolhits <- sample(1:ceiling(n/2), size = ncol, replace = T)
nrowhits <- sample(ceiling(n/100):ceiling(n/20), size = nrow, replace = T)
correlated_cols <- T
col_nblocks <- 1
col_baseline_correlation <- 0.1
col_block_correlation <- 0.95

correlated_rows <- F
row_nblocks <- 1
row_baseline_correlation <- 0.1
row_block_correlation <- 0.95

#simulate data
colhits <- lapply(ncolhits, function(nh) sample(1:n, size = nh, replace = F))
rowhits <- lapply(nrowhits, function(nh) sample(1:n, size = nh, replace = F))

#simulate correlated colhits?
if(correlated_cols){
  jacmat_col_1 <- sapply(colhits, function(x1) sapply(colhits, function(x2) jaccard(x1, x2)))
  if(col_nblocks > 1){
    rcol <- diag(ncol) + col_baseline_correlation - diag(ncol) * col_baseline_correlation
    block_cuts <- round(cumsum(c(0,gtools::rdirichlet(1, rep(1, col_nblocks-1)))) * (ncol-1))
    for(i in 1:(col_nblocks-1)){
      dims <- (block_cuts[i+1]) - (block_cuts[i]+1) + 1
      if(dims == 0){next()}
      rcol[(block_cuts[i]+1):(block_cuts[i+1]),(block_cuts[i]+1):(block_cuts[i+1])] <- diag(dims) + col_block_correlation - diag(dims) * col_block_correlation
    }
  } else {
    rcol <- diag(ncol) + col_block_correlation - diag(ncol) * col_block_correlation
  }
  liabilities_col <- rmvnorm(n, sigma = rcol)
  ncolhits <- ceiling(invlogit(rmvnorm(1, sigma = rcol * 2)) * n / 3)
  thresholds <- sapply(1:ncol, function(ci) sort(liabilities_col[,ci])[n-ncolhits[ci]])
  colhits <- lapply(1:ncol, function(ci) which(liabilities_col[,ci] > thresholds[ci]))
  jacmat_col_2 <- sapply(colhits, function(x1) sapply(colhits, function(x2) jaccard(x1, x2)))
  par(mfrow = c(2,1))
  hist(jacmat_col_1[upper.tri(jacmat_col_1)], breaks = c(0:100/100), xlim = c(0,1))
  hist(jacmat_col_2[upper.tri(jacmat_col_2)], breaks = c(0:100/100), xlim = c(0,1))
}

#simulate correlated rowhits?
if(correlated_rows){
  jacmat_row_1 <- sapply(rowhits, function(x1) sapply(rowhits, function(x2) jaccard(x1, x2)))
  if(row_nblocks > 1){
    rrow <- diag(nrow) + row_baseline_correlation - diag(nrow) * row_baseline_correlation
    block_cuts <- round(cumsum(c(0,gtools::rdirichlet(1, rep(1, row_nblocks-1)))) * (nrow-1))
    for(i in 1:(row_nblocks-1)){
      dims <- (block_cuts[i+1]) - (block_cuts[i]+1) + 1
      if(dims == 0){next()}
      rrow[(block_cuts[i]+1):(block_cuts[i+1]),(block_cuts[i]+1):(block_cuts[i+1])] <- diag(dims) + row_block_correlation - diag(dims) * row_block_correlation
    }
  } else {
    rrow <- diag(nrow) + row_block_correlation - diag(nrow) * row_block_correlation
  }
  liabilities_row <- rmvnorm(n, sigma = rrow)
  nrowhits <- ceiling(invlogit(rmvnorm(1, sigma = rrow * 2)) * n / 3)
  thresholds <- sapply(1:nrow, function(ci) sort(liabilities_row[,ci])[n-nrowhits[ci]])
  rowhits <- lapply(1:nrow, function(ci) which(liabilities_row[,ci] > thresholds[ci]))
  jacmat_row_2 <- sapply(rowhits, function(x1) sapply(rowhits, function(x2) jaccard(x1, x2)))
  par(mfrow = c(2,1))
  hist(jacmat_row_1[upper.tri(jacmat_row_1)], breaks = c(0:100/100), xlim = c(0,1))
  hist(jacmat_row_2[upper.tri(jacmat_row_2)], breaks = c(0:100/100), xlim = c(0,1))
}

bothhits <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], rowhits[[ri]]))))

data <- data.frame(do.call(rbind, lapply(1:ncol, function(ci) do.call(rbind, lapply(1:nrow, function(ri) 
  c(count = length(intersect(colhits[[ci]], rowhits[[ri]])),
    total = n,
    row_count = nrowhits[ri],
    col_count = ncolhits[ci],
    row_index = ri,
    col_index = ci,
    row_n = nrow,
    col_n = ncol)))
)))

#sample hits for each population
d <- list(cell_count = data$count,
          total = data$total,
          row_count = data$row_count,
          col_count = data$col_count,
          row_index = data$row_index,
          col_index = data$col_index,
          row_n = nrow,
          col_n = ncol)

base = paste0("simple_test")
stan_program <- '
data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=0> total[row_n * col_n];
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=0> row_count[row_n * col_n];
    int<lower=0> col_count[row_n * col_n];
    int<lower=0> cell_count[row_n * col_n];
}
transformed data {
    int<lower=1> n = row_n * col_n;
    int<lower=0> cell_count_compl[n];
    int<lower=0> row_count_compl[n];
    for(i in 1:n){
      cell_count_compl[i] = col_count[i] - cell_count[i];
      row_count_compl[i] = total[i] - row_count[i];
    }
}
parameters {
    //col params
    real col_mean;
    real<lower=0> col_sd;
    vector[col_n] raw_col_logodds;
    real<lower=0> cell_sd;
    vector[n] raw_cell_logodds;

    //differences in deviations terms
    real overall_difference;
    vector[row_n] raw_row_difference;
    vector[col_n] raw_col_difference;
    vector[n] raw_cell_difference;
    
    real<lower=0> row_difference_sd;
    real<lower=0> col_difference_sd;
    real<lower=0> cell_difference_sd;
}
transformed parameters {
    //recenter params
    vector[col_n] col_logodds = raw_col_logodds * col_sd;
    vector[n] cell_logodds = raw_cell_logodds * cell_sd + col_logodds[col_index];

    //incorporate difference
    vector[col_n] col_difference = raw_col_difference * col_difference_sd;
    vector[row_n] row_difference = raw_row_difference * row_difference_sd;
    vector[n] cell_difference = raw_cell_difference * cell_difference_sd;
    vector[n] cell_logodds_focal = cell_logodds +
              (overall_difference + row_difference[row_index] + col_difference[col_index] + cell_difference) / 2;
    vector[n] cell_logodds_compl = cell_logodds -
              (overall_difference + row_difference[row_index] + col_difference[col_index] + cell_difference) / 2;
}
model {
    //priors and hyperpriors

    //marginal params
    col_mean ~ normal(0,2);
    col_sd ~ std_normal();
    raw_col_logodds ~ std_normal();
    raw_cell_logodds ~ std_normal();
    cell_sd ~ std_normal();

    //difference params
    overall_difference ~ std_normal();

    raw_col_difference ~ std_normal();
    col_difference_sd ~ std_normal();

    raw_row_difference ~ std_normal();
    row_difference_sd ~ std_normal();

    raw_cell_difference ~ std_normal();
    cell_difference_sd ~ std_normal();

    //likelihood
    cell_count ~ binomial_logit(row_count, cell_logodds_focal);
    cell_count_compl ~ binomial_logit(row_count_compl, cell_logodds_compl);

}
generated quantities {
    vector[n] cell_total_prob_difference = inv_logit(cell_logodds_focal) - inv_logit(cell_logodds_compl);
}
'

stan_program <- '
data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=0> total[row_n * col_n];
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=0> row_count[row_n * col_n];
    int<lower=0> col_count[row_n * col_n];
    int<lower=0> cell_count[row_n * col_n];
}
transformed data {
    int<lower=1> n = row_n * col_n;
    int<lower=0> cell_count_compl[n];
    int<lower=0> row_count_compl[n];
    for(i in 1:n){
      cell_count_compl[i] = col_count[i] - cell_count[i];
      row_count_compl[i] = total[i] - row_count[i];
    }
}
parameters {
    //col params
    real col_mean;
    real<lower=0> col_sd;
    vector[col_n] raw_col_logodds;
    real<lower=0> cell_sd;
    vector[n] raw_cell_logodds;

    //differences in deviations terms
    vector[col_n] raw_col_difference;
    vector[n] raw_cell_difference;
    
    real<lower=0> col_difference_sd;
    real<lower=0> cell_difference_sd;
}
transformed parameters {
    //recenter params
    vector[col_n] col_logodds = raw_col_logodds * col_sd;
    vector[n] cell_logodds = raw_cell_logodds * cell_sd + col_logodds[col_index];

    //incorporate difference
    vector[col_n] col_difference = raw_col_difference * col_difference_sd;
    vector[n] cell_difference = raw_cell_difference * cell_difference_sd;
    vector[n] cell_logodds_focal = cell_logodds +
              (col_difference[col_index] + cell_difference) / 2;
    vector[n] cell_logodds_compl = cell_logodds -
              (col_difference[col_index] + cell_difference) / 2;
}
model {
    //priors and hyperpriors

    //marginal params
    col_mean ~ normal(0,2);
    col_sd ~ std_normal();
    raw_col_logodds ~ std_normal();
    raw_cell_logodds ~ std_normal();
    cell_sd ~ std_normal();

    //difference params
    raw_col_difference ~ std_normal();
    col_difference_sd ~ std_normal();

    raw_cell_difference ~ std_normal();
    cell_difference_sd ~ std_normal();

    //likelihood
    cell_count ~ binomial_logit(row_count, cell_logodds_focal);
    cell_count_compl ~ binomial_logit(row_count_compl, cell_logodds_compl);

}
generated quantities {
    vector[n] cell_total_prob_difference = inv_logit(cell_logodds_focal) - inv_logit(cell_logodds_compl);
}
'

#andrew's model
# colhits
# rowhits
# n
# ncol
# nrow

data2 <- data.frame(i = rep(1:n, times = ncol), h = rep(1:ncol, each = n))
data2$y <- 0
data2$y[unlist(lapply(0:(ncol-1), function(x) x * n + colhits[[x+1]]))] <- 1
smat <- do.call(cbind, lapply(1:nrow, function(x){
  foo <- rep(0, n)
  foo[rowhits[[x]]] <- 1
  foo
}))
colnames(smat) <- paste0("s", 1:nrow)
data2 <- cbind(data2, smat[rep(1:n, times = ncol),])

d <- list(
  I = n,
  H = ncol,
  S = nrow,
  y = data2$y,
  i = data2$i,
  h = data2$h,
  s = smat,
  t = as.numeric(apply(smat, 1, sum) > 0.1)
)

cor(do.call(cbind, split(data2$y, data2$h)))

base = paste0("gelman_model")
stan_program <- '
data {
    int<lower=1> I; //# of individuals
    int<lower=1> H; //# of hobbies
    int<lower=1> S; //# of subjects
    int<lower=0,upper=1> y[I * H]; //whether individual i partakes in hobby h
    int<lower=1,upper=I> i[I * H]; //individual index for y
    int<lower=1,upper=H> h[I * H]; //hobby index for y
    matrix<lower=0.0,upper=1.0>[I, S] s; //subject indicators for 1:I
    vector<lower=0.0,upper=1.0>[I] t; //teacher indicator for 1:I
}
transformed data {
    
}
parameters {
    //main parameters
    real b_0; //avg hobby frequency
    real b_T; //average teacher effect
    vector[H] raw_b_H; //non-centered hobby prob
    vector[H] raw_b_TxH; //non-centered teacher x hobby effect
    vector[S] raw_b_S; //non-centered average subject effect
    matrix[S,H] raw_b_SxH; //non-centered subject x hobby effect
    
    //hyperparameters
    real<lower=0> sd_H;
    real<lower=0> sd_TxH;
    real<lower=0> sd_S;
    real<lower=0> sd_SxH;
}
transformed parameters {
    vector[H] b_H = raw_b_H * sd_H;
    vector[H] b_TxH = raw_b_TxH * sd_TxH;
    vector[S] b_S = raw_b_S * sd_S;
    matrix[S,H] b_SxH = raw_b_SxH * sd_SxH;
}
model {
    //priors
    b_0 ~ normal(0,2);
    b_T ~ std_normal();
    raw_b_H ~ std_normal();
    raw_b_TxH ~ std_normal();
    raw_b_S ~ std_normal();
    to_vector(raw_b_SxH) ~ std_normal();
     
    //hyperpriors
    sd_H ~ std_normal();
    sd_TxH ~ std_normal();
    sd_S ~ std_normal();
    sd_SxH ~ std_normal();
    
    //likelihood
    y ~ bernoulli_logit(b_0 + 
                  b_T * t[i] + 
                  b_H[h] + 
                  b_TxH[h] .* t[i] + 
                  to_vector(s * b_S)[i] + 
                  to_vector(s * b_SxH)); // outputs I x H matrix, and then function vectorizes it in column-major order, matching i & h
}
generated quantities {
    
}
'

# stan_program <- '
# data {
#     int<lower=1> I; //# of individuals
#     int<lower=1> H; //# of hobbies
#     int<lower=1> S; //# of subjects
#     int<lower=0,upper=1> y[I * H]; //whether individual i partakes in hobby h
#     int<lower=1,upper=I> i[I * H]; //individual index for y
#     int<lower=1,upper=H> h[I * H]; //hobby index for y
#     vector<lower=0.0,upper=1.0>[S] s[I * H]; //subject indicators for y
#     vector<lower=0.0,upper=1.0>[I * H] t; //teacher indicator for y
# }
# transformed data {
#     
# }
# parameters {
#     //main parameters
#     real b_0; //avg hobby frequency
#     real b_T; //average teacher effect
#     vector[H] raw_b_H; //non-centered hobby prob
#     vector[H] raw_b_TxH; //non-centered teacher x hobby effect
#     vector[S] raw_b_S; //non-centered average subject effect
#     vector[H] raw_b_SxH[S]; //non-centered subject x hobby effect
#     
#     //hyperparameters
#     real<lower=0> sd_H;
#     real<lower=0> sd_TxH;
#     real<lower=0> sd_S;
#     real<lower=0> sd_SxH;
# }
# transformed parameters {
#     vector[H] b_H = raw_b_H * sd_H;
#     vector[H] b_TxH = raw_b_TxH * sd_TxH;
#     vector[S] b_S = raw_b_S * sd_S;
#     vector[H] b_SxH[S];
#     for(si in 1:S){
#       b_SxH[si,] = raw_b_SxH[si,] * sd_SxH;
#     }
# }
# model {
#     //priors
#     b_0 ~ normal(0,2);
#     b_T ~ std_normal();
#     raw_b_H ~ std_normal();
#     raw_b_TxH ~ std_normal();
#     raw_b_S ~ std_normal();
#     for(si in 1:S){
#       raw_b_SxH[si,] ~ std_normal();
#     }
#      
#     //hyperpriors
#     sd_H ~ std_normal();
#     sd_TxH ~ std_normal();
#     sd_S ~ std_normal();
#     sd_SxH ~ std_normal();
#     
#     //likelihood
#     vector[I * H] theta = rep_vector(0, I * H);
#     for(ind in 1:(I*H)){
#       theta[ind] = b_0 + b_T * t[ind] + b_H[h[ind]] + b_TxH[h[ind]] * t[ind] + sum(b_S .* s[ind,]) + sum(b_SxH[,h[ind]] .* s[ind,]);
#     }
#     y ~ bernoulli_logit(theta);
# }
# generated quantities {
#     
# }
# '


# 
# base = paste0("gelman_model_centered")
# stan_program <- '
# data {
#     int<lower=1> I; //# of individuals
#     int<lower=1> H; //# of hobbies
#     int<lower=1> S; //# of subjects
#     int<lower=0,upper=1> y[I * H]; //whether individual i partakes in hobby h
#     int<lower=1,upper=I> i[I * H]; //individual index for y
#     int<lower=1,upper=H> h[I * H]; //hobby index for y
#     matrix<lower=0.0,upper=1.0>[I * H, S] s; //subject indicators for y
#     vector<lower=0.0,upper=1.0>[I * H] t; //teacher indicator for y
# }
# transformed data {
#     
# }
# parameters {
#     //main parameters
#     real b_0; //avg hobby frequency
#     real b_T; //average teacher effect
#     vector[H] b_H; //centered hobby prob
#     vector[H] b_TxH; //centered teacher x hobby effect
#     vector[S] b_S; //centered average subject effect
#     matrix[S,H] b_SxH; //centered subject x hobby effect
#     
#     //hyperparameters
#     real<lower=0> sd_H;
#     real<lower=0> sd_TxH;
#     real<lower=0> sd_S;
#     real<lower=0> sd_SxH;
# }
# transformed parameters {
# 
# }
# model {
#     //priors
#     b_0 ~ normal(0,2);
#     b_T ~ std_normal();
#     b_H ~ normal(0, sd_H);
#     b_TxH ~ normal(0, sd_TxH);
#     b_S ~ normal(0, sd_S);
#     to_vector(b_SxH) ~ normal(0, sd_SxH)
#     //for(si in 1:S){
#     //  b_SxH[si,] ~ normal(0, sd_SxH);
#     //}
#      
#     //hyperpriors
#     sd_H ~ std_normal();;
#     sd_TxH ~ std_normal();
#     sd_S ~ std_normal();
#     sd_SxH ~ std_normal();
#     
#     //likelihood
#     y ~ bernoulli_logit(b_0 + 
#                   b_T * t + 
#                   b_H[h] + 
#                   b_TxH[h] .* t + 
#                   to_vector(s * b_S) + 
#                   to_vector(s * b_SxH)[h]);
# }
# generated quantities {
#     
# }
# '



#bob's model 1
# base = paste0("simple_test")
# 
# hits_in_row <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], rowhits[[ri]]))))
# hits_not_in_row <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], setdiff(1:n, rowhits[[ri]])))))
# y <- abind::abind(hits_in_row, hits_not_in_row, along = 3)
# N <- rbind(nrowhits, n - nrowhits)
# 
# data <- data.frame(do.call(rbind, lapply(1:ncol, function(ci) do.call(rbind, lapply(1:nrow, function(ri)
#   c(count = length(intersect(colhits[[ci]], rowhits[[ri]])),
#     total = n,
#     row_count = nrowhits[ri],
#     col_count = ncolhits[ci],
#     row_index = ri,
#     col_index = ci,
#     row_n = nrow,
#     col_n = ncol)))
# )))
# 
# #sample hits for each population
# d <- list(row_n = nrow,
#           col_n = ncol,
#           slice_n = 2)
# 
# stan_program <- "
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=1> slice_n;
#     
#     int<lower=0> N[slice_n, row_n]; \\number of teachers in school s and subject a
#     int<lower=0> y[col_n, row_n, slice_n]; \\number of teachers in school s and subject a who have hobby h, in 0:N[slice, row]
# 
# }
# parameters {
#     //col_n hobbies, slice_n schools, and row_n subject areas
#     alpha[1:col_n]; \\intercept
#     beta[1:2, 1:col_n];  \\school effects
#     gamma[1:row_n, 1:col_n];  \\subject effects
# }
# model {
#     //priors and hyperpriors
#     alpha ~ std_normal();
#     beta[s, 1:col_n] ~ multi_normal_cholesky(mu_beta, L_Sigma_beta);
#     gamma[a, 1:col_n] ~ multi_normal_cholesky(mu_gamma, L_Sigma_gamma);
#     
#     //likelihood
#     y[h, s, a] ~ binomial_logit(N[s, a], alpha + beta[s] + gamma[a]);
# }"
# 
# 
# bothhits <- sapply(1:ncol, function(ci) sapply(1:nrow, function(ri) length(intersect(colhits[[ci]], rowhits[[ri]]))))
# 
# data <- data.frame(do.call(rbind, lapply(1:ncol, function(ci) do.call(rbind, lapply(1:nrow, function(ri) 
#   c(count = length(intersect(colhits[[ci]], rowhits[[ri]])),
#     total = n,
#     row_count = nrowhits[ri],
#     col_count = ncolhits[ci],
#     row_index = ri,
#     col_index = ci,
#     row_n = nrow,
#     col_n = ncol)))
# )))
# 
# #sample hits for each population
# d <- list(cell_count = data$count,
#           total = data$total,
#           row_count = data$row_count,
#           col_count = data$col_count,
#           row_index = data$row_index,
#           col_index = data$col_index,
#           row_n = nrow,
#           col_n = ncol)

# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=0> total[row_n * col_n];
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> row_count[row_n * col_n];
#     int<lower=0> col_count[row_n * col_n];
#     int<lower=0> cell_count[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
#     int<lower=0> cell_count_compl[n];
#     int<lower=0> row_count_compl[n];
#     for(i in 1:n){
#       cell_count_compl[i] = col_count[i] - cell_count[i];
#       row_count_compl[i] = total[i] - row_count[i];
#     }
# }
# parameters {
#     //col params
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
#     real<lower=0> cell_sd;
#     vector[n] raw_cell_logodds;
# 
#     //biases in deviations terms
#     real overall_difference;
#     vector[row_n] raw_row_difference;
#     vector[col_n] raw_col_difference;
#     vector[n] raw_cell_difference;
#     real<lower=0> row_difference_sd;
#     real<lower=0> col_difference_sd;
#     real<lower=0> cell_difference_sd;
#     
#     //correlation params pre-multiply raw bias terms
#     cholesky_factor_corr[col_n] L_col_difference;
#     cholesky_factor_corr[row_n] L_row_difference;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd;
#     vector[n] cell_logodds = raw_cell_logodds * cell_sd + col_logodds[col_index];
# 
#     //incorporate bias
#     vector[col_n] col_difference = L_col_difference * raw_col_difference * col_difference_sd;
#     vector[row_n] row_difference = L_row_difference * raw_row_difference * row_difference_sd;
#     vector[n] cell_difference = raw_cell_difference * cell_difference_sd;
#     vector[n] cell_logodds_focal = cell_logodds +
#               (overall_difference + row_difference[row_index] + col_difference[col_index] + cell_difference) / 2;
#     vector[n] cell_logodds_compl = cell_logodds -
#               (overall_difference + row_difference[row_index] + col_difference[col_index] + cell_difference) / 2;
# }
# model {
#     //priors and hyperpriors
# 
#     //marginal params
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     raw_col_logodds ~ std_normal();
#     raw_cell_logodds ~ std_normal();
#     cell_sd ~ std_normal();
# 
#     //bias params
#     overall_difference ~ std_normal();
# 
#     raw_col_difference ~ std_normal();
#     col_difference_sd ~ std_normal();
# 
#     raw_row_difference ~ std_normal();
#     row_difference_sd ~ std_normal();
# 
#     raw_cell_difference ~ std_normal();
#     cell_difference_sd ~ std_normal();
#     
#     //bias correlation params
#     L_col_difference ~ lkj_corr_cholesky(1.0);
#     L_row_difference ~ lkj_corr_cholesky(1.0);
#     
#     //likelihood
#     cell_count ~ binomial_logit(row_count, cell_logodds_focal);
#     cell_count_compl ~ binomial_logit(row_count_compl, cell_logodds_compl);
# 
# }
# generated quantities {
#     vector[n] cell_total_prob_difference = inv_logit(cell_logodds_focal) - inv_logit(cell_logodds_compl);
# }
# '

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#write model
write_stan_file(stan_program, dir = "~/Desktop/", basename = base)
write_stan_json(d, paste0("~/Desktop/", base, ".json"))

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E2, iter_warmup = 2E2, data = d, parallel_chains = 4, 
                   adapt_delta = 0.85, refresh = 1, init = 0.1, max_treedepth = 15, thin = 2)
out <- mod$variational(data = d, output_samples = 1000)

check_mcmc_diagnostics <- F
if(check_mcmc_diagnostics){
  summ <- out$summary()
  summ[order(summ$ess_bulk),]
  summ[order(summ$rhat, decreasing = T),]
}

samps <- data.frame(as_draws_df(out$draws()))

#evaluate marginal posteriors of Gelman's model
par(mfrow = c(3,2))
CI_range <- 0:100/100
q.b_S <- apply(subset_samps("b_S", c("raw", "sd", "b_SxH"), samps = samps), 2, prop_greater_than_0)
hist(q.b_S, xlim = c(0,1), breaks = 0:20/20, freq = F, xlab = "quantile of true value in marginal posterior", main = "b_S")
CI_coverage.b_S <- sapply(CI_range, function(CIr) sum((sapply(q.b_S, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.b_S))
plot(CI_range, CI_coverage.b_S, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value (0)", main = "b_S")
abline(0,1,lty=2,lwd=2,col=2)
legend("topleft", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")

q.b_SxH <- apply(subset_samps("b_SxH", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
hist(q.b_SxH, xlim = c(0,1), breaks = 0:20/20, freq = F, xlab = "quantile of true value in marginal posterior", main = "b_SxH")
CI_coverage.b_SxH <- sapply(CI_range, function(CIr) sum((sapply(q.b_SxH, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.b_SxH))
plot(CI_range, CI_coverage.b_SxH, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value (0)", main = "b_SxH")
abline(0,1,lty=2,lwd=2,col=2)
legend("bottomright", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")

q.b_TxH <- apply(subset_samps("b_TxH", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
hist(q.b_TxH, xlim = c(0,1), breaks = 0:20/20, freq = F, xlab = "quantile of true value in marginal posterior", main = "b_TxH")
CI_coverage.b_TxH <- sapply(CI_range, function(CIr) sum((sapply(q.b_TxH, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.b_TxH))
plot(CI_range, CI_coverage.b_TxH, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value (0)", main = "b_TxH")
abline(0,1,lty=2,lwd=2,col=2)
legend("bottomright", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")

mean(samps$b_T > 0)

hist(samps$b_T)


dev.off()
hist(samps$overall_difference)
mean(samps$overall_difference < 0)

prop_greater_than_0 <- function(x) mean(x>0)

cell_difference <- apply(subset_samps("cell_difference", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
sum(cell_difference > 0.9)
sum(cell_difference < 0.1)
hist(cell_difference, breaks = 100)

cell_total_difference <- apply(subset_samps("cell_total_prob_difference", c("raw", "sd", "logodds"), samps = samps), 2, prop_greater_than_0)
sum(cell_total_difference > 0.9)
sum(cell_total_difference < 0.1)
hist(cell_total_difference, breaks = 10)

row_difference <- apply(subset_samps("row_difference", c("raw", "sd", "L"), samps = samps), 2, prop_greater_than_0)
sum(row_difference > 0.9)
sum(row_difference < 0.1)
hist(row_difference, breaks = 100)

col_difference <- apply(subset_samps("col_difference", c("raw", "sd", "L"), samps = samps) + samps$overall_difference, 2, prop_greater_than_0)
col_difference <- apply(subset_samps("col_difference", c("raw", "sd"), samps = samps), 2, prop_greater_than_0)
sum(col_difference > 0.9)
sum(col_difference < 0.1)
hist(col_difference, breaks = 100)

#reconstruct correlation matrices
#columns
L_col_difference_samps <- subset_samps("L_col_difference", c("raw"), samps = samps)
col_difference_corrs_samps <- do.call(abind::abind, list(lapply(1:nrow(samps), function(i){L <- matrix(unlist(L_col_difference_samps[i,]), ncol, ncol); L %*% t(L)}), along = 3))
col_difference_corrs_mean <- apply(col_difference_corrs_samps, c(1,2), mean)
col_difference_corrs_gr0.5 <- apply(col_difference_corrs_samps, c(1,2), prop_greater_than_0)
hist(col_difference_corrs_gr0.5[upper.tri(col_difference_corrs_gr0.5)])
plot(col_difference_corrs_mean[upper.tri(col_difference_corrs_mean)], rcol[upper.tri(rcol)])

#rows
L_row_difference_samps <- subset_samps("L_row_difference", c("raw"), samps = samps)
row_difference_corrs_samps <- do.call(abind::abind, list(lapply(1:nrow(samps), function(i){L <- matrix(unlist(L_row_difference_samps[i,]), nrow, nrow); L %*% t(L)}), along = 3))
row_difference_corrs_mean <- apply(row_difference_corrs_samps, c(1,2), mean)
row_difference_corrs_gr0.5 <- apply(row_difference_corrs_samps, c(1,2), prop_greater_than_0)
hist(row_difference_corrs_gr0.5[upper.tri(row_difference_corrs_gr0.5)])

