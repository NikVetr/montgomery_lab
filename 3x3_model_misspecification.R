library(cmdstanr)
library(posterior)

#overall parameters
a <- 0
sigma_y <- 2
sigma_eff <- 2
n <- 10
n_cat <- 5
x <- rep(1:n_cat, each = n / n_cat)
n <- length(x)

###################
######## A ########
###################

# Parameters and simulation for Model A
model_a <- list(
  a = a,
  sigma_y = sigma_y,
  sigma_beta = sigma_eff,
  n = n,
  n_cat = n_cat
)

# Simulate x
model_a$x <- x

# Simulate b values independently for each state of x
model_a$b <- rnorm(n_cat, 0, model_a$sigma_beta)

# Compute mu for each observation
model_a$mu <- model_a$a + model_a$b[model_a$x]

# Simulate y
model_a$y <- rnorm(model_a$n, model_a$mu, model_a$sigma_y)

###################
######## B ########
###################

# Parameters and simulation for Model B
model_b <- list(
  a = a,
  sigma_y = sigma_y,
  sigma_e = sigma_eff,
  n = n,
  n_cat = n_cat
)

# Use x from Model A (assuming the same x for all models)
model_b$x <- model_a$x 

# Simulate e values and compute b values cumulatively
model_b$e <- rnorm(n_cat-1, 0, model_b$sigma_e)
model_b$b <- cumsum(c(0, model_b$e)) # Assume b[1] = 0 for simplicity, adjust as needed

# Compute mu for each observation
model_b$mu <- model_b$a + model_b$b[model_b$x]

# Simulate y
model_b$y <- rnorm(model_b$n, model_b$mu, model_b$sigma_y)

###################
######## C ########
###################

# Parameters and simulation for Model C
model_c <- list(
  a = a,
  sigma_y = sigma_y,
  sigma_e = sigma_eff,
  mu_e = 2,
  n = n,
  n_cat = n_cat
)

# Use x from Model A/B
model_c$x <- model_a$x 

# Simulate e values with pooled mean
model_c$e <- rnorm(n_cat-1, model_c$mu_e, model_c$sigma_e)
model_c$b <- cumsum(c(0, model_c$e))

# Compute mu for each observation
model_c$mu <- model_c$a + model_c$b[model_c$x]

# Simulate y
model_c$y <- rnorm(model_c$n, model_c$mu, model_c$sigma_y)


#package into list
data <- list(A = model_a, B = model_b, C = model_c)


stan_models <- list(
  
  A =
    
  "data {
    int<lower=0> N;         // Number of observations  
    int<lower=0> N_CAT;         // Number of categories  
    vector[N] y;            // Response variable
    array[N] int<lower=1, upper=N_CAT> x; // Predictor variable (ordinal, N_CAT states)
  }
  parameters {
    real a;                 // Intercept now a parameter
    real<lower=0> sigma_y;  // Standard deviation of y
    vector[N_CAT] b;            // Slope for each state of x
    real<lower=0> sigma_b;  // Standard deviation of b
  }
  model {
    b ~ normal(0, sigma_b);
    for (i in 1:N)
      y[i] ~ normal(a + b[x[i]], sigma_y);
  }
  ", 
  
  B =
    
  "data {
    int<lower=0> N;         // Number of observations  
    int<lower=0> N_CAT;         // Number of categories  
    vector[N] y;            // Response variable
    array[N] int<lower=1, upper=N_CAT> x; // Predictor variable (ordinal, N_CAT states)
  }
  parameters {
    real a;                 // Intercept now a parameter
    real<lower=0> sigma_y;
    vector[N_CAT-1] e;            // Differences between consecutive b values
    real<lower=0> sigma_e;  // Standard deviation of e
  }
  transformed parameters {
    vector[N_CAT] b;            // Cumulative sum to obtain b values
    b[1] = 0;               // Set the first value of b to 0 for identification
    for (i in 2:N_CAT)
      b[i] = b[i-1] + e[i-1];
  }
  model {
    e ~ normal(0, sigma_e);
    for (i in 1:N)
      y[i] ~ normal(a + b[x[i]], sigma_y);
  }
  ",
  
  C =
    
  "data {
    int<lower=0> N;         // Number of observations  
    int<lower=0> N_CAT;         // Number of categories  
    vector[N] y;            // Response variable
    array[N] int<lower=1, upper=N_CAT> x; // Predictor variable (ordinal, N_CAT states)
  }
  parameters {
    real a;                 // Intercept now a parameter
    real<lower=0> sigma_y;
    vector[N_CAT-1] e;
    real<lower=0> sigma_e;
    real mu_e;              // Mean of e for pooling
  }
  transformed parameters {
    vector[N_CAT] b;
    b[1] = 0;
    for (i in 2:N_CAT)
      b[i] = b[i-1] + e[i-1];
  }
  model {
    e ~ normal(mu_e, sigma_e);
    for (i in 1:N)
      y[i] ~ normal(a + b[x[i]], sigma_y);
  }
  "
)


# compile models
mod <- list(A = cmdstan_model(write_stan_file(stan_models$A), cpp_options = list(STAN_THREADS = T, stan_threads = TRUE)),
            B = cmdstan_model(write_stan_file(stan_models$B), cpp_options = list(STAN_THREADS = T, stan_threads = TRUE)),
            C = cmdstan_model(write_stan_file(stan_models$C), cpp_options = list(STAN_THREADS = T, stan_threads = TRUE))
            )

# Prepare data for Model A (using simulated data from above)
dat <- list(A = list(N = data$A$n, N_CAT = data$A$n_cat, y = data$A$y, x = as.integer(data$A$x)),
            B = list(N = data$B$n, N_CAT = data$B$n_cat, y = data$B$y, x = as.integer(data$B$x)),
            C = list(N = data$C$n, N_CAT = data$C$n_cat, y = data$C$y, x = as.integer(data$C$x)))


# Fit models to each others' data
fits <- list()
for(sim_model in c("A", "B", "C")){
  fits[[sim_model]] = list()
  for(inf_model in c("A", "B", "C")){
    fits[[sim_model]][[inf_model]]$fit <- mod[[inf_model]]$sample(data = dat[[sim_model]], adapt_delta = 0.95, 
                                                                  refresh = 500, max_treedepth = 15, 
                                                                  threads_per_chain = 1, parallel_chains = 4)
    fits[[sim_model]][[inf_model]]$samps <- data.frame(as_draws_df(fits[[sim_model]][[inf_model]]$fit))
  }  
}

#### plot 3 x 3 grid? or 3 x 1 grid ####
model_names <- setNames(c("A", "B", "C"), c("A", "B", "C"))
b_range <- lapply(model_names, function(sim_model){
  range(unlist(sapply(c("A", "B", "C"), function(inf_model){
    
    rout <- range(c(apply(fits[[sim_model]][[inf_model]]$samps[,paste0("b.", 1:n_cat, ".")] + 
                          fits[[sim_model]][[inf_model]]$samps$a, 
                        2, quantile, prob = c(0.04, 0.96))), data[[sim_model]]$b + data[[sim_model]]$a)
    return(rout)
    
  })))
})


model_key <- c(A = "B$_{i}$ ~ normal(0, \\sigma$_B$)", 
               B = "B$_{i}$ = B$_{i-1}$ + e$_{i}$, e$_{\\mu}$ = 0", 
               C = "B$_{i}$ = B$_{i-1}$ + e$_{i}$, e$_{\\mu}$ ~ normal(0, \\sigma$_e$)")
model_cols <- c(A = "darkred", B = "purple", C = "darkblue")
model_disp <- setNames(c(0:2) / 6, names(model_cols))

grDevices::cairo_pdf(filename = paste0("~/repos/egtex-ase/figures/3x3 comparison.pdf"), 
                     width = 1250 / 72, height = 600 / 72, family="Arial Unicode MS", pointsize = 19)
par(mfrow = c(1,4))


for(sim_model in c("A", "B", "C")){
  
  for(inf_model in c("A", "B", "C")){
    
    focal_samps <- fits[[sim_model]][[inf_model]]$samps[,paste0("b.", 1:n_cat, ".")] + fits[[sim_model]][[inf_model]]$samps$a
    CE90 <- apply(focal_samps, 2, quantile, prob = c(0.05, 0.95))
    true_vals <- data[[sim_model]]$b + data[[sim_model]]$a
    
    if(inf_model == "A"){
      plot(true_vals, 1:n_cat, xlim = b_range[[sim_model]], ylim = c(0.5, n_cat+0.5), 
           ylab = "effect category", yaxt = "n", xlab = "category EV",
           main = "", pch = 19, col = model_cols[sim_model], cex = 1.5)  
      sim_text <- latex2exp::TeX(paste0("Sim. Model: ", model_key[sim_model], ""))
      text(sim_text, pos = 3, col = model_cols[sim_model], y = par("usr")[4], 
           x = mean(par("usr")[1:2]), xpd = NA, cex = 1.25)
      axis(2, 1:n_cat, 1:n_cat)
      segments(x0 = true_vals, x1 = apply(CE90, 2, mean), y0 = 1:n_cat, lwd = 2,
               y1 = 1:n_cat - 0.25 - model_disp[inf_model], col = model_cols[sim_model])
      
    }
    
    
    
    # inf_text <- latex2exp::TeX(paste0("Inf. Model: ", model_key[inf_model]))
    # text(inf_text, pos = 3, col = model_cols[inf_model], y = par("usr")[4] + strheight(sim_text), 
    #      x = mean(par("usr")[1:2]), xpd = NA)
          
    
    segments(x0 = CE90[1,], x1 = CE90[2,], y0 = 1:n_cat - 0.25 - model_disp[inf_model], 
             y1 = 1:n_cat - 0.25 - model_disp[inf_model], col = model_cols[inf_model], lwd = 3)
    
    
        
  }
  
}

legend(x = 12, y = par("usr")[4], lwd = c(NA, 3), pch = c(19, NA), 
       legend = c("True Value", "90% CI"), xpd = NA, box.lty = 2, pt.cex = c(1.5,NA), seg.len = 2)

dev.off()