#check stepwise regression with random noise taking out random noise each time -- residualizing noise each time, see deflated p-value

simulate_and_fit <- function(){
  n <- 1E2
  x = rbinom(n, 2, prob = 0.5)
  b = 0
  e <- rnorm(n,0,2)
  a <- 0
  y <- a + b * x + e 
  summary(lm(y ~ x))$coefficients[2,4] < 0.05
}

mean(replicate(1000, simulate_and_fit()))

residualize_noise <- function(y, x, n){
  noise <- rbinom(n, 2, prob = 0.5)
  fit <- lm(y ~ 1 + noise + x)
  y - fit$coefficients["noise"]*noise
}

recursive_residualize_noise <- function(y,x,n,nr){
  if(nr > 1){
    recursive_residualize_noise(residualize_noise(y,x,n),x,n,nr-1)
  } else {
    residualize_noise(y,x,n)
  }
}

simulate_and_fit_residualize <- function(){
  n <- 1E2
  x = rbinom(n, 2, prob = 0.5)
  b = 0
  e <- rnorm(n,0,2)
  a <- 0
  y <- a + b * x + e 
  y_resid <- recursive_residualize_noise(y,x,n,50)
  # summary(lm(y_resid ~ x))$coefficients[2,4] < 0.05
  # (1 - pt(summary(lm(y_resid ~ x))$coefficients[2,3], df = n-2)) * 2 < 0.05 #p-val from t-dist, no correction
  (1 - pt(summary(lm(y_resid ~ x))$coefficients[2,3], df = n-2-50)) * 2 < 0.05 #p-val from t-dist, with correction
}

mean(replicate(500, simulate_and_fit_residualize()))

#interpretation -- initial noise gets protected and eventually there ends up very little residual error around it, leading to inflated FPR


#now try an even simpler model -- normal random variable with mean 0

simulate_and_fit <- function(){
  n <- 1E2
  x = rnorm(n)
  summary(lm(x ~ 1))$coefficients[1,4] < 0.05
}

mean(replicate(1000, simulate_and_fit()))

residualize_noise <- function(x, n){
  noise <- rbinom(n, 1, prob = 0.5)
  fit <- lm(x ~ 1 + noise)
  x - fit$coefficients["noise"]*noise
}

recursive_residualize_noise <- function(x,n,nr){
  if(nr > 1){
    recursive_residualize_noise(residualize_noise(x,n),n,nr-1)
  } else {
    residualize_noise(x,n)
  }
}

simulate_and_fit_residualize <- function(){
  n <- 1E2
  x = rnorm(n)
  x_resid <- recursive_residualize_noise(x,n,50)
  summary(lm(x_resid ~ 1))$coefficients[1,4] < 0.05
}

mean(replicate(500, simulate_and_fit_residualize()))

#interpretation -- initial noise gets protected and eventually there ends up very little residual error around it, leading to inflated FPR