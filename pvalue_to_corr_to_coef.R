#specify generative model parameters
n = 20000
a = 3
b = 0.001

#simulate data
x = rnorm(n)
y = a + b * x + rnorm(n)
d <- data.frame(x = x, y = y)

d$xs <- (d$x - mean(d$x)) / sd(d$x)
d$ys <- (d$y - mean(d$y)) / sd(d$y)

#from correlation to pval
summary(lm(ys~xs, d))$coefficients[2,4]
(1 - pnorm(abs(atanh(cor(d$x,d$y))) / (1 / sqrt(n-3))))*2

#from pval to correlation
pval <- summary(lm(y~x, d))$coefficients[2,4]
tanh(qnorm(1-pval/2) / sqrt(n-3)) * sign(summary(lm(y~x, d))$coefficients[2,1])
cor(d$x,d$y)

#from pval to partly standardized coefficient
x_prob <- 0.23
x = rbinom(n = n, prob = x_prob, size = 2)
y = a + b * x + rnorm(n)
d <- data.frame(x = x, y = y)
d$ys <- (d$y - mean(d$y)) / sd(d$y)
d$xs <- (d$x - mean(d$x)) / sd(d$x)

pval <- summary(lm(y~x, d))$coefficients[2,4]

summary(lm(ys~x, d))$coefficients[2,1]
tanh(qnorm(1-pval/2) / sqrt(n-3)) / sqrt(x_prob * (1-x_prob) * 2) * sign(summary(lm(ys~x, d))$coefficients[2,1])

#can point estimates and p-vals / SEs be used to reconstruct a likelihood surface?
n = 500
a = 2
b = 2.5

x = rnorm(n)
y = a + b * x + rnorm(n, sd = 3)

d <- data.frame(x = x, y = y)
d$xs <- (d$x - mean(d$x)) / sd(d$x)
d$ys <- (d$y - mean(d$y)) / sd(d$y)

summary(lm(y~x, d))
summary(lm(ys~xs, d))
var(summary(lm(ys~xs, d))$residuals)

# likelihood_surface_abs <- expand.grid(a = -100:100/10, b = -100:100/10, s = 0:100/10)
b_range <- -1000:1000/1000
likelihood_surface_b <- data.frame(b = b_range, ll = sapply(b_range, function(bi) sum(dnorm(x = d$ys, mean = 0 + d$xs*bi, sd = 1, log = T))))
# likelihood_surface_b$ll <- likelihood_surface_b$ll - max(likelihood_surface_b$ll) + 3
likelihood_surface_b$b[which.max(likelihood_surface_b$ll)]
plot(likelihood_surface_b$b, exp(likelihood_surface_b$ll), type = "l")
plot(likelihood_surface_b$b, (likelihood_surface_b$ll), type = "l")

moments_from_grid <- function(x, d, m = c(1:2), ll = F){
  
  #prevent underflow
  if(ll){
    d <- d - max(d) + 10
    d <- exp(d)
  } 
  
  m1 <- sum(x*d) / sum(d)
  m_eval <- setdiff(m, 1)
  ms <- sapply(m_eval, function(m_i) sum((x-m1)^m_i*d) / sum(d))
  if(1 %in% m){return(unlist(c(m1, ms)))}else{return(ms)}  
  
}

# moments_from_grid(-1000:1000/100, dnorm(-1000:1000/100, mean = 2, sd = 2), m = 1:5)
# moments_from_grid(-1000:1000/100, dnorm(-1000:1000/100, mean = 2, sd = 2, log = T), m = 1:5, ll = T)

c(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 1, T), 
  sqrt(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 2, T)))

summary(lm(ys~xs, d))$coefficients[2,1:2]

run_sim <- function(){
  
  n = 500
  a = runif(1, -2, 2)
  b = runif(1, -2, 2)
  
  x_prob <- 0.23
  x = rbinom(n = n, prob = x_prob, size = 2)  
  y = a + b * x + rnorm(n, sd = 3)
  
  d <- data.frame(x = x, y = y)
  d$xs <- (d$x - mean(d$x)) / sd(d$x)
  d$ys <- (d$y - mean(d$y)) / sd(d$y)  
  
  b_range <- -1000:1000/1000
  likelihood_surface_b <- data.frame(b = b_range, ll = sapply(b_range, function(bi) sum(dnorm(x = d$ys, mean = 0 + d$xs*bi, sd = 1, log = T))))
  
  c(summary(lm(ys~xs, d))$coefficients[2,2], sqrt(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 2, T)))
}

SEs_vs_LL <- t(replicate(1E2, run_sim()))
SEs_vs_LL <- SEs_vs_LL[order(SEs_vs_LL[,1]),]
plot(SEs_vs_LL, type = "l")

#### using the hessian? ####
library(calculus)

# simulate some data
n = 20
a = 2
b = 2.5

x = rnorm(n)
y = a + b * x + rnorm(n, sd = 3)
d <- data.frame(x = x, y = y)

#fit a linear model
fit = lm(y~x, d = d)


#specify the likelihood
func_string <- sapply(1:n, function(i) paste0("(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                              y[i], "-(a+b*", x[i], "))/sig)^2))"))
func <- parse(text = paste0(func_string, collapse = " * "))
vars <- c(sig = summary(fit)$sigma, 
                a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
                b =  summary(fit)$coefficients["x","Estimate"])

#compute the hessian
hess_mat <- hessian(f = func, var = vars)

#compare the standard errors compute the normal way (ie analytically) 
#to the inverse of the negative rescaled hessian
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat / prod(dnorm(x = y, 
                                       mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * x,
                                       sd = summary(fit)$sigma)))))[-1]

jacobian(f = func, var = vars)

-solve(diag(summary(fit)$coefficients[,"Std. Error"])^2) * prod(dnorm(x = y, 
                                                                      mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * x,
                                                                      sd = summary(fit)$sigma))
hess_mat[-1,-1]


#standardize output variable
# simulate some data
n = 20
a = 2
b = 2.5

x = rnorm(n)
y = a + b * x + rnorm(n, sd = 3)
d <- data.frame(x = x, y = y)
d$xs <- (d$x - mean(d$x)) / sd(d$x)
d$ys <- (d$y - mean(d$y)) / sd(d$y)


###ohhh wait let's go back to the moments thing
#fit a linear model
fit = lm(ys~x, d = d)

#specify the likelihood
func_string <- sapply(1:n, function(i) paste0("(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                              d$ys[i], "-(a+b*", d$x[i], "))/sig)^2))"))
func <- parse(text = paste0(func_string, collapse = " * "))
vars <- c(sig = summary(fit)$sigma, 
          a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
          b =  summary(fit)$coefficients["x","Estimate"])

#compute the hessian
hess_mat <- hessian(f = func, var = vars)

#compare the standard errors compute the normal way (ie analytically) 
#to the inverse of the negative rescaled hessian
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat / prod(dnorm(x = d$ys, 
                                       mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * d$x,
                                       sd = summary(fit)$sigma)))))[-1]

jacobian(f = func, var = vars)

-solve(diag(summary(fit)$coefficients[,"Std. Error"])^2) * prod(dnorm(x = d$ys, 
                                                                      mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * d$x,
                                                                      sd = summary(fit)$sigma))
hess_mat[-1,-1]

b_range <- -500:500/100
likelihood_surface_b <- data.frame(b = b_range, ll = sapply(b_range, function(bi) sum(dnorm(x = d$ys, mean = 0 + d$x*bi, sd = summary(fit)$sigma, log = T))))
# likelihood_surface_b$ll <- likelihood_surface_b$ll - max(likelihood_surface_b$ll) + 3
likelihood_surface_b$b[which.max(likelihood_surface_b$ll)]
plot(likelihood_surface_b$b, exp(likelihood_surface_b$ll), type = "l")
c(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 1, T), 
  sqrt(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 2, T)))
summary(fit)$coefficients[2,1:2]


#### now with multiple predictors?
n = 20
a = 2
b = 2.5
c = 1.5

x1_prob <- 0.23
x1 = rbinom(n = n, prob = x1_prob, size = 2)
x2 = rnorm(n) + x1 * rnorm(n, x1)
y = a + b * x1 + c * x2 + rnorm(n, sd = 3)
d <- data.frame(x1 = x1, x2 = x2, y = y)
d$ys <- (d$y - mean(d$y)) / sd(d$y)
#fit a linear model
fit = lm(ys~x1+x2, d = d)

#specify the likelihood
func_string <- sapply(1:n, function(i) paste0("(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                              d$ys[i], "-(a+b*", d$x1[i], "+ c*", d$x2[i], "))/sig)^2))"))
func <- parse(text = paste0(func_string, collapse = " * "))
vars <- c(sig = summary(fit)$sigma, 
          a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
          b =  summary(fit)$coefficients["x1","Estimate"],
          c =  summary(fit)$coefficients["x2","Estimate"])

#compute the hessian
hess_mat <- hessian(f = func, var = vars)

#compare the standard errors compute the normal way (ie analytically) 
#to the inverse of the negative rescaled hessian
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat / prod(dnorm(x = d$ys, 
                                       mean = vars["a"] + vars["b"] * d$x1 + vars["c"] * d$x2,
                                       sd = vars["sig"])))))[-1]




jacobian(f = func, var = vars)

-solve(diag(summary(fit)$coefficients[,"Std. Error"])^2) * prod(dnorm(x = d$ys, 
                                                                      mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * d$x,
                                                                      sd = summary(fit)$sigma))
hess_mat[-1,-1]

b_range <- -500:500/100
likelihood_surface_b <- data.frame(b = b_range, ll = sapply(b_range, function(bi) sum(dnorm(x = d$ys, mean = 0 + d$x*bi, sd = summary(fit)$sigma, log = T))))
# likelihood_surface_b$ll <- likelihood_surface_b$ll - max(likelihood_surface_b$ll) + 3
likelihood_surface_b$b[which.max(likelihood_surface_b$ll)]
plot(likelihood_surface_b$b, exp(likelihood_surface_b$ll), type = "l")
c(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 1, T), 
  sqrt(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 2, T)))
summary(fit)$coefficients[2,1:2]


#### ok, let's bring it all together ####

## specify generative model parameters
n = 20
a = 3
b = 0.1

## simulate data
x_prob <- 0.23
x = rbinom(n = n, prob = x_prob, size = 2)
y = a + b * x + rnorm(n)
d <- data.frame(x = x, y = y)
d$ys <- (d$y - mean(d$y)) / sd(d$y)
d$xs <- (d$x - mean(d$x)) / sd(d$x)
fit <- lm(y~x, d)

## from correlation to pval
summary(fit)$coefficients[2,4]
(1 - pnorm(abs(atanh(cor(d$x,d$y))) / (1 / sqrt(n-3))))*2 #using Fisher transf. / approx.
tail_prop <- function(p) min(p, 1-p) * 2
tail_prop(pt(q = cor(d$x,d$y) * sqrt(n-2) / sqrt(1-cor(d$x,d$y)^2), df = n - 2))
cor.test(x, y)$p.value

## from pval to correlation
pval <- summary(fit)$coefficients[2,4]
tanh(qnorm(1-pval/2) / sqrt(n-3)) * sign(summary(lm(y~x, d))$coefficients[2,1]) #using Fisher transf. / approx.
inv_tail_prop_right <- function(p) 1 - p / 2

a = qt(p = inv_tail_prop_right(pval), df = n - 2) / sqrt(n-2)
abs(a / sqrt(a^2 +1)) * sign(summary(fit)$coefficients["x", "Estimate"])
cor(d$x,d$y)

## from pval to partly standardized coefficient
pval <- summary(fit)$coefficients[2,4]
tanh(qnorm(1-pval/2) / sqrt(n-3)) / sqrt(x_prob * (1-x_prob) * 2) * sign(summary(lm(ys~x, d))$coefficients[2,1]) #using Fisher transf. / approx.
summary(lm(ys~x, d))$coefficients[2,1]
a = qt(p = inv_tail_prop_right(pval), df = n - 2) / sqrt(n-2)
abs(a / sqrt(a^2 +1)) * sign(summary(fit)$coefficients["x", "Estimate"]) / sqrt(x_prob * (1-x_prob) * 2) 
#last bit is the variance of a binomial random variable, which would assume pop at HW-equilibrium. For exact equality to lm(), divide by sd(x), which we don't know
abs(a / sqrt(a^2 +1)) * sign(summary(fit)$coefficients["x", "Estimate"]) / sd(x)

## can use standard errors for quadratic approx of marginal posterior
library(calculus)

#specify the likelihood
likelihood_string <- sapply(1:n, function(i) paste0("(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                              d$ys[i], "-(a+b*", d$x[i], "))/sig)^2))"))
likelihood <- parse(text = paste0(likelihood_string, collapse = " * "))
vars <- c(sig = summary(fit)$sigma, 
          a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
          b =  summary(fit)$coefficients["x","Estimate"])

#compute & compare the hessian
hess_mat <- hessian(f = likelihood, var = vars)
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat / prod(dnorm(x = d$ys, 
                                       mean = summary(fit)$coefficients[1,1] + summary(fit)$coefficients[2,1] * d$x,
                                       sd = summary(fit)$sigma)))))[-1]

#now try the log-likelihood
loglikelihood_string <- sapply(1:n, function(i) paste0("log(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                                    d$ys[i], "-(a+b*", d$x[i], "))/sig)^2))"))
loglikelihood <- parse(text = paste0(loglikelihood_string, collapse = " + "))
vars <- c(sig = summary(fit)$sigma, 
          a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
          b =  summary(fit)$coefficients["x","Estimate"])
hess_mat <- hessian(f = loglikelihood, var = vars)
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat)))[-1]


## manually find first and second moment of marginal posterior
moments_from_grid <- function(x, d, m = c(1:2), ll = F){
  
  #prevent underflow -- gets normalized out later
  if(ll){
    d <- d - max(d) + 10
    d <- exp(d)
  } 
  
  m1 <- sum(x*d) / sum(d)
  m_eval <- setdiff(m, 1)
  ms <- sapply(m_eval, function(m_i) sum((x-m1)^m_i*d) / sum(d))
  if(1 %in% m){return(unlist(c(m1, ms)))}else{return(ms)}  
  
}

b_range <- -500:500/100
likelihood_surface_b <- data.frame(b = b_range, ll = sapply(b_range, function(bi) sum(dnorm(x = d$ys, mean = 0 + d$xs*bi, sd = summary(lm(ys~xs, d))$sigma, log = T))))
# likelihood_surface_b$ll <- likelihood_surface_b$ll - max(likelihood_surface_b$ll) + 3
plot(likelihood_surface_b$b, exp(likelihood_surface_b$ll), type = "l")
c(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 1, T), 
  sqrt(moments_from_grid(likelihood_surface_b$b, likelihood_surface_b$ll, 2, T)))
summary(lm(ys~xs, d))$coefficients[2,1:2]

# get residual sd
(1 - summary(lm(ys~0+xs, d))$sigma)*2
