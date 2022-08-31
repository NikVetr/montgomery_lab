#specify a model with multiple parameters
n = 20
a = 2
b = 2.5
c = 1.5

#binomial x1
x1_prob <- 0.23
x1 = rbinom(n = n, prob = x1_prob, size = 2)

#x2 is collinear with x1
x2 = rnorm(n) + x1 * rnorm(n, x1)

#y is the outcome variable, so simulate y
y = a + b * x1 + c * x2 + rnorm(n, sd = 3)
d <- data.frame(x1 = x1, x2 = x2, y = y)
d$ys <- (d$y - mean(d$y)) / sd(d$y)

#fit the linear model the usual way
fit = lm(ys~x1+x2, d = d)

#specify the likelihood function
func_string <- sapply(1:n, function(i) paste0("(1 / sig / sqrt(2 * pi) * exp(-1/2 * ((", 
                                              d$ys[i], "-(a+b*", d$x1[i], "+ c*", d$x2[i], "))/sig)^2))"))
func <- parse(text = paste0(func_string, collapse = " * "))
vars <- c(sig = summary(fit)$sigma, 
          a =  summary(fit)$coefficients["(Intercept)","Estimate"], 
          b =  summary(fit)$coefficients["x1","Estimate"],
          c =  summary(fit)$coefficients["x2","Estimate"])

#compute the hessian at the MLE 
#(can also get numerically from stats::optim(hessian=T) or optimx::optimx(hessian=T))
#cos calculus::hessian() uses a symbolic calclulus library iirc
hess_mat <- hessian(f = func, var = vars)

#compare the standard errors you got the normal way (ie analytically, via ols) 
#to the inverse of the negative rescaled hessian
summary(fit)$coefficients[,"Std. Error"]
sqrt(diag(solve(-hess_mat / prod(dnorm(x = d$ys, 
                                       mean = vars["a"] + vars["b"] * d$x1 + vars["c"] * d$x2,
                                       sd = vars["sig"])))))[-1]

#bada bing bada boom you can now use these for hypothesis testing same as you would anything else
summary(fit)$coefficients[,"Pr(>|t|)"]
1 - abs(0.5 - pt(q = summary(fit)$coefficients[,"Estimate"] / summary(fit)$coefficients[,"Std. Error"], df = n-3)) * 2

