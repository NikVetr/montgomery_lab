
simfit <- function(n = 1E1, a1 = 3, b1= 2, a2 = 5, b2 = 6, s1 = 2, s2 = 3, prop_s_differ = 0.3, prop_zero_b = 0.1){
  x1 = rnorm(n)
  x2 = rnorm(n)
  s_differ <- sample(c(T,F), 1, prob = c(prop_s_differ, 1-prop_s_differ))
  zero_b <- sample(c(T,F), 1, prob = c(prop_zero_b, 1-prop_zero_b))
  if(zero_b){b1 <- 0; b2 <- 0}
  y1 = a1 + b1 * x1 + rnorm(n, sd = s1)
  y2 = a2 + b2 * x2 + rnorm(n, sd = ifelse(s_differ, s1, s2))
  d <- data.frame(x = c(x1, x2), y = c(y1, y2), sex = c(rep(0, n), rep(1, n)))
  m_inter <- lm(y ~ 1 + x + sex + x*sex, data = d)
  # m_sex0 <- lm(y ~ 1 + x, data = d[d$sex == 0,])
  # m_sex1 <- lm(y ~ 1 + x, data = d[d$sex == 1,])
  # return(c(summary(m_inter)$coefficients[2,4], summary(m_sex0)$coefficients[2,4], summary(m_sex1)$coefficients[2,4]))
  vt <- var.test(m_inter$residuals[1:n], m_inter$residuals[(n+1):(2*n)])
  return(vt$p.value)
}

n_rep <- 4E3
prop_s_differ = 0.4
pvals <- t(replicate(n_rep, simfit(n = 1E3, prop_zero_b = 1, prop_s_differ = prop_s_differ, b1 = 2, b2 = 2)))
1-qvalue::pi0est(pvals)$pi0
