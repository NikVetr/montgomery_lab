n <- 1E2
g <- 1000
gs <- 1:g
p0 <- 0.1
nrep <- 1E4

sim_drift <- function(n, g, p0){
  ps <- rep(NA, g)
  ps[1] <- p0
  
  x <- rbinom(n, 2, ps[1])
  for(i in 2:g){
    ps[i] <- sum(x) / n / 2
    x <- rbinom(n, 2, ps[i])
  }
  
  return(ps)
}

ps <- sim_drift(n, g, p0)
plot(ps, ylim = c(0,1), type = "l")

# pmat <- replicate(nrep, sim_drift(n, g, p0))
pmat <- do.call(cbind, parallel::mclapply(1:nrep, function(i) sim_drift(n, g, p0), mc.cores = 12))
pvar <- apply(pmat, 1, var)
pdmat <- apply(pmat, 2, diff)
pdvar <- apply(pdmat, 1, var)

par(mfrow = c(2,1), mar = c(4,4,1,2))
plot(gs, pvar, type = "l", ylim = range(c(0, p0 * (1-p0))))
abline(h = p0 * (1-p0), lty = 2, col = adjustcolor(2,0.5), lwd = 2)
plot(gs[-1], pdvar, type = "l")
points(0, p0 * (1-p0) / 2 / n, pch = 17, col = 2, cex = 2)

