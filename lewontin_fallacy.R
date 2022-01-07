#define some functions
var_mix_2norm <- function(p1, p2, m1, m2, v1, v2) p1*v1 + p2*v2 + p1*m1^2 + p2*m2^2 - (p1*m1 + p2*m2)^2
dist_vec <- function(x, inds, pop2) sqrt(sum((x[inds[1],,1] - x[inds[2],,pop2])^2))

#check to make sure formula works lol
n = 5E5
var(c(rnorm(n), rnorm(n / 3, sd = 2) + 1.25))
var_mix_2norm(0.75, 0.25, 0, 1.25, 1, 4)

#specify simulation parameters
n_dim <- 2E4
n_obs <- 1E3
mean_separation <- 0.45

#simulate data
x <- abind::abind(t(matrix(rnorm(n_obs*n_dim),n_dim, n_obs)), t(matrix(rnorm(n_obs*n_dim),n_dim, n_obs)) + mean_separation, along = 3)

#find distances
n_rep <- 5E3
within_dists <- replicate(1E3, dist_vec(x, sample(1:n_obs, 2, F), 1))
between_dists <- replicate(1E3, dist_vec(x, sample(1:n_obs, 2, F), 2))
difference_dists <- between_dists - within_dists

#find prop total variance
prop_variance_within_groups <- 1 / var_mix_2norm(0.5, 0.5, 0, mean_separation, 1, 1)

#do plotting
breaks <- seq(from = min(c(within_dists, between_dists)), to = max(c(within_dists, between_dists)), length.out = 30)
hist_w <- hist(within_dists, breaks = breaks, plot = F)
hist_b <- hist(between_dists, breaks = breaks, plot = F)
par(mfrow = c(1,2), mar = c(5,5,4,3.5))
plot(hist_w, col = adjustcolor("orange", 0.5), ylim = c(0,max(c(hist_w$density, hist_b$density))), freq = F,
     main = paste0(round(prop_variance_within_groups*100, 2), "% of variance within groups"),
     xlab = "distances between individuals")
hist(between_dists, col = adjustcolor("blue", 0.35), breaks = breaks, add = T, freq = F)
legend(x = mean(par("usr")[1:2]) - diff(par("usr")[1:2])/5, y = par("usr")[4], cex = 1,
       col = c(adjustcolor("orange", 0.5), adjustcolor("blue", 0.35)), bg = adjustcolor(1, 0.02),
       pch = 15, legend = c("within-groups", "between-groups"))
hist(difference_dists, xlim = c(0, max(difference_dists)), freq = F, xlab = "difference in distances",
     main = "distribution of differences between distances", 
     breaks = seq(0, max(difference_dists), length.out = 30))


