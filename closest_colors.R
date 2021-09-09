library(farver)
library(MotrpacBicQC)

#functions
polar2cart <- function(t, r){c(r*cos(t), r * sin(t))}
logit <- function(p){log(p/(1-p))}
invlogit <- function(x){exp(x)/(1+exp(x))}
cyl2cart <- function(t,r,z){c(polar2cart(t,r),z)}
cart2polar <- function(x,y){c(atan(y/x), sqrt(x^2 + y^2))}
dist_kd <- function(c1, c2){sqrt(sum((c1-c2)^2))}
shrink <- function(x, a = 0.1){c(x[1] + diff(x)/2*a, x[2] - diff(x)/2*a)}

#starting color spec
motrpac_cols = T
brewer_cols = F
rainbow_circle_cols = F
n_rainbow = 6
n_optim_runs <- 20

compress_saturation_by = 0
compress_luminance_by = 0
n_cols_to_add = 4


if(brewer_cols){
  cols <- RColorBrewer::brewer.pal(8, "Set1")
  dat <- orig_cols <- farver::decode_colour(cols, to = "hsl")
}

if(motrpac_cols){
  cols <- unique(MotrpacBicQC::bic_animal_tissue_code$tissue_hex_colour)
  cols <- cols[!is.na(cols)]
  dat <- orig_cols <- farver::decode_colour(cols, to = "hsl")
}

if(rainbow_circle_cols){
  dat <- do.call(rbind, lapply(1:n_rainbow, function(x) farver::decode_colour("red", to = "hsl")) )
  dat[,"h"] <- seq(0, 360, length.out = nrow(dat))  
}

dat[,"s"] <- dat[,"s"] / 100
dat[,"h"] <- dat[,"h"] / 360 
dat[,"l"] <- dat[,"l"] / 100

dat_orig <- dat
dat_orig[,"h"] <- dat_orig[,"h"] * 360
dat_orig[,"s"] <- dat_orig[,"s"] * 100
dat_orig[,"l"] <- dat_orig[,"l"] * 100
cols <- encode_colour(dat_orig, from = "hsl")

#set maximum and minimum bounds to saturtation and luminance
bounds_saturation = quantile((dat[,"s"]), probs = c(0.05,0.95))
bounds_luminance = quantile((dat[,"l"]), probs = c(0.1,0.7))

bounds_saturation = c(0.1,0.9)
bounds_luminance = c(0.3,0.6)

dat <- list(old_cols = dat, 
            bounds_saturation = bounds_saturation, bounds_luminance = bounds_luminance, 
            compress_saturation_by = compress_saturation_by, compress_luminance_by = compress_luminance_by)

par <- lapply(1:n_optim_runs, function(x) rnorm(n = n_cols_to_add*3, sd = 2))

mean_dist <- function(par, dat){
  
  #extract and process dat info
  bounds_saturation = dat$bounds_saturation
  bounds_luminance = dat$bounds_luminance
  compress_saturation_by = dat$compress_saturation_by
  compress_luminance_by = dat$compress_luminance_by
  dat = dat$old_cols
  dat[,"h"] <- dat[,"h"] * 2 * pi
  
  #configure params to appropriate values
  #chroma & luminance bounded between 0 and 1, hue between 0 and 2pi
  pars <- pars_logit_scale <- matrix(par, ncol = 3, byrow = T)
  colnames(pars) <- colnames(dat)
  pars[2:nrow(pars),"l"] <- exp(pars[2:nrow(pars),"l"]) #for identifiability
  pars[,"l"] <- cumsum(pars[,"l"]) #for identifiability
  pars[,"l"] <- invlogit(pars[,"l"]) * diff(bounds_luminance) + bounds_luminance[1]
  pars[,"s"] <- invlogit(pars[,"s"]) * diff(bounds_saturation) + bounds_saturation[1] 
  pars[,"h"] <- invlogit(pars[,"h"]) * 2 * pi 
  
  #now compress luminance and saturation
  pars[,"l"] <- invlogit(pars[,"l"]) * (1-compress_luminance_by) + compress_luminance_by/2
  pars[,"s"] <- invlogit(pars[,"s"]) * (1-compress_saturation_by) + compress_saturation_by/2
  dat[,"l"] <- invlogit(dat[,"l"]) * (1-compress_luminance_by) + compress_luminance_by/2
  dat[,"s"] <- invlogit(dat[,"s"]) * (1-compress_saturation_by) + compress_saturation_by/2
  
  #compute distances
  dists_dat <- c(sapply(1:nrow(pars), function(ri) sapply(1:nrow(dat), function(ci) 
    dist_kd(cyl2cart(pars[ri,1], pars[ri,2], pars[ri,3]), cyl2cart(dat[ci,1], dat[ci,2], dat[ci,3])) 
            )))
  
  if(nrow(pars) > 1){
    dists_par <- t(combn(x = 1:nrow(pars), m = 2))
    dists_par <- sapply(1:nrow(dists_par), function(ri) 
      dist_kd(cyl2cart(pars[dists_par[ri,1],1], pars[dists_par[ri,1],2], pars[dists_par[ri,1],3]), 
              cyl2cart(pars[dists_par[ri,2],1], pars[dists_par[ri,2],2], pars[dists_par[ri,2],3])))
    dists <- c(dists_par, dists_dat)
  } else {
    dists <- dists_dat
  }
  
  # print(min(dists))
  
  #compute objective function
  mean_dists <- (1/mean((1/dists)^0.5)) #harmonic mean
  # mean_dists <- mean(dists) #arithmetic mean
  # mean_dists <- exp(mean(log(dists))) #geometric mean
  target <- -mean_dists
  # target <- -min(dists)
  
  #regularize at logit-scale boundaries w/ double exponential
  rate_reg <- 1E-3
  target <- target - sum(dexp(abs(c(pars_logit_scale[,-1])), rate = rate_reg, log = T)) / 10
  
  return(target)
}

mean_dist(par[[1]], dat)


# output <- output_orig <- optimx::optimx(par = par, fn = mean_dist, dat = dat, method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "Rvmmin", "CG"), control=list(kkt=FALSE))
# output <- optim(par = par, fn = mean_dist, dat = dat, method = c("Nelder-Mead"), hessian = T)
output <- output_orig <- do.call(rbind, lapply(1:n_optim_runs, function(x) optimx::optimx(par = par[[x]], fn = mean_dist, dat = dat, method = c("nlm"), control=list(kkt=FALSE))))
# solve(-output$hessian)
output <- output[which.min(output$value),]
raw_pars <- unlist(output[paste0("p", 1:length(par[[1]]))])

mean_dist(par = raw_pars, dat = dat)

par_est <- matrix(raw_pars, ncol = 3, byrow = T)
colnames(par_est) <- colnames(dat$old_cols)
par_est[2:nrow(par_est),"l"] <- exp(par_est[2:nrow(par_est),"l"]) #for identifiability
par_est[,"l"] <- cumsum(par_est[,"l"]) #for identifiability
par_est[,"l"] <- (invlogit(par_est[,"l"]) * diff(bounds_luminance) + bounds_luminance[1]) * 100
par_est[,"h"] <- invlogit(par_est[,"h"]) * 360
par_est[,"s"] <- (invlogit(par_est[,"s"]) * diff(bounds_saturation) + bounds_saturation[1]) * 100
new_cols <- farver::encode_colour(par_est, from = "hsl")


#plot colorspace w/ old and new colors
par(mfrow = c(2,1), mar = c(0,0,0,0))
plot(x = NULL, y = NULL, pch = 15, col = cols, cex = 10, xlim = c(0,5), ylim = c(0, 1),
     frame.plot = F, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
max_col_num <- max(length(cols), length(new_cols))
col_height <- 0.9 / max_col_num
rect(xleft = rep(0.5, length(cols)), xright = rep(1.5, length(cols)), 
     ybottom = cumsum(c(0, rep(col_height, length(cols)-1))),
     ytop = cumsum(rep(col_height, length(cols))), col = cols)
rect(xleft = rep(2, length(new_cols)), xright = rep(3, length(new_cols)), 
     ybottom = cumsum(c(0, rep(col_height, length(new_cols)-1))),
     ytop = cumsum(rep(col_height, length(new_cols))), col = new_cols)

#color hex values
text(x = 1, y = cumsum(c(0, rep(col_height, length(cols)-1))) + col_height / 2, labels = cols, col = "white")
text(x = 2.5, y = cumsum(c(0, rep(col_height, length(new_cols)-1))) + col_height / 2, labels = new_cols, col = "white")

#color labels
text(1, length(cols) * col_height, cex = 1.5, "old colors", srt = 0, pos = 3)
text(2.5, length(new_cols) * col_height, cex = 1.5, "new colors", srt = 0, pos = 3)
points(1.75, length(cols) * col_height + 0.04, cex = 2, pch = 1)
points(3.325, length(new_cols) * col_height + 0.04, cex = 1.75, pch = 0, xpd = NA)

#luminosity
nsegs <- 50
rect(xleft = rep(3.9, nsegs), xright = rep(4.1, nsegs), ybottom = seq(0,0.5, length.out = nsegs), ytop = seq(0,0.5, length.out = nsegs) + 0.5 / nsegs, 
         border = farver::encode_colour(matrix(nrow = nsegs, ncol = 3, data = c(rep(0, nsegs), rep(0,nsegs), seq(0,100, length.out = nsegs))), from = "hsl"),
         col = farver::encode_colour(matrix(nrow = nsegs, ncol = 3, data = c(rep(0, nsegs), rep(0,nsegs), seq(0,100, length.out = nsegs))), from = "hsl"))
rect(xleft = 3.9, xright = 4.1, ybottom = 0, ytop = 0.5 + 1 / nsegs)
points(x = rep(3.8, 10) - seq(0,0.1, length.out = 10), y = seq(0,0.5, length.out = 10), pch = 21, cex = seq(0,100, length.out = 10) / 50 + 0.5)
points(x = rep(4.2, 10) + seq(0,0.1, length.out = 10), y = seq(0,0.5, length.out = 10), pch = 22, cex = seq(0,100, length.out = 10) / 50 + 0.5)
text(4, 0.5 + 1 / nsegs, cex = 1.5, "luminosity", srt = 0, pos = 3)


#plot color wheel
plot(1,1,col = "white", xlim = c(-1,1), ylim = c(-1,1), frame.plot = F, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
n_slices_t <- 60
n_slices_r <- 20
degrees <- seq(0,2*pi, length.out = n_slices_t+1)
dd <- diff(degrees)[1]/2
radii <- seq(0,1, length.out = n_slices_r+1)
rd <- diff(radii)[1]/2
for(ts in 1:n_slices_t){
  for(rs in 1:n_slices_r){
    coords <- rbind(polar2cart(t = degrees[ts], r = radii[rs]),
                polar2cart(t = degrees[ts], r = radii[rs+1]),
                polar2cart(t = degrees[ts+1], r = radii[rs+1]),
                polar2cart(t = degrees[ts+1], r = radii[rs]))
    coords <- rbind(coords, coords[1,])
    polygon(x = coords[,1], y = coords[,2], border = NA, col = encode_colour(t(c((degrees[ts] + dd) / 2 / pi * 360, (radii[rs] + rd) * 100, 50)), from = "hsl"))
  }
}

coords_old_cols <- t(sapply(1:nrow(dat$old_cols), function(ri) polar2cart(t = dat$old_cols[ri,"h"] * 2 * pi, r = dat$old_cols[ri,"s"])))
points(x = coords_old_cols[,1], y = coords_old_cols[,2], pch = 21, bg = cols, col = "black", cex = dat$old_cols[,"l"] * 2 + 0.5)

coords_new_cols <- t(sapply(1:nrow(par_est), function(ri) polar2cart(t = par_est[ri,"h"] / 360 * 2 * pi, r = par_est[ri,"s"] / 100)))
points(x = coords_new_cols[,1], y = coords_new_cols[,2], pch = 22, bg = new_cols, col = "black", cex = par_est[,"l"] / 50 + 0.5)

output$value
