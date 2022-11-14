#functions
addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}

## for debugging
# group_diff = 1
# outlier_threshold = 0.1

#this function computes tail probabilities for a mixture of two normals, both numerically and analytically
#the latter is done in terms of the 'inner', 'middle', and 'outer' integrals of the consitutent normals, to ease of further investigation
prop_outliers <- function(group_diff, outlier_threshold, numerical_integration = F, range = seq(-10, 10, by = 1E-2), sd1 = 1, sd2 = 1, w1 = 0.5, returnAll = F){
  
  w2 <- 1 - 0.5 #weight of second normal, deterministically computed
  # mu1 <- 0 - group_diff / 2 #the first normal's mean, to the left
  # mu2 <- 0 + group_diff / 2 #the second normal's mean, to the right
  # alternatively, for unequal weights to ensure E(mixture) = 0 (for convenience, as an arbitrary translation step)
  mu1 <- -group_diff / (1 + w1 / w2) #the first normal's mean, to the left
  mu2 <- mu1 + group_diff #the second normal's mean, to the right
  
  
  #from https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians/16609#16609
  sd_comb <- sqrt(w1*sd1^2 + w2*sd2^2 + w1*mu1^2 + w2*mu2^2 - (w1*mu1 + w2*mu2)^2) #the standard deviation of the mixture distribution
  
  if(numerical_integration){
    
    dens1 <- dnorm(range, mu1, sd1)
    dens2 <- dnorm(range, mu2, sd2)
    dens_comb <- dens1 * w1 + dens2 * w2
    # sum(dens_comb * binsize) #sanity check
    cumulDens <- cumsum(dens_comb) * diff(range)[1] #assumes equally spaced bins
    prop_outliers_mix <- cumulDens[max(which(range < (-outlier_threshold * sd_comb)))] + (1 - cumulDens[min(which(range > (outlier_threshold * sd_comb)))])
    #equivalently, due to symmetry: 2 * prop_outliers <- cumulDens[max(which(range < (-outlier_threshold * sd_comb)))], but we may want to vary sds etc. later
    
  } else { #alternatively, let's work with the cdfs
    
    #take individual probabilities for funsies
    inner_prop1 <- (1-pnorm(mu1 + outlier_threshold*sd1, mu1, sd1))*w1
    outer_prop1 <- pnorm(mu1 - outlier_threshold*sd1, mu1, sd1)*w1
    inner_prop2 <- pnorm(mu2 - outlier_threshold*sd2, mu2, sd2)*w2
    outer_prop2 <- (1-pnorm(mu2 + outlier_threshold*sd2, mu2, sd2))*w2
    
    #integrate parts of normal #2
    inner_prop_mix2 <- pnorm(min(-outlier_threshold * sd_comb, mu2 - outlier_threshold * sd2), mu2, sd2) * w2 + 
      max(pnorm(mu2 - outlier_threshold * sd2, mu2, sd2) - pnorm(outlier_threshold * sd_comb, mu2, sd2), 0) * w2
    outer_prop_mix2 <- min((1-pnorm(mu2 + outlier_threshold * sd2, mu2, sd2))*w2, (1-pnorm(outlier_threshold * sd_comb, mu2, sd2))*w2)
    middle_prop_mix2 <- max(
      (pnorm(mu2 + outlier_threshold * sd2, mu2, sd2) - pnorm(max(outlier_threshold * sd_comb, mu2 - outlier_threshold * sd2), mu2, sd2))*w2, 
      0)
    
    #integrate parts of normal #1
    inner_prop_mix1 <- (1 - pnorm(max(mu1 + outlier_threshold * sd1, sd_comb * outlier_threshold), mu1, sd1)) * w1 + 
      max(pnorm(-outlier_threshold * sd_comb, mu1, sd1) - pnorm(mu1 + outlier_threshold*sd1, mu1, sd1), 0) * w1
    outer_prop_mix1 <- min((pnorm(mu1 - outlier_threshold * sd1, mu1, sd1))*w1, (pnorm(-outlier_threshold * sd_comb, mu1, sd1))*w1)
    middle_prop_mix1 <- max(
      (pnorm(min(-outlier_threshold * sd_comb, mu1 + outlier_threshold * sd1), mu1, sd1) - pnorm(mu1 - outlier_threshold * sd1, mu1, sd1))*w1, 
      0)
    
    #similar symmetry here as before w/ equal variances etc.
    
    #sum the probabilities (or return them as components)
    prop_outliers_subpop <- inner_prop1 + outer_prop1 + inner_prop2 + outer_prop2
    prop_outliers_mix <- inner_prop_mix1 + outer_prop_mix1 + middle_prop_mix1 + inner_prop_mix2 + outer_prop_mix2 + middle_prop_mix2
    
  }
  if(!returnAll){
    return(prop_outliers_mix)
  } else {
    #these are all the contributions of the subpops to the tails of the mixture pop
    return(list(prop_outliers_mix = prop_outliers_mix, 
                inner_prop_mix1 = inner_prop_mix1, 
                outer_prop_mix1 = outer_prop_mix1, 
                middle_prop_mix1 = middle_prop_mix1, 
                inner_prop_mix2 = inner_prop_mix2, 
                outer_prop_mix2 = outer_prop_mix2, 
                middle_prop_mix2 = middle_prop_mix2,
                inner_prop1 = inner_prop1,
                outer_prop1 = outer_prop1,
                inner_prop2 = inner_prop2,
                outer_prop2 = outer_prop2,
                prop_outliers_subpop = prop_outliers_subpop))
  }
}

## test to confirm numerical vs analytical match
prop_outliers(group_diff = group_diff, outlier_threshold = outlier_threshold, numerical_integration = T, range = seq(-10, 10, by = 1E-5))
prop_outliers(group_diff = group_diff, outlier_threshold = outlier_threshold, numerical_integration = F, range = seq(-10, 10, by = 1E-5))
# test <- lapply(runif(10, 0, 5), function(ot) t(sapply(runif(10, 0, 5), function(gd) c(
#   num = prop_outliers(group_diff = gd, outlier_threshold = ot, numerical_integration = T, range = seq(-10, 10, by = 1E-4)),
#   ana = prop_outliers(group_diff = gd, outlier_threshold = ot, numerical_integration = F, range = seq(-10, 10, by = 1E-4))
# ))))
# test <- do.call(rbind, test)
# test <- test[!is.na(apply(test, 1, sum)),]
# plot(test); abline(0,1)
## dat speedup
# microbenchmark::microbenchmark(prop_outliers(group_diff = 2, outlier_threshold = 1, numerical_integration = T, range = seq(-10, 10, by = 1E-3)),
#                                prop_outliers(group_diff = 2, outlier_threshold = 1, numerical_integration = F, range = seq(-10, 10, by = 1E-3))
#                                )

#### make a static figure ####

#explore range of conditions
# binsize <- 5E-3
# range <- seq(-20, 20, by = binsize)
group_diffs <- seq(0, 8, length.out = 40)
outlier_thresholds <- 0:1600/200 #in z-scores
everything_together <- F
inners_only <- F
outers_only <- T
redo_number_crunching <- T
if(redo_number_crunching | !exists("d")){
  if(everything_together){
    d <- lapply(group_diffs, function(gd) t(sapply(outlier_thresholds, function(ot) c(
        indiv = (1 - pnorm(ot, mean = 0, sd = 1))*2,
        comb = prop_outliers(group_diff = gd, outlier_threshold = ot, range = range, numerical_integration = F)))))
    extra = ""
  }
  if(inners_only){
    d <- lapply(group_diffs, function(gd) t(sapply(outlier_thresholds, function(ot) c(
      indiv = (1 - pnorm(ot, mean = 0, sd = 1))*2,
      comb = prop_outliers(group_diff = gd, outlier_threshold = ot, range = range, numerical_integration = F, returnAll = T)$inner_prop_mix1 * 4))))
    extra <- "_innersOnly"
  }
  if(outers_only){
    d <- lapply(group_diffs, function(gd) t(sapply(outlier_thresholds, function(ot) c(
      indiv = (1 - pnorm(ot, mean = 0, sd = 1))*2,
      comb = prop_outliers(group_diff = gd, outlier_threshold = ot, range = range, numerical_integration = F, returnAll = T)$outer_prop_mix1 * 4))))
    extra <- "_outersOnly"
  }
}

#do the actual plotting
cols <- viridis::plasma(n = length(d), end = 0.95)
grDevices::cairo_pdf(filename = paste0("~/Documents/prop_outliers", extra, ".pdf"), width = 800 / 72, height = 800 / 72, family="Arial Unicode MS")
par(mar = c(3.5,4,1,1), xpd = F)

plot(1,1,xlim = c(0,1), ylim = c(0,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
mtext(text = "Subpopulation Proportion Outside Given Z-Score", cex = 1.75, line = 2, side = 1) #horiz axis label
mtext(text = "Mixture Proportion Outside Given Z-Score", cex = 1.75, line = 2, side = 2) #vert axis label
axis(1, at = 0:10/10, labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.9)
mtext(text = 0:10/10, side = 1, at = 0:10/10, cex = 1.25, line = -0.5)
axis(2, at = 0:10/10, labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.9)
mtext(text = 0:10/10, side = 2, at = 0:10/10, cex = 1.25, line = -0.5)

#plot color legend
xl <- 0.01; yb <- 0.5; xr <- 0.05; yt <- 0.955;
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/length(cols)),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/length(cols)),-1),
  col=cols
)
text(labels = round(seq(from = range(group_diffs)[1], to = range(group_diffs)[2], length.out = 9)), y = seq(yb, yt - 0.01, length.out = 9), x = xr, pos = 4, las=2, cex=1.5)
text(labels = latex2exp::TeX("$\\Delta_{\\mu$"), pos = 4, x = xl - 0.0075, y = yt + 0.015, cex = 2, font = 2)

#subplot corners & connecting lines
if(inners_only){
  xl <- 0.125; yb <- 0.55; xr <- 0.53; yt <- 0.955;
} else if(outers_only){
  xl <- 0.55; yb <- 0.025; xr <- 0.975; yt <- 0.45;
} else {
  xl <- 0.55; yb <- 0.025; xr <- 0.975; yt <- 0.45; 
}
x0 <- 0; y0 <- 0; x1 <- 0.1; y1 <- 0.1;
# segments(x0 = x0, y0 = y1, x1 = xl, y1 = yt, lty = 2, lwd = 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y0, x1 = xr, y1 = yb, lty = 2, lwd = 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y1, x1 = xr, y1 = yt, lty = 2, lwd = 1, col = rgb(0,0,0,0.4))
# segments(x0 = x0, y0 = y0, x1 = xl, y1 = yb, lty = 2, lwd = 1, col = rgb(0,0,0,0.4))

for(i in 1:length(d)){
  lines(d[[i]], col = cols[i], lwd = 1.1)
}

#make subplot on log-log scale
rect(x0, y0, x1, y1, lwd = 1, col = rgb(0,0,0,0.05))
bounds_of_subplot <- c(-1, -10)
for(i in 1:length(d)){
  pts2plot <- which(log10(d[[i]][,1]) < bounds_of_subplot[1] & log10(d[[i]][,1]) > bounds_of_subplot[2] & 
                      log10(d[[i]][,2]) < bounds_of_subplot[1] & log10(d[[i]][,2]) > bounds_of_subplot[2])
  xs <- log10(d[[i]][pts2plot,1])
  ys <- log10(d[[i]][pts2plot,2])
  #translate to origin
  xs <- xs - bounds_of_subplot[2]
  ys <- ys - bounds_of_subplot[2]
  #rescale to subplot
  xs <- xs / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (xr-xl)
  ys <- ys / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (yt-yb)
  #translate to subplot region
  xs <- xs + xl
  ys <- ys + yb
  #plot
  lines(xs, y = ys, col = cols[i], lwd = 2)
}

rect(xl, yb, xr, yt, lwd = 2, col = rgb(0,0,0,0.05))
segments(x0 = xl, y0 = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1), x1 = xl - 0.0075, lwd = 2,
         y1 = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1)) #horiz
segments(x0 = seq(xl, xr, length.out = abs(diff(bounds_of_subplot)) + 1), y0 = yt, 
         x1 = seq(xl, xr, length.out = abs(diff(bounds_of_subplot)) + 1), lwd = 2,
         y1 = yt + 0.00735) #vert
text(x = xl - 0.0075, y = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1) + 0.02, pos = 2, srt = 90,
     labels = c(sapply(bounds_of_subplot[2]:(bounds_of_subplot[1]), function(i) as.expression(bquote(10^ .(i)))))) #vert 
text(x = seq(xl, xr, length.out = abs(diff(bounds_of_subplot))+1), y = yt, pos = 3,
     labels = (c(sapply(bounds_of_subplot[2]:(bounds_of_subplot[1]), function(i) as.expression(bquote(10^ .(i))))))) #horiz

# segments(0,0,1,1, lwd = 2, col = rgb(0,0,0,0.15))
text(x = 0.9, y = 0.9, labels = "1-to-1 line", pos = ifelse(inners_only, 3, 1), srt = 45, cex = 1.3)
rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, lwd = 3)

#title 
if(inners_only){
  title(main = "inner probabilities only", cex.main = 2, line = -1)
} else if(outers_only){
  title(main = "outer probabilities only", cex.main = 2, line = -1)
}

dev.off()


#### let's draw an explanatory figure too ####
grDevices::cairo_pdf(filename = paste0("~/Documents/zscore_terminology.pdf"), width = 800 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mar = c(4,0,0,2), xpd = F)
cols2 <- c(cols[1], cols[21], cols[11])
xlims <- c(-4,4); ylims <- c(0,2)
plot(1,1,xlim = xlims, ylim = ylims, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sd1 = 0.5
sd2 = 0.65
z = 2
mu1 = -0.75
mu2 = 0.75
range <- seq(mu1 - 4*sd1, mu2 + 4*sd2, by = 0.01)
w1 <- w2 <- 0.5
mu3 <- mu1 * w1 + mu2 * w2
sd3 <- sqrt(w1*sd1^2 + w2*sd2^2 + w1*mu1^2 + w2*mu2^2 - (w1*mu1 + w2*mu2)^2) #the standard deviation of the mixture distribution
dist1 <- dnorm(range, mu1, sd1)
dist2 <- dnorm(range, mu2, sd2)
dist3 <- (dist1 + dist2)
if(max(dist3) > 0.9){dist3 <- dist3 / max(dist3) * 0.9; text("objects in figure may be \nlarger than they appear",x = -3, y = 0, pos = 3)}
lines(range, dist1 + 1, lwd = 2, col = cols2[1])
lines(range, dist2 + 1, lwd = 2, col = cols2[2])
lines(range, dist3, lwd = 2, col = cols2[3])
polygon(x = range, y = dist1 + 1, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2))
polygon(x = range, y = dist2 + 1, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2))
polygon(x = range, y = dist3, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2))
text(labels = latex2exp::TeX("$\\mu_{1}$"), pos = 3, x = range[which.max(dist1)], y = max(dist1) + 1, cex = 2, font = 2, col = cols2[1])
text(labels = latex2exp::TeX("$\\mu_{2}$"), pos = 3, x = range[which.max(dist2)], y = max(dist2) + 1, cex = 2, font = 2, col = cols2[2])
text(labels = latex2exp::TeX("$\\mu_{mix}$"), pos = 3, x = mu1 / 2 + mu2 / 2, y = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, cex = 2, font = 2, col = cols2[3])
#bounds for dist1
segments(x0 = mu1 + sd1 * z, x1 = mu1 + sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 + sd1 * z), decreasing = F)[1]], col = cols2[1], lwd = 2, lty = 2)
segments(x0 = mu1 - sd1 * z, x1 = mu1 - sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 - sd1 * z), decreasing = F)[1]], col = cols2[1], lwd = 2, lty = 2)
polygon(x = c(range[range < mu1 - sd1 * z], mu1 - sd1 * z), y = c(dist1[range < mu1 - sd1 * z] + 1, 1), col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu1 + sd1 * z], mu1 + sd1 * z), y = c(dist1[range > mu1 + sd1 * z] + 1, 1), col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2), border = NA)
#bounds for dist2
segments(x0 = mu2 + sd2 * z, x1 = mu2 + sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 + sd2 * z), decreasing = F)[1]], col = cols2[2], lwd = 2, lty = 2)
segments(x0 = mu2 - sd2 * z, x1 = mu2 - sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 - sd2 * z), decreasing = F)[1]], col = cols2[2], lwd = 2, lty = 2)
polygon(x = c(range[range < mu2 - sd2 * z], mu2 - sd2 * z), y = c(dist2[range < mu2 - sd2 * z] + 1, 1), col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu2 + sd2 * z], mu2 + sd2 * z), y = c(dist2[range > mu2 + sd2 * z] + 1, 1), col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2), border = NA)
#bounds for dist3
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = cols2[3], lwd = 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = cols2[3], lwd = 2, lty = 2)
polygon(x = c(range[range < mu3 - sd3 * z], mu3 - sd3 * z), y = c(dist3[range < mu3 - sd3 * z], 0), col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu3 + sd3 * z], mu3 + sd3 * z), y = c(dist3[range > mu3 + sd3 * z], 0), col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2), border = NA)
#mark means with lines?
segments(x0 = mu1, x1 = mu1, y0 = dist2[order(abs(range - mu1), decreasing = F)[1]] + 1, y1 = max(dist1) + 1, col = cols2[1], lwd = 1, lty = 3)
segments(x0 = mu2, x1 = mu2, y0 = dist1[order(abs(range - mu2), decreasing = F)[1]] + 1, y1 = max(dist2) + 1, col = cols2[2], lwd = 1, lty = 3)
segments(x0 = mu3, x1 = mu3, y0 = 0, y1 = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, col = cols2[3], lwd = 1, lty = 3)

#labels for bounds
oid <- max(min((abs((mu1 + sd1*z) - (mu2 - sd2*z)) - 0.5) / 1.5, 0), -0.055) #overlapping_inner_disp
# strwidth(latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), font = 2, units = "figure") * par("pin")[1] / par("fin")[1] * diff(xlims) / 2
text(labels = latex2exp::TeX(paste0("$\\mu_{1} - $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 - sd1*z, y = 1, cex = 1, font = 2, col = cols2[1])
text(labels = latex2exp::TeX(paste0("$\\mu_{1} + $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 + sd1*z, y = 1 + oid, cex = 1, font = 2, col = cols2[1])
text(labels = latex2exp::TeX(paste0("$\\mu_{2} - $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 - sd2*z, y = 1, cex = 1, font = 2, col = cols2[2])
text(labels = latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 + sd2*z, y = 1, cex = 1, font = 2, col = cols2[2])
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} - $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 - sd3*z, y = 0, cex = 1, font = 2, col = cols2[3])
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} + $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 + sd3*z, y = 0, cex = 1, font = 2, col = cols2[3])
#add some bars to separate the mus
# segments(x0 = mu1, x1 = mu2, y0 = max(dist1) + 1.15, y1 = max(dist2) + 1.15, lwd = 1.5)
# segments(x0 = mu1, x1 = mu1, y0 = max(dist1) + 1.125, y1 = max(dist1) + 1.175, lwd = 1.5)
# segments(x0 = mu2, x1 = mu2, y0 = max(dist2) + 1.125, y1 = max(dist2) + 1.175, lwd = 1.5)
# text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = (max(dist2) + max(dist1)) / 2 + 1.1275, 
#      srt = atan(((max(dist2) - max(dist1)) / (diff(ylims)) * par("pin")[2]) / ((mu2 - mu1) / diff(xlims) * par("pin")[1])) * 180 / pi)
#make it a horizontal bar only
segments(x0 = mu1, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.15, y1 = max(c(dist1, dist2)) + 1.15, lwd = 1.5)
segments(x0 = mu1, x1 = mu1, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd = 1.5)
segments(x0 = mu2, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd = 1.5)
text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = max(c(dist1, dist2)) + 1.1275, 
     srt = 0)
#show inclusion / exclusion zone
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = "grey50", lwd = 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = "grey50", lwd = 2, lty = 2)
text(x = mu3 + sd3 * z, y = 0.65, labels = "excluded", srt = 90, pos = 2, col = "grey50")
text(x = mu3 + sd3 * z, y = 0.65, labels = "included", srt = 270, pos = 4, col = "grey50")
text(x = mu3 - sd3 * z, y = 0.65, labels = "included", srt = 90, pos = 2, col = "grey50")
text(x = mu3- sd3 * z, y = 0.65, labels = "excluded", srt = 270, pos = 4, col = "grey50")

#label regions of subpopulation kdes
text(x = mu1 + sd1 * z, y = 1.05, labels = "m", srt = 90, pos = 2, col = cols2[1], cex = 0.75)
text(x = mu1 + sd1 * z, y = 1.05, labels = "in", srt = 270, pos = 4, col = cols2[1], cex = 0.75)
text(x = mu1 - sd1 * z, y = 1.05, labels = "ot", srt = 90, pos = 2, col = cols2[1], cex = 0.75)
text(x = mu1- sd1 * z, y = 1.05, labels = "m", srt = 270, pos = 4, col = cols2[1], cex = 0.75)
text(x = mu2 + sd2 * z, y = 1.05, labels = "m", srt = 90, pos = 2, col = cols2[2], cex = 0.75)
text(x = mu2 + sd2 * z, y = 1.05, labels = "ot", srt = 270, pos = 4, col = cols2[2], cex = 0.75)
text(x = mu2 - sd2 * z, y = 1.05, labels = "in", srt = 90, pos = 2, col = cols2[2], cex = 0.75)
text(x = mu2- sd2 * z, y = 1.05, labels = "m", srt = 270, pos = 4, col = cols2[2], cex = 0.75)

#draw little bargraph for tail proportions
po <- prop_outliers(group_diff = mu2-mu1, outlier_threshold = z, sd1 = sd1, sd2 = sd2, returnAll = T)
segments(x0 = 3.5, x1 = 4, y0 = 0, y1 = 0, lwd= 2)
rect(xleft = 3.5, xright = 3.75, ybottom = 0, ytop = (po$inner_prop1 + po$outer_prop1) / 0.25, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.36))
rect(xleft = 3.5, xright = 3.75, ybottom = (po$inner_prop1 + po$outer_prop1) / 0.25, ytop = (po$inner_prop1 + po$outer_prop1) / 0.25 + (po$inner_prop2 + po$outer_prop2) / 0.25, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.36))
rect(xleft = 3.75, xright = 4, ybottom = 0, ytop = po$prop_outliers_mix / 0.25, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.36))
# text(x = 3.75, pos = 4, y = 1, labels = expression(integral()))
# addImg(png::readPNG("~/Documents/integral_symbol.png"), x = 3.65, y = -0.05, width = 0.075)
text(x = 3.625, pos = 1, y = 0.0125, labels = "sub")
text(x = 3.875, pos = 1, y = 0.0125, labels = "mix")
dev.off()

#### static figure redux per rachel's suggestion ####
out <- prop_outliers(group_diff = 1, outlier_threshold = 2.5, returnAll = T)

#gained outliers... not accommodating dist2's contribution to dist1's tail :/
out$inner_prop1 - out$inner_prop_mix1
out$inner_prop2 - out$inner_prop_mix2

#lost outliers... not accommodating 
out$outer_prop1 - out$outer_prop_mix1



#### animation ####
# ok, and now instead of just an explanatory figure, let's make an explanatory animation -- first, let's try making the base figure
###

png(filename = paste0("~/Documents/zscore_outlier_animation/zscore_cartoon.png"), width = 2600, height = 1000, family="Arial Unicode MS")
layout((t(c(1,1,1,2,2))))
par(mar = c(4,0,0,2), xpd = F)
cols2 <- c(cols[1], cols[21], cols[11])
xlims <- c(-4,4); ylims <- c(0,2)
plot(1,1,xlim = xlims, ylim = ylims, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sd1 = 0.5
sd2 = 0.65
z = 2
mu1 = -0.75
mu2 = 0.75
range <- seq(mu1 - 4*sd1, mu2 + 4*sd2, by = 0.01)
segments(x0 = min(range), x1 = max(range), y0 = c(0,1), lwd = 2)
w1 <- w2 <- 0.5
mu3 <- mu1 * w1 + mu2 * w2
sd3 <- sqrt(w1*sd1^2 + w2*sd2^2 + w1*mu1^2 + w2*mu2^2 - (w1*mu1 + w2*mu2)^2) #the standard deviation of the mixture distribution
dist1 <- dnorm(range, mu1, sd1)
dist2 <- dnorm(range, mu2, sd2)
dist3 <- (dist1 + dist2)
if(max(dist3) > 0.9){dist3 <- dist3 / max(dist3) * 0.9; text("objects in figure may be \nlarger than they appear",x = -3, y = 0, pos = 3)}
lines(range, dist1 + 1, lwd =1.75* 2, col = cols2[1])
lines(range, dist2 + 1, lwd =1.75* 2, col = cols2[2])
lines(range, dist3, lwd =1.75* 2, col = cols2[3])
polygon(x = range, y = dist1 + 1, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2))
polygon(x = range, y = dist2 + 1, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2))
polygon(x = range, y = dist3, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2))
text(labels = latex2exp::TeX("$\\mu_{1}$"), pos = 3, x = range[which.max(dist1)], y = max(dist1) + 1, cex =4, font = 2, col = cols2[1])
text(labels = latex2exp::TeX("$\\mu_{2}$"), pos = 3, x = range[which.max(dist2)], y = max(dist2) + 1, cex =4, font = 2, col = cols2[2])
text(labels = latex2exp::TeX("$\\mu_{mix}$"), pos = 3, x = mu1 / 2 + mu2 / 2, y = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, cex =4, font = 2, col = cols2[3])
#bounds for dist1
segments(x0 = mu1 + sd1 * z, x1 = mu1 + sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 + sd1 * z), decreasing = F)[1]], col = cols2[1], lwd =1.75* 2, lty = 2)
segments(x0 = mu1 - sd1 * z, x1 = mu1 - sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 - sd1 * z), decreasing = F)[1]], col = cols2[1], lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu1 - sd1 * z], mu1 - sd1 * z), y = c(dist1[range < mu1 - sd1 * z] + 1, 1), col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu1 + sd1 * z], mu1 + sd1 * z), y = c(dist1[range > mu1 + sd1 * z] + 1, 1), col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2), border = NA)
#bounds for dist2
segments(x0 = mu2 + sd2 * z, x1 = mu2 + sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 + sd2 * z), decreasing = F)[1]], col = cols2[2], lwd =1.75* 2, lty = 2)
segments(x0 = mu2 - sd2 * z, x1 = mu2 - sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 - sd2 * z), decreasing = F)[1]], col = cols2[2], lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu2 - sd2 * z], mu2 - sd2 * z), y = c(dist2[range < mu2 - sd2 * z] + 1, 1), col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu2 + sd2 * z], mu2 + sd2 * z), y = c(dist2[range > mu2 + sd2 * z] + 1, 1), col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2), border = NA)
#bounds for dist3
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = cols2[3], lwd =1.75* 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = cols2[3], lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu3 - sd3 * z], mu3 - sd3 * z), y = c(dist3[range < mu3 - sd3 * z], 0), col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu3 + sd3 * z], mu3 + sd3 * z), y = c(dist3[range > mu3 + sd3 * z], 0), col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2), border = NA)
#mark means with lines?
segments(x0 = mu1, x1 = mu1, y0 = dist2[order(abs(range - mu1), decreasing = F)[1]] + 1, y1 = max(dist1) + 1, col = cols2[1], lwd =1.75* 1, lty = 3)
segments(x0 = mu2, x1 = mu2, y0 = dist1[order(abs(range - mu2), decreasing = F)[1]] + 1, y1 = max(dist2) + 1, col = cols2[2], lwd =1.75* 1, lty = 3)
segments(x0 = mu3, x1 = mu3, y0 = 0, y1 = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, col = cols2[3], lwd =1.75* 1, lty = 3)

#labels for bounds
oid <- max(min((abs((mu1 + sd1*z) - (mu2 - sd2*z)) - 0.5) / 1.5, 0), -0.055) #overlapping_inner_disp
# strwidth(latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), font = 2, units = "figure") * par("pin")[1] / par("fin")[1] * diff(xlims) / 2
text(labels = latex2exp::TeX(paste0("$\\mu_{1} - $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 - sd1*z, y = 0.9875, cex = 2.75, font = 2, col = cols2[1])
text(labels = latex2exp::TeX(paste0("$\\mu_{1} + $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 + sd1*z, y = 0.9875 + oid, cex = 2.75, font = 2, col = cols2[1])
text(labels = latex2exp::TeX(paste0("$\\mu_{2} - $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 - sd2*z, y = 0.9875, cex = 2.75, font = 2, col = cols2[2])
text(labels = latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 + sd2*z, y = 0.9875, cex = 2.75, font = 2, col = cols2[2])
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} - $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 - sd3*z, y = -0.0125, cex = 2.75, font = 2, col = cols2[3])
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} + $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 + sd3*z, y = -0.0125, cex = 2.75, font = 2, col = cols2[3])
#add some bars to separate the mus
# segments(x0 = mu1, x1 = mu2, y0 = max(dist1) + 1.15, y1 = max(dist2) + 1.15, lwd =1.75* 1.5)
# segments(x0 = mu1, x1 = mu1, y0 = max(dist1) + 1.125, y1 = max(dist1) + 1.175, lwd =1.75* 1.5)
# segments(x0 = mu2, x1 = mu2, y0 = max(dist2) + 1.125, y1 = max(dist2) + 1.175, lwd =1.75* 1.5)
# text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = (max(dist2) + max(dist1)) / 2 + 1.1275, 
#      srt = atan(((max(dist2) - max(dist1)) / (diff(ylims)) * par("pin")[2]) / ((mu2 - mu1) / diff(xlims) * par("pin")[1])) * 180 / pi)
#make it a horizontal bar only
segments(x0 = mu1, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.15, y1 = max(c(dist1, dist2)) + 1.15, lwd =1.75* 1.5)
segments(x0 = mu1, x1 = mu1, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd =1.75* 1.5)
segments(x0 = mu2, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd =1.75* 1.5)
text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = max(c(dist1, dist2)) + 1.145, 
     srt = 0, cex = 3)
#show inclusion / exclusion zone
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = "grey50", lwd =1.75* 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = "grey50", lwd =1.75* 2, lty = 2)
text(x = mu3 + sd3 * z - 0.025, y = 0.65, labels = "excluded", srt = 90, pos = 2, col = "grey50", cex = 2.5)
text(x = mu3 + sd3 * z + 0.025, y = 0.65, labels = "included", srt = 270, pos = 4, col = "grey50", cex = 2.5)
text(x = mu3 - sd3 * z - 0.025, y = 0.65, labels = "included", srt = 90, pos = 2, col = "grey50", cex = 2.5)
text(x = mu3- sd3 * z + 0.025, y = 0.65, labels = "excluded", srt = 270, pos = 4, col = "grey50", cex = 2.5)

#label regions of subpopulation kdes
text(x = mu1 + sd1 * z  - 0.025, y = 1.05, labels = "m", srt = 90, pos = 2, col = cols2[1], cex =2)
text(x = mu1 + sd1 * z  + 0.025, y = 1.05, labels = "in", srt = 270, pos = 4, col = cols2[1], cex =2)
text(x = mu1 - sd1 * z  - 0.025, y = 1.05, labels = "ot", srt = 90, pos = 2, col = cols2[1], cex =2)
text(x = mu1- sd1 * z  + 0.025, y = 1.05, labels = "m", srt = 270, pos = 4, col = cols2[1], cex =2)
text(x = mu2 + sd2 * z  - 0.025, y = 1.05, labels = "m", srt = 90, pos = 2, col = cols2[2], cex =2)
text(x = mu2 + sd2 * z  + 0.025, y = 1.05, labels = "ot", srt = 270, pos = 4, col = cols2[2], cex =2)
text(x = mu2 - sd2 * z  - 0.025, y = 1.05, labels = "in", srt = 90, pos = 2, col = cols2[2], cex =2)
text(x = mu2- sd2 * z  + 0.025, y = 1.05, labels = "m", srt = 270, pos = 4, col = cols2[2], cex =2)

#draw little bargraph for tail proportions
po <- prop_outliers(group_diff = mu2-mu1, outlier_threshold = z, sd1 = sd1, sd2 = sd2, returnAll = T)
segments(x0 = 3.5, x1 = 4, y0 = 0, y1 = 0, lwd=1.75* 2)
rect(xleft = 3.5, xright = 3.75, ybottom = 0, ytop = (po$inner_prop1 + po$outer_prop1) / 0.25, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.36), lwd = 2)
rect(lwd = 2, xleft = 3.5, xright = 3.75, ybottom = (po$inner_prop1 + po$outer_prop1) / 0.25, ytop = (po$inner_prop1 + po$outer_prop1) / 0.25 + (po$inner_prop2 + po$outer_prop2) / 0.25, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.36))
rect(xleft = 3.75, xright = 4, ybottom = 0, ytop = po$prop_outliers_mix / 0.25, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.36), lwd = 2)
# text(x = 3.75, pos = 4, y = 1, labels = expression(integral()))
# addImg(png::readPNG("~/Documents/integral_symbol.png"), x = 3.65, y = -0.05, width = 0.075)
text(x = 3.625, pos = 1, y = -0.00625, labels = "sub", cex = 2.5)
text(x = 3.875, pos = 1, y = -0.00625, labels = "mix", cex = 2.5)

#now for the second part, my fancy feather figure
par(mar = c(6,3,0,1), xpd = F)

plot(1,1,xlim = c(0,1), ylim = c(0,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
mtext(text = "Subpopulation Proportion Outside Given Z-Score", cex =1.75* 1.75, line = 4, side = 1) #horiz axis label
mtext(text = "Mixture Proportion Outside Given Z-Score", cex =1.75* 1.75, line = 2.5, side = 2) #vert axis label
axis(1, at = 0:10/10, labels = rep("", 11), lwd =1.75* 2, cex.axis = 2, tck = -0.015, line = -3.7)
mtext(text = 0:10/10, side = 1, at = 0:10/10, cex =1.75* 1.25, line = -0.25)
axis(2, at = 0:10/10, labels = rep("", 11), lwd =1.75* 2, cex.axis = 2, tck = -0.015, line = -3.9)
mtext(text = 0:10/10, side = 2, at = 0:10/10, cex =1.75* 1.25, line = -1.5)

#plot color legend
xl <- 0.01; yb <- 0.5; xr <- 0.05; yt <- 0.955;
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/length(cols)),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/length(cols)),-1),
  col=cols
)
text(labels = round(seq(from = range(group_diffs)[1], to = range(group_diffs)[2], length.out = 9)), y = seq(yb, yt - 0.01, length.out = 9), x = xr, pos = 4, las=2, cex=1.75*1.5)
text(labels = latex2exp::TeX("$\\Delta_{\\mu$"), pos = 4, x = xl - 0.0015, y = yt + 0.015, cex =4, font = 2)

#annotate the fancy graph w/ arrow for delta mu
arrows(x0 = 0.125, x1 = 0.0875, y0 = (mu2-mu1) * (yt-yb) / 8 + yb, y1 = (mu2-mu1) * (yt-yb) / 8 + yb, lwd = 5, col = "darkred")
text(x = 0.125, pos = 4, y = (mu2-mu1) * (yt-yb) / 8 + yb, cex = 3, font = 3, labels = "currently", col = "darkred")

#subplot corners & connecting lines
if(inners_only){
  xl <- 0.125; yb <- 0.55; xr <- 0.53; yt <- 0.955;
} else if(outers_only){
  xl <- 0.55; yb <- 0.025; xr <- 0.975; yt <- 0.45;
} else {
  xl <- 0.55; yb <- 0.025; xr <- 0.975; yt <- 0.45; 
}
x0 <- 0; y0 <- 0; x1 <- 0.1; y1 <- 0.1;
# segments(x0 = x0, y0 = y1, x1 = xl, y1 = yt, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y0, x1 = xr, y1 = yb, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y1, x1 = xr, y1 = yt, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x0, y0 = y0, x1 = xl, y1 = yb, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))

for(i in 1:length(d)){
  lines(d[[i]], col = cols[i], lwd =1.75* 1.1)
}

#make subplot on log-log scale
rect(x0, y0, x1, y1, lwd =1.75* 1, col = rgb(0,0,0,0.05))
bounds_of_subplot <- c(-1, -10)
for(i in 1:length(d)){
  pts2plot <- which(log10(d[[i]][,1]) < bounds_of_subplot[1] & log10(d[[i]][,1]) > bounds_of_subplot[2] & 
                      log10(d[[i]][,2]) < bounds_of_subplot[1] & log10(d[[i]][,2]) > bounds_of_subplot[2])
  xs <- log10(d[[i]][pts2plot,1])
  ys <- log10(d[[i]][pts2plot,2])
  #translate to origin
  xs <- xs - bounds_of_subplot[2]
  ys <- ys - bounds_of_subplot[2]
  #rescale to subplot
  xs <- xs / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (xr-xl)
  ys <- ys / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (yt-yb)
  #translate to subplot region
  xs <- xs + xl
  ys <- ys + yb
  #plot
  lines(xs, y = ys, col = cols[i], lwd =1.75* 2)
}

rect(xl, yb, xr, yt, lwd =1.75* 2, col = rgb(0,0,0,0.05))
segments(x0 = xl, y0 = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1), x1 = xl - 0.0075, lwd =1.75* 2,
         y1 = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1)) #horiz
segments(x0 = seq(xl, xr, length.out = abs(diff(bounds_of_subplot)) + 1), y0 = yt, 
         x1 = seq(xl, xr, length.out = abs(diff(bounds_of_subplot)) + 1), lwd =1.75* 2,
         y1 = yt + 0.00735) #vert
text(x = xl - 0.015, y = seq(yb, yt, length.out = abs(diff(bounds_of_subplot)) + 1) + 0.025, pos = 2, srt = 90, cex = 2.25,
     labels = c(sapply(bounds_of_subplot[2]:(bounds_of_subplot[1]), function(i) as.expression(bquote(10^ .(i)))))) #vert 
text(x = seq(xl, xr, length.out = abs(diff(bounds_of_subplot))+1), y = yt + 0.01, pos = 3, cex = 2.25,
     labels = (c(sapply(bounds_of_subplot[2]:(bounds_of_subplot[1]), function(i) as.expression(bquote(10^ .(i))))))) #horiz

# segments(0,0,1,1, lwd =1.75* 2, col = rgb(0,0,0,0.15))
text(x = 0.9, y = 0.895, labels = "1-to-1 line", pos = ifelse(inners_only, 3, 1), srt = 45, cex = 3)
rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, lwd =4)

#annotate the fancy graph w/ dot for point along line
points(x = po$prop_outliers_subpop, y = po$prop_outliers_mix, pch = 19, cex = 3.25)
points(x = po$prop_outliers_subpop, y = po$prop_outliers_mix, col = "#FAED27", pch = 19, cex = 2.5)
if(log10(po$prop_outliers_subpop) < bounds_of_subplot[1] & log10(po$prop_outliers_subpop) > bounds_of_subplot[2] & 
   log10(po$prop_outliers_mix) < bounds_of_subplot[1] & log10(po$prop_outliers_mix) > bounds_of_subplot[2]){

  xsp <- (log10(po$prop_outliers_subpop) - bounds_of_subplot[2]) / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (xr-xl) + xl
  ysp <- (log10(po$prop_outliers_mix) - bounds_of_subplot[2]) / -(bounds_of_subplot[2] - bounds_of_subplot[1]) * (yt-yb) + yb
  
  points(x = xsp, y = ysp, pch = 19, cex = 4.25)
  points(x = xsp, y = ysp, col = "#FAED27", pch = 19, cex = 3.25)
  
  
}

dev.off()
