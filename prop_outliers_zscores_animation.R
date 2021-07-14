library(foreach)
library(doParallel)
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
    
    #symilar symmetry here as before w/ equal variances etc.
    
    #sum the probabilities (or return them as components)
    prop_outliers_subpop <- inner_prop1 + outer_prop1 + inner_prop2 + outer_prop2
    prop_outliers_mix <- inner_prop_mix1 + outer_prop_mix1 + middle_prop_mix1 + inner_prop_mix2 + outer_prop_mix2 + middle_prop_mix2
    
  }
  if(!returnAll){
    return(prop_outliers_mix)
  } else {
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
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  
  
  
  xy <- xy.coords(x,y)
  
  xo <- r*strwidth('A')
  
  yo <- r*strheight('A')
  
  for (i in theta) {
    
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
          
          labels, col=bg, ... )
    
  }
  
  text(xy$x, xy$y, labels, col=col, ... )
  
}
group_diffs <- seq(0, 4, length.out = 40)
outlier_thresholds <- 0:800/100 #in z-scores
# if(!exists("d")){
  d <- lapply(group_diffs, function(gd) t(sapply(outlier_thresholds, function(ot) c(
    indiv = (1 - pnorm(ot))*2,
    comb = prop_outliers(group_diff = gd, outlier_threshold = ot, range = range, sd1 = 0.5, sd2 = 0.5, numerical_integration = F)))))
# }

#let's parallelize these operations
# if(exists("cl")){stopCluster(cl = cl); rm(cl)}
if(!exists("cl")){
  cl <- makeCluster(16, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

renderVideo <- T
file.remove(paste0("~/Documents/zscore_outlier_animation/frames/", list.files(path = "~/Documents/zscore_outlier_animation/frames/", pattern = ".png")))
dir.create("~/Documents/zscore_outlier_animation/")
dir.create("~/Documents/zscore_outlier_animation/frames/")
framethin <- 1
nfps <- 60
cols <- viridis::plasma(n = length(d), end = 0.95) #color scheme to use

#phase #1, set the scene by combining pops

#timing of events in phase #1 of animation
#first number is start time, second is duration
makeD1 <- c(0,1) * nfps 
makeD2 <- c(sum(makeD1) / nfps,1) * nfps
raiseD2 <- c(sum(makeD2) / nfps,1.5) * nfps
makeD3 <- c(sum(raiseD2) / nfps,1) * nfps
dropD3 <- c(sum(makeD3) / nfps,1.5) * nfps
nsec_1 <- sum(dropD3) / nfps
nf <- (nsec_1 * nfps)

foreach(f1=seq(1, nf, framethin), .packages = c("latex2exp")) %dopar% {
  
cat(paste0(f1, " "))
png(filename = paste0("~/Documents/zscore_outlier_animation/frames/zscore_cartoon_", 
                      paste0(rep(0, 5-nchar(((f1 - 1) / framethin) + 1)), collapse = ""), ((f1 - 1) / framethin) + 1,".png"), 
    width = 2600, height = 1000, family="Arial Unicode MS")
layout((t(c(1,1,1,2,2))))
par(mar = c(4,0,0,2), xpd = F)
cols2 <- c(cols[1], cols[21], cols[11])
xlims <- c(-4,4); ylims <- c(0,2)
plot(1,1,xlim = xlims, ylim = ylims, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sd1 = 0.5
sd2 = 0.5
z = 2
mu1 = -0.75
mu2 = 0.75
range <- seq(mu1 - 4*sd1, mu2 + 4*sd2, by = 0.01)
segments(x0 = min(range), x1 = max(range), y0 = c(1 - min(max(0,(f1-dropD3[1])/dropD3[2]),1),1), lwd = 2, 
         col = c(grDevices::adjustcolor(1, alpha.f = min(max(0,(f1-makeD3[1])/makeD3[2]),1)), 
                 grDevices::adjustcolor(1, alpha.f = min(max(0,(f1-makeD1[1])/makeD1[2]),1))))
w1 <- w2 <- 0.5
mu3 <- mu1 * w1 + mu2 * w2
sd3 <- sqrt(w1*sd1^2 + w2*sd2^2 + w1*mu1^2 + w2*mu2^2 - (w1*mu1 + w2*mu2)^2) #the standard deviation of the mixture distribution
dist1 <- dnorm(range, mu1, sd1)
dist2 <- dnorm(range, mu2, sd2)
dist3 <- (dist1 + dist2)
polygon(x = range, y = dist1 + 1, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f1 - makeD1[1]) / makeD1[2]),0.2)), border = NA)
if(f1 > (dropD3[1]-1)){
  polygon(x = range, y = dist2 + 1, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f1 - makeD2[1]) / makeD2[2]),0.2)), border = NA)
} else {
  polygon(x = c(range, rev(range)), y = c(dist2 + dist1 * min(max(0,(f1-raiseD2[1])/raiseD2[2]),1) + 1, rev(dist1 * min(max(0,(f1-raiseD2[1])/raiseD2[2]),1) + 1)), 
          col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f1-makeD2[1])/makeD2[2]),0.2)), border = NA)
}
lines(range, dist1 + 1, lwd =1.75* 2, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f1 - makeD1[1]) / makeD1[2]),1)))
if(f1 > (dropD3[1]-1)){
  lines(range, dist2 + 1, lwd =1.75* 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f1 - makeD2[1]) / makeD2[2]),1)))
} else {
  lines(range, dist2 + dist1 * min(max(0,(f1-raiseD2[1])/raiseD2[2]),1) + 1, lwd =1.75* 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f1 - makeD2[1]) / makeD2[2]),1)))
}
polygon(x = range, y = dist3 + 1 - min(max(0,(f1-dropD3[1])/dropD3[2]),1), 
        col = grDevices::adjustcolor("white", alpha.f = min(max(0,(f1-makeD3[1])/makeD3[2]),1)), border = NA)
polygon(x = range, y = dist3 + 1 - min(max(0,(f1-dropD3[1])/dropD3[2]),1), 
        col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f1-makeD3[1])/makeD3[2]),0.2)), border = NA)
lines(range, dist3 + 1 - min(max(0,(f1-dropD3[1])/dropD3[2]),1), lwd =1.75* 2, 
      col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f1-makeD3[1])/makeD3[2]),1)))


dev.off()


}

f1 = nf #holdover from single-core implementation lol




#phase #2, plot all the extra info and vary it

#timing of events in phase #2 of animation
#first number is start time, second is duration
makeMus <- c(nsec_1,1.25) * nfps
makeZs <- c(sum(makeMus) / nfps,1.25) * nfps
makeExcl <- c(sum(makeZs) / nfps,1.25) * nfps
makeSPlot <- c(sum(makeExcl) / nfps,1.25) * nfps
changeZs <- c(sum(makeSPlot) / nfps,8) * nfps
changeMus <- c(sum(changeZs) / nfps, 8) * nfps
changeSDs <- c(sum(changeMus) / nfps, 8) * nfps
changeEverything <- c(sum(changeSDs) / nfps, 13) * nfps

#here is the range of model parameters to be varied
zs <- c(rep(2, changeZs[1]), seq(2.001,1,length.out = changeZs[2] / 6), seq(1.001,4,length.out = changeZs[2] * 3 / 6), seq(4.001,2.001,length.out = changeZs[2] * 2 / 6), rep(2, 5000))
mu1s <- c(rep(-0.75, changeMus[1]), seq(-0.7501,0.25,length.out = changeMus[2] / 6), seq(0.25001,-2.25,length.out = changeMus[2] * 3 / 6), seq(-2.25001,-0.75001,length.out = changeMus[2] * 2 / 6), rep(-0.75, 5000))
sd1s <- c(rep(0.5, changeSDs[1]), seq(0.5001, 1, length.out = changeSDs[2] / 3), seq(1.001, 0.25, length.out = changeSDs[2] / 2), seq(0.25001, 0.5001, length.out = changeSDs[2] / 6), rep(0.5, 500))
zs[changeEverything[1]:(sum(changeEverything)-1)] <- c(seq(2,1.25,length.out = changeEverything[2] / 4), seq(1.25001,2.5,length.out = changeEverything[2] / 2), seq(2.5001,2,length.out = changeEverything[2] / 4))
mu1s[changeEverything[1]:(sum(changeEverything)-1)] <- c(seq(-0.75,-1.25,length.out = changeEverything[2] / 2), seq(-1.25001,-0.25,length.out = changeEverything[2] / 4), seq(-0.25001,-0.75,length.out = changeEverything[2] / 4))
sd1s[changeEverything[1]:(sum(changeEverything)-1)] <- c(seq(0.5,1,length.out = changeEverything[2] / 4), seq(1.001,0.25,length.out = changeEverything[2] / 4), seq(0.25001,0.5,length.out = changeEverything[2] / 2))

#iterate through images
foreach(f2=seq(f1+1, sum(changeEverything), framethin), .packages = c("latex2exp")) %dopar% {
# for(f2 in 2038){ #for troubleshooting one image at a time
cat(paste0(f2, " "))
  
png(filename = paste0("~/Documents/zscore_outlier_animation/frames/zscore_cartoon_", paste0(rep(0, 5-nchar(((f2 - 1) / framethin) + 1)), collapse = ""),((f2 - 1) / framethin) + 1,".png"), width = 2600, height = 1000, family="Arial Unicode MS")
layout((t(c(1,1,1,2,2))))
par(mar = c(4,0,0,2), xpd = F)
cols2 <- c(cols[1], cols[21], cols[11])
xlims <- c(-4,4); ylims <- c(0,2)
plot(1,1,xlim = xlims, ylim = ylims, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sd1 = sd1s[f2]
sd2 = 0.5
z = zs[f2]
mu1 = mu1s[f2]
mu2 = 0.75
if(sd1s[f2]!= sd1s[f2-1] | f2 == (f1+1)){
  d <- lapply(group_diffs, function(gd) t(sapply(outlier_thresholds, function(ot) c(
    indiv = (1 - pnorm(ot))*2,
    comb = prop_outliers(group_diff = gd, outlier_threshold = ot, range = range, sd1 = sd1, sd2 = sd2, numerical_integration = F)))))
}
range <- seq(mu1 - 4*sd1, mu2 + 4*sd2, by = 0.01)
segments(x0 = min(range), x1 = max(range), y0 = c(0,1), lwd = 2)
w1 <- w2 <- 0.5
mu3 <- mu1 * w1 + mu2 * w2
sd3 <- sqrt(w1*sd1^2 + w2*sd2^2 + w1*mu1^2 + w2*mu2^2 - (w1*mu1 + w2*mu2)^2) #the standard deviation of the mixture distribution
dist1 <- dnorm(range, mu1, sd1)
dist2 <- dnorm(range, mu2, sd2)
dist3 <- (dist1 + dist2)
aboveDistRat <- 0.825 / max(dist3)
if(max(dist3) > 0.825){
  dist3 <- dist3 / max(dist3) * 0.825
  plot_warning <- T
} else{
  plot_warning <- F
}
if(max(dist1) > 0.8 | max(dist2) > 0.8){
  to_div_by <- max(c(dist1,dist2))
  dist1 <- dist1 / to_div_by * 0.8
  dist2 <- dist2 / to_div_by * 0.8
}
lines(range, dist1 + 1, lwd =1.75* 2, col = cols2[1])
lines(range, dist2 + 1, lwd =1.75* 2, col = cols2[2])
lines(range, dist3, lwd =1.75* 2, col = cols2[3])
polygon(x = range, y = dist1 + 1, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.2))
polygon(x = range, y = dist2 + 1, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.2))
polygon(x = range, y = dist3, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.2))
text(labels = latex2exp::TeX("$\\mu_{1}$"), pos = 3, x = range[which.max(dist1)], y = max(dist1) + 1, cex =4, font = 2, 
     col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))
text(labels = latex2exp::TeX("$\\mu_{2}$"), pos = 3, x = range[which.max(dist2)], y = max(dist2) + 1, cex =4, font = 2, 
     col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))
text(labels = latex2exp::TeX("$\\mu_{mix}$"), pos = 3, x = mu1 / 2 + mu2 / 2, y = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, cex =4, font = 2, 
     col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))

#bounds for dist1
segments(x0 = mu1 + sd1 * z, x1 = mu1 + sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 + sd1 * z), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
segments(x0 = mu1 - sd1 * z, x1 = mu1 - sd1 * z, y0 = 1, y1 = 1 + dist1[order(abs(range - mu1 - sd1 * z), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu1 - sd1 * z], mu1 - sd1 * z), y = c(dist1[range < mu1 - sd1 * z] + 1, 1), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu1 + sd1 * z], mu1 + sd1 * z), y = c(dist1[range > mu1 + sd1 * z] + 1, 1), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)

#bounds for dist2
segments(x0 = mu2 + sd2 * z, x1 = mu2 + sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 + sd2 * z), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
segments(x0 = mu2 - sd2 * z, x1 = mu2 - sd2 * z, y0 = 1, y1 = 1 + dist2[order(abs(range - mu2 - sd2 * z), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu2 - sd2 * z], mu2 - sd2 * z), y = c(dist2[range < mu2 - sd2 * z] + 1, 1), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu2 + sd2 * z], mu2 + sd2 * z), y = c(dist2[range > mu2 + sd2 * z] + 1, 1), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)

#bounds for dist3
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y0 = 0, y1 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), lwd =1.75* 2, lty = 2)
polygon(x = c(range[range < mu3 - sd3 * z], mu3 - sd3 * z), y = c(dist3[range < mu3 - sd3 * z], 0), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)
polygon(x = c(range[range > mu3 + sd3 * z], mu3 + sd3 * z), y = c(dist3[range > mu3 + sd3 * z], 0), col = grDevices::adjustcolor(grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), alpha.f = 0.2), border = NA)

#mark means with lines?
segments(x0 = mu1, x1 = mu1, y0 = dist2[order(abs(range - mu1), decreasing = F)[1]] + 1, y1 = max(dist1) + 1, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)), lwd =1.75* 1, lty = 3)
segments(x0 = mu2, x1 = mu2, y0 = dist1[order(abs(range - mu2), decreasing = F)[1]] + 1, y1 = max(dist2) + 1, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)), lwd =1.75* 1, lty = 3)
segments(x0 = mu3, x1 = mu3, y0 = 0, y1 = dist3[order(abs((mu1 / 2 + mu2 / 2) - range), decreasing = F)[1]] + 0.1, col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)), lwd =1.75* 1, lty = 3)

#labels for bounds
oid <- max(min((abs((mu1 + sd1*z) - (mu2 - sd2*z)) - 0.5 - ifelse(((round(z, 1) * 10) %% 10) != 0, 0.15, 0)) / 1.5, 0), -0.055) #overlapping_inner_disp
# strwidth(latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), font = 2, units = "figure") * par("pin")[1] / par("fin")[1] * diff(xlims) / 2
text(labels = latex2exp::TeX(paste0("$\\mu_{1} - $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 - sd1*z, y = 0.9875, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))
text(labels = latex2exp::TeX(paste0("$\\mu_{1} + $", round(z, 1),"$\\sigma_{1}$")), pos = 1, x = mu1 + sd1*z, y = 0.9875 + oid, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))
text(labels = latex2exp::TeX(paste0("$\\mu_{2} - $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 - sd2*z, y = 0.9875, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))
text(labels = latex2exp::TeX(paste0("$\\mu_{2} + $", round(z, 1),"$\\sigma_{2}$")), pos = 1, x = mu2 + sd2*z, y = 0.9875, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} - $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 - sd3*z, y = -0.0125, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))
text(labels = latex2exp::TeX(paste0("$\\mu_{mix} + $", round(z, 1),"$\\sigma_{mix}$")), pos = 1, x = mu3 + sd3*z, y = -0.0125, cex = 2.75, font = 2, col = grDevices::adjustcolor(cols2[3], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)))

#add some bars to separate the mus
# segments(x0 = mu1, x1 = mu2, y0 = max(dist1) + 1.15, y1 = max(dist2) + 1.15, lwd =1.75* 1.5)
# segments(x0 = mu1, x1 = mu1, y0 = max(dist1) + 1.125, y1 = max(dist1) + 1.175, lwd =1.75* 1.5)
# segments(x0 = mu2, x1 = mu2, y0 = max(dist2) + 1.125, y1 = max(dist2) + 1.175, lwd =1.75* 1.5)
# text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = (max(dist2) + max(dist1)) / 2 + 1.1275,
#      srt = atan(((max(dist2) - max(dist1)) / (diff(ylims)) * par("pin")[2]) / ((mu2 - mu1) / diff(xlims) * par("pin")[1])) * 180 / pi)
#actually, make it a horizontal bar only
segments(x0 = mu1, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.15, y1 = max(c(dist1, dist2)) + 1.15, lwd =1.75* 1.5, col = grDevices::adjustcolor(1, alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))
segments(x0 = mu1, x1 = mu1, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd =1.75* 1.5, col = grDevices::adjustcolor(1, alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))
segments(x0 = mu2, x1 = mu2, y0 = max(c(dist1, dist2)) + 1.125, y1 = max(c(dist1, dist2)) + 1.175, lwd =1.75* 1.5, col = grDevices::adjustcolor(1, alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))
text(latex2exp::TeX("$\\Delta_{\\mu} = \\mu_{2} - \\mu_{1}$"), x = (mu2 + mu1) / 2, pos = 3, y = max(c(dist1, dist2)) + 1.145,
     srt = 0, cex = 3, col = grDevices::adjustcolor(1, alpha.f = min(max(0,(f2-makeMus[1])/makeMus[2]),1)))

#show inclusion / exclusion zone
segments(x0 = mu3 + sd3 * z, x1 = mu3 + sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 + sd3 * z)), decreasing = F)[1]], col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), lwd =1.75* 2, lty = 2)
segments(x0 = mu3 - sd3 * z, x1 = mu3 - sd3 * z, y1 = 0.9, y0 = dist3[order(abs(range - (mu3 - sd3 * z)), decreasing = F)[1]], col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), lwd =1.75* 2, lty = 2)
text(x = mu3 + sd3 * z - 0.025, y = 0.65, labels = "excluded", srt = 90, pos = 2, col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), cex = 2.5)
text(x = mu3 + sd3 * z + 0.025, y = 0.65, labels = "included", srt = 270, pos = 4, col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), cex = 2.5)
text(x = mu3 - sd3 * z - 0.025, y = 0.65, labels = "included", srt = 90, pos = 2, col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), cex = 2.5)
text(x = mu3- sd3 * z + 0.025, y = 0.65, labels = "excluded", srt = 270, pos = 4, col = grDevices::adjustcolor("grey50", alpha.f = min(max(0,(f2-makeExcl[1])/makeExcl[2]),1)), cex = 2.5)

#label regions of subpopulation kdes
text(x = mu1 + sd1 * z  - 0.025, y = 1.05, labels = "m", srt = 90, pos = 2, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu1 + sd1 * z  + 0.025, y = 1.05, labels = "in", srt = 270, pos = 4, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu1 - sd1 * z  - 0.025, y = 1.05, labels = "ot", srt = 90, pos = 2, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu1- sd1 * z  + 0.025, y = 1.05, labels = "m", srt = 270, pos = 4, col = grDevices::adjustcolor(cols2[1], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu2 + sd2 * z  - 0.025, y = 1.05, labels = "m", srt = 90, pos = 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu2 + sd2 * z  + 0.025, y = 1.05, labels = "ot", srt = 270, pos = 4, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu2 - sd2 * z  - 0.025, y = 1.05, labels = "in", srt = 90, pos = 2, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)
text(x = mu2- sd2 * z  + 0.025, y = 1.05, labels = "m", srt = 270, pos = 4, col = grDevices::adjustcolor(cols2[2], alpha.f = min(max(0,(f2-makeZs[1])/makeZs[2]),1)), cex =2)

#draw little bargraph for tail proportions
po <- prop_outliers(group_diff = mu2-mu1, outlier_threshold = z, sd1 = sd1, sd2 = sd2, returnAll = T)
addX <- -0.2 + max((mu3 + z*sd3) - 3.035, 0)
addY <- 0 + max((mu3 + z*sd3) - 2.5, 0) / 8
segments(x0 = 3.5 + addX, x1 = 4 + addX, y0 = 0 + addY, y1 = 0 + addY, lwd=1.75* 2)
rect(xleft = 3.5 + addX, xright = 3.75 + addX, ybottom = 0 + addY, ytop = (po$inner_prop1 + po$outer_prop1) / 0.25 + addY, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.36), lwd = 2)
rect(lwd = 2, xleft = 3.5 + addX, xright = 3.75 + addX, ybottom = (po$inner_prop1 + po$outer_prop1) / 0.25 + addY, 
     ytop = (po$inner_prop1 + po$outer_prop1) / 0.25 + (po$inner_prop2 + po$outer_prop2) / 0.25 + addY, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.36))
rect(xleft = 3.75 + addX, xright = 4 + addX, ybottom = 0 + addY, ytop = po$prop_outliers_mix / 0.25 + addY, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.36), lwd = 2)
# text(x = 3.75, pos = 4, y = 1, labels = expression(integral()))
# addImg(png::readPNG("~/Documents/integral_symbol.png"), x = 3.65, y = -0.05, width = 0.075) #ehh maybe not
text(x = 3.625 + addX, pos = 1, y = -0.00625 + addY, labels = "sub", cex = 2.5)
text(x = 3.875 + addX, pos = 1, y = -0.00625 + addY, labels = "mix", cex = 2.5)
text(x = 3.75 + addX, pos = 1, y = -0.05 + addY, labels = "integrals", cex = 2.5, font = 4, col = "darkred", xpd = NA)
rect(xleft = 3.49 + addX, xright = 4.01 + addX, yb = -0.4 + addY, yt = 1 + addY, col = grDevices::adjustcolor("white", alpha.f = 1 - min(max(0,(f2-makeSPlot[1])/makeSPlot[2]),1)), border = NA, xpd = NA)

#put title for what's a-changing
if(sd1s[f2]!= sd1s[f2-1] | mu1s[f2]!= mu1s[f2-1] | zs[f2]!= zs[f2-1] | f2 > sum(makeSPlot)){
  text(x = -3.375, y = 1.9125, labels = "now changing", xpd = NA, srt = 20, cex = 4)
  if(sd1s[f2]!= sd1s[f2-1] & mu1s[f2]== mu1s[f2-1] & zs[f2]== zs[f2-1]){
    text(x = -3.2, y = 1.85, labels = "population variances", xpd = NA, srt = 20, cex = 5, font = 2, col = "darkred")
  } else if(sd1s[f2]== sd1s[f2-1] & mu1s[f2]!= mu1s[f2-1] & zs[f2]== zs[f2-1]){
    text(x = -3.2, y = 1.85, labels = "population means", xpd = NA, srt = 20, cex = 5, font = 2, col = "darkred")
  } else if(sd1s[f2] == sd1s[f2-1] & mu1s[f2] == mu1s[f2-1] & zs[f2]!= zs[f2-1]){
    text(x = -3.2, y = 1.85, labels = "exclusion threshold", xpd = NA, srt = 20, cex = 5, font = 2, col = "darkred")
  } else if((sd1s[f2]!= sd1s[f2-1] & mu1s[f2]!= mu1s[f2-1] & zs[f2]!= zs[f2-1]) & f2 > changeEverything[1] | f2 > changeEverything[1]){
    cols3 <- c(cols[seq(1,length(cols), by = 2)], rev(cols[seq(1,length(cols), by = 2)])[-1])
    shadowtext(x = -4.2, y = 1.65, labels = "EVERYTHING", xpd = NA, srt = 20, cex = 8, font = 2, col = cols3[(f2 %% length(cols3)) %% length(cols3)], bg = "black", r = 0.5, pos = 4)
    for(char in 1:(nchar("EVERYTHING")-1)){
      text(x = -4.2, y = 1.65, labels = substr("EVERYTHING", 1, (nchar("EVERYTHING") - char)), xpd = NA, srt = 20, 
       cex = 8, font = 2, col = cols3[(f2 %% length(cols3) + char) %% length(cols3)], pos = 4)
    }
  }
}

#make left hand side be white
rect(xleft = -4.5, xright = -4.295, ybottom = -1, ytop = 1.2, col = "white", border = NA)

#plot warning
if(plot_warning){
  ## tall warning
  # rect(xleft = -5, xright = -3.225, ybottom = -2, ytop = 0.5, xpd = NA, lwd = 3, col = "white", border = NA)
  # rect(xleft = -4.3, xright = -3.25, ybottom = -0.16, ytop = 0.02, xpd = NA, lwd = 3, col = "grey92")
  # rect(xleft = -4.3, xright = -3.25, ybottom = -0.16, ytop = -0.05, xpd = NA, lwd = 1.5, col = "white")
  # text("objects in figure may be \nlarger than they appear",x = -3.7875, y = -0.15, pos = 3, cex = 2, xpd = NA)
  # shadowtext(labels = "WARNING:", x = -3.7875, y = -0.045, pos = 3, cex = 3, xpd = NA, col = "red", bg = "black", r = 0.2)
  # long warning
  rect(xleft = -4.3, xright = -2.425, ybottom = -0.16, ytop = -0.08, xpd = NA, lwd = 2.5, col = "grey98")
  rect(xleft = -4.3, xright = -3.41, ybottom = -0.16, ytop = -0.08, xpd = NA, lwd = 3, col = "#fdca17")
  text("Axes Rescaled",x = -3.3975, y = -0.1275, pos = 4, cex = 3, xpd = NA, font = 3, col = "black")
  text(labels = "WARNING:", x = -4.15, y = -0.1275, pos = 4, cex = 3, xpd = NA, col = "black", font = 2)
  points(x = -4.2, y = -0.1275, pch = 17, cex = 5.5, xpd = NA)
  points(x = -4.2, y = -0.1275, pch = 17, cex = 4, xpd = NA, col = "#fdca17")
  points(x = -4.2, y = -0.1275, pch = "!", cex = 2.5, xpd = NA)
}

#plot grey version of warning
rect(xleft = -4.3, xright = -2.425, ybottom = -0.16, ytop = -0.08, xpd = NA, lwd = 2.5, col = grDevices::adjustcolor("grey95", alpha.f = min(aboveDistRat^16, 1)),
     border = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)))
rect(xleft = -4.3, xright = -3.41, ybottom = -0.16, ytop = -0.08, xpd = NA, lwd = 3, col = grDevices::adjustcolor("grey95", alpha.f = min(aboveDistRat^16, 1)),
     border = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)))
text("Axes Rescaled",x = -3.3975, y = -0.1275, pos = 4, cex = 3, xpd = NA, font = 3, col = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)))
text(labels = "WARNING:", x = -4.15, y = -0.1275, pos = 4, cex = 3, xpd = NA, col = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)), font = 2)
points(x = -4.2, y = -0.1275, pch = 17, cex = 5.5, xpd = NA, col = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)))
points(x = -4.2, y = -0.1275, pch = 17, cex = 4, xpd = NA, col = grDevices::adjustcolor("grey95", alpha.f = min(aboveDistRat^16, 1)))
points(x = -4.2, y = -0.1275, pch = "!", cex = 2.5, xpd = NA, col = grDevices::adjustcolor("grey60", alpha.f = min(aboveDistRat^16, 1)))
# fade in at appropriate time
rect(xleft = -4.4, xright = -2.4, ybottom = -0.2, ytop = -0.07, xpd = NA, lwd = 2.5, col = grDevices::adjustcolor("white", alpha.f = 1 - min(max(0,(f2-makeSPlot[1])/makeSPlot[2]),1)),
     border = grDevices::adjustcolor("white", alpha.f = 1 - min(max(0,(f2-makeSPlot[1])/makeSPlot[2]),1)))

#now for the second part, my fancy feather figure!
par(mar = c(6,3,0,1), xpd = F)

plot(1,1,xlim = c(0,1), ylim = c(0,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
# mtext(text = "Subpopulation Proportion Outside Given Z-Score", cex =1.75* 1.75, line = 4, side = 1) #horiz axis label
text(x = 0.055, y = -0.08, labels = "Subpopulation Proportion Outside Given Z-Score (     +     )", cex = 4, srt = 0, xpd = NA, pos = 4) #horiz axis label

# mtext(text = "Mixture Proportion Outside Given Z-Score", cex =1.75* 1.75, line = 2.5, side = 2) #vert axis label
text(x = -0.08, y = 0.125, labels = "Mixture Proportion Outside Given Z-Score (    )", cex = 4, srt = 90, xpd = NA, pos = 4) #vert axis label
points(x = -0.075, y = 0.84225, cex = 9, pch = 15, col = grDevices::adjustcolor(cols2[3], alpha.f = 0.36), xpd = NA)
points(x = -0.075, y = 0.84225, cex = 9.2, pch = 22, col = "black", xpd = NA)
axis(1, at = 0:10/10, labels = rep("", 11), lwd =1.75* 2, cex.axis = 2, tck = -0.015, line = -3.7)
points(x = 0.8475, y = -0.0795, cex = 9, pch = 15, col = grDevices::adjustcolor(cols2[1], alpha.f = 0.36), xpd = NA)
points(x = 0.8475, y = -0.0795, cex = 9.2, pch = 22, col = "black", xpd = NA)
points(x = 0.915, y = -0.0795, cex = 9, pch = 15, col = grDevices::adjustcolor(cols2[2], alpha.f = 0.36), xpd = NA)
points(x = 0.915, y = -0.0795, cex = 9.2, pch = 22, col = "black", xpd = NA)

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
text(labels = round(seq(from = range(group_diffs)[1], to = range(group_diffs)[2], length.out = diff(range(group_diffs)) + 1)), 
     y = seq(yb, yt - 0.01, length.out = diff(range(group_diffs)) + 1), x = xr, pos = 4, las=2, cex=1.75*1.5)
text(labels = latex2exp::TeX("$\\Delta_{\\mu$"), pos = 4, x = xl - 0.0015, y = yt + 0.015, cex =4, font = 2)

#annotate the fancy graph w/ arrow for delta mu
arrows(x0 = 0.125, x1 = 0.0875, y0 = (mu2-mu1) * (yt-yb) / 8 + yb, y1 = (mu2-mu1) * (yt-yb) / 8 + yb, lwd = 5, col = "darkred")
text(x = 0.125, pos = 4, y = (mu2-mu1) * (yt-yb) / 8 + yb, cex = 3, font = 3, labels = "currently", col = "darkred")

#subplot corners & connecting lines
xl <- 0.55 + abs(sd1 - sd2) / 15; yb <- 0.025; xr <- 0.975; yt <- 0.45 - abs(sd1 - sd2) / 15;
x0 <- 0; y0 <- 0; x1 <- 0.1; y1 <- 0.1;
# segments(x0 = x0, y0 = y1, x1 = xl, y1 = yt, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y0, x1 = xr, y1 = yb, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x1, y0 = y1, x1 = xr, y1 = yt, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))
# segments(x0 = x0, y0 = y0, x1 = xl, y1 = yb, lty = 2, lwd =1.75* 1, col = rgb(0,0,0,0.4))

#plot the actual lines
for(i in 1:length(d)){
  lines(d[[i]], col = cols[i], lwd =1.75* 1.1)
}

#make subplot on log-log scale
rect(x0, y0, x1, y1, lwd =1.75* 1, col = rgb(0,0,0,0.05))
bounds_of_subplot <- c(-1, -10)
# rect(xl-0.04, yb, xr, yt+0.04, border = NA, col = "white") #white frame for unequal variances
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
# text(x = 0.9, y = 0.895, labels = "1-to-1 line", pos = 1, srt = 45, cex = 3)
shadowtext(x = 0.9, y = 0.895, labels = "1-to-1 line", pos = 1, srt = 45, cex = 3, col = 1, bg = "white", r = 0.15)
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

#legend for dot
points(x = 0.86125, y = 1.0125, pch = 19, cex = 3.25)
points(x = 0.86125, y = 1.0125, col = "#FAED27", pch = 19, cex = 2.5)
text(1, 1.01125, "= where we at", pos = 2, cex = 2.25)

#fade in side plot
rect(xleft = -0.15, xright = 1.5, yb = -0.75, yt = 1.5, col = grDevices::adjustcolor("white", alpha.f = 1 - min(max(0,(f2-makeSPlot[1])/makeSPlot[2]),1)), border = NA, xpd = NA)

dev.off()
}

# system("cd Documents/zscore_outlier_animation; zip -r frames.zip frames") #zip the frames for export
# now actually render the video!
if(renderVideo){
  if(file.exists("~/Documents/zscore_outlier_animation/zscore_animation.mp4")){file.remove("~/Documents/zscore_outlier_animation/zscore_animation.mp4")}
  system(paste0("cd Documents/zscore_outlier_animation; ffmpeg -r ", nfps / framethin," -f image2 -s 2600x1000 -i frames/zscore_cartoon_%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p zscore_animation.mp4"))
}
