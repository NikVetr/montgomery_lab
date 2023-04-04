
#### functions ####
intersect_intervals <- function(l1, l2){
  # create matrix for overlaps
  start <- cbind(c(l1[,1], l2[,1]), 1)
  end <- cbind(c(l1[,2], l2[, 2]), -1)
  over <- rbind(start, end)
  # order
  over <- over[order(over[,1]),]
  # get overlap count
  over <- cbind(over, overlap=cumsum(over[,2]))
  # create the overlap matrix
  inter <- cbind(L=over[(over[,2] == 1) & (over[,3] == 2), 1],
                 U=over[(over[,2] == -1) & (over[, 3] == 1), 1])
  inter
}

pnorm_bounds_old <- function(bounds, mu, sd, w){
  sapply(1:nrow(bounds), function(ri){
    pnorm(bounds[ri,"U"], mu, sd) - pnorm(bounds[ri,"L"], mu, sd)
  }) * w
}

pnorm_bounds <- function(bounds, mu, sd, w, return_log = F){
  sapply(1:nrow(bounds), function(ri){
    la <- pnorm(bounds[ri,"U"], mu, sd, log.p = T)
    lb <- pnorm(bounds[ri,"L"], mu, sd, log.p = T)
    ldiff <- (la + log(1 - exp(lb - la))) + log(w)
    if(return_log){
      return(ldiff)
    } else {
      return(exp(ldiff))
    }
  })
}

prop_outliers <- function(group_diff = 1, zscore = 2.5, sd_1 = 1, sd_2 = 1, w_1 = 0.5, stratify_direction = T, eps = 1E-15){
  
  #check for comparisons w/ zero or negative sd
  sd_1 <- max(c(sd_1, eps))
  sd_2 <- max(c(sd_2, eps))
  
  #turn high-level parameters into lower level parameters
  w_2 <- 1 - w_1 #weight of second normal, deterministically computed
  mu_1 <- -group_diff / (1 + w_1 / w_2) #the first normal's mean to one side (left if group_diff is +)
  mu_2 <- mu_1 + group_diff #the second normal's mean to the other side (right if group_diff is +)
  
  #from https://stats.stackexchange.com/questions/16608/w_hat-is-the-variance-of-the-w_eighted-mixture-of-tw_o-gaussians/16609#16609
  sd_3 <- sqrt(w_1*sd_1^2 + w_2*sd_2^2 + w_1*mu_1^2 + w_2*mu_2^2 - (w_1*mu_1 + w_2*mu_2)^2) #the standard deviation of the mixture distribution
  mu_3 <- mu_1 * w_1 + mu_2 * w_2
  
  #find upper and lower bounds in all three dists
  bounds_1 <- rbind(c(-Inf, mu_1 - sd_1 * zscore), c(mu_1 + sd_1 * zscore, Inf)); colnames(bounds_1) <- c("L", "U")
  bounds_2 <- rbind(c(-Inf, mu_2 - sd_2 * zscore), c(mu_2 + sd_2 * zscore, Inf)); colnames(bounds_2) <- c("L", "U")
  bounds_3 <- rbind(c(-Inf, mu_3 - sd_3 * zscore), c(mu_3 + sd_3 * zscore, Inf)); colnames(bounds_3) <- c("L", "U")
  
  #compute their overall integrals
  integ_sep <- c(sum(pnorm_bounds(bounds_1, mu_1, sd_1, w_1)), sum(pnorm_bounds(bounds_2, mu_2, sd_2, w_2)))
  integ_mix <- c(sum(pnorm_bounds(bounds_3, mu_1, sd_1, w_1)), sum(pnorm_bounds(bounds_3, mu_2, sd_2, w_2)))
  
  #find overlaps in bounds to decompose integrals
  integ_shared <- c(sum(pnorm_bounds(intersect_intervals(bounds_1, bounds_3), mu_1, sd_1, w_1)), 
                    sum(pnorm_bounds(intersect_intervals(bounds_2, bounds_3), mu_2, sd_2, w_2)))
  
  #compute changes to outlier prop from mixing
  integ_lost <- integ_sep - integ_shared #extra frac in sep not shared ij mix
  integ_gained <- integ_mix - integ_shared #extra frac in mix not shared in sep
  integ_net <- integ_gained - integ_lost
  
  #also relative to existing outlier proportion
  if(stratify_direction){
    prop_integ_lost <- sum(integ_gained) / sum(integ_mix)
    prop_integ_gained <- sum(integ_lost) / sum(integ_mix)
    prop_integ_net <- -sum(integ_net) / sum(integ_mix)
  } else {
    prop_integ_lost <- sum(integ_lost) / sum(integ_sep)
    prop_integ_gained <- sum(integ_gained) / sum(integ_sep)
    prop_integ_net <- sum(integ_net) / sum(integ_sep)
  }
  
  #find overall result and return -- specifically, individuals lost by mixing, not by stratifying (flip lost and gained in that case)
  if(stratify_direction){
    return(c(lost = sum(integ_gained), gained = sum(integ_lost), net = -sum(integ_net),
             p_lost = prop_integ_lost, p_gained = prop_integ_gained, p_net = prop_integ_net))
  } else {
    return(c(lost = sum(integ_lost), gained = sum(integ_gained), net = sum(integ_net),
             p_lost = prop_integ_lost, p_gained = prop_integ_gained, p_net = prop_integ_net))
  }
}

nice_ticks <- function(d, nt = 5, buffer = c(L = T, U = F)){
  d <- round(d, max(c(0, -floor(log10(diff(range(d)) / 1000)))))
  magvar <- 10^floor(log10(diff(range(d)) / nt))
  vals <- seq(from = min(d) - min(d) %% magvar, 
              to = max(d) - max(d) %% magvar, 
              by = magvar)
  tick_sub <- floor(length(vals) / ny)
  tick0 <-  which.min(abs(vals))
  new_vals <- vals[c(rev(seq(from = tick0, to = 1, by = -tick_sub)[-1]), 
                     tick0, 
                     seq(from = tick0, to = length(vals), by = tick_sub)[-1])]
  if(min(new_vals) > min(d) & buffer["L"]){
    new_vals <- c(new_vals[1] - diff(new_vals[1:2]), new_vals)
  }
  
  if(max(new_vals) < max(d) & buffer["U"]){
    new_vals <- c(new_vals, new_vals[length(new_vals)] + diff(new_vals[1:2]))
  }
  
  new_vals
}


#### obtain input data ####
zscore <- 2.5
dmu <- F
use_factor_changes <- F #if T, represent changes as factor differences
#so going from 1% outliers -> 2% outliers is not a change of 0.01, but of 2

if(dmu){
  delta_mus <- 0:400/100
  prop_gained_lost <- data.frame(t(sapply(delta_mus, function(dmu) prop_outliers(group_diff = dmu, zscore = zscore, sd_1 = 1, sd_2 = 1, w_1 = 0.5))))
  input_var <- delta_mus
} else {
  sds <- 0:400/100
  prop_gained_lost <- data.frame(t(sapply(sds, function(this_sd) prop_outliers(group_diff = 0, zscore = zscore, sd_1 = 1, sd_2 = this_sd, w_1 = 0.5))))
  input_var <- sds
}

if(use_factor_changes){
  d <- prop_gained_lost[,c("p_lost", "p_gained", "p_net")]
  d <- asinh(d) #could also do pseudolog?
  infinds <- unique(which(d == Inf | is.na(d), arr.ind = T)[,1])
  if(length(infinds) > 0){
    d <- d[-infinds,]
    input_var <- input_var[-infinds]  
  }
  colnames(d) <- gsub("p_", "", colnames(d))
} else {
  d <- prop_gained_lost[,c("lost", "gained", "net")]
}

nid <- length(input_var)

#simulate dummy empirical data for up top
empirical_obs <- rbeta(1E4, 1, 6) * diff(range(input_var)) + input_var[1]

#### plot the fig ####
par(mar = c(4,7,17,4), xpd = T)

#plot main plot
lwd <- 2
cols <- RColorBrewer::brewer.pal(3, "Dark2")
plot(input_var, d$net, type = "l", lwd = lwd, col = cols[1], 
     ylim = range(d), xlim = range(input_var) + c(0,1), 
     xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#get a few useful graphical params
lwh <- lwd / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
lww <- lwd / 96 / par("pin")[1] * (par("usr")[2] - par("usr")[1])

#figure axes
nx <- 5
xvals <- seq(0, max(input_var), by = floor(max(input_var) / nx / 10^floor(log10(max(input_var) / nx))) * 10^floor(log10(max(input_var) / nx)))

ny <- 6
if(use_factor_changes){
  yvals <- nice_ticks(d, ny)
  ytxt <- sinh(yvals)
  trailing_mag <- (10^(floor(log10(abs(ytxt)))-1) + 1E-30)
  ytxt[!is.na(trailing_mag)] <- ytxt[!is.na(trailing_mag)] - 
    ytxt[!is.na(trailing_mag)] %% trailing_mag[!is.na(trailing_mag)]
  yvals <- asinh(ytxt)
  if(any(log10(abs(ytxt)) > 2)){
    # ytxt <- gsub("e\\+|e", "\u00D710^", formatC(ytxt, digits = 1, format = "e"))
    ytxt <- gsub("e\\+|e", "E", formatC(ytxt, digits = 1, format = "e"))  
  }
  text(y = yvals, x = xvals[1] - diff(range(input_var)) / 20, labels = ytxt, pos = 2)
} else {
  yvals <- nice_ticks(d, ny)
  text(y = yvals, x = xvals[1] - diff(range(input_var)) / 20, labels = yvals, pos = 2)
  
}

#horiz axis
segments(x0 = xvals[1], x1 = max(xvals), y0 = yvals[1] - diff(range(d)) / 30, y1 = yvals[1] - diff(range(d)) / 30, lwd = 2)
segments(x0 = xvals, x1 = xvals, y0 = yvals[1] - diff(range(d)) / 30, y1 = yvals[1] - diff(range(d)) / 20, lwd = 2)
text(x = xvals, y = yvals[1] - diff(range(d)) / 20, lwd = 2, labels = xvals, pos = 1)
text(x = mean(range(xvals)), y = yvals[1] - diff(range(d)) / 8, lwd = 2, 
     labels = latex2exp::TeX(paste0("", ifelse(dmu, "\\Delta\\mu", "\\sigma$_2$ when \\sigma$_1$ = 1"))), pos = 1, cex = 1.5)

#vert axis
segments(y0 = yvals[1], y1 = max(yvals), x0 = xvals[1] - diff(range(input_var)) / 30, x1 = xvals[1] - diff(range(input_var)) / 30, lwd = 2)
segments(y0 = yvals, y1 = yvals, x0 = xvals[1] - diff(range(input_var)) / 30, x1 = xvals[1] - diff(range(input_var)) / 20, lwd = 2)
text(y = mean(range(yvals)), x = xvals[1] - diff(range(input_var)) / 5, lwd = 2, 
     labels = paste0("change in proportion of outliers\ngiven a z-score threshold of ", zscore), 
     pos = 3, srt = 90, cex = 1.25)

#label the lines
segments(y0 = d$gained[nid] + lwh / sqrt(2), y1 = d$gained[nid] + lwh / sqrt(2) + strheight(".") * 2,
         x0 = input_var[nid] + diff(range(input_var)) / 100, x1 = input_var[nid] + diff(range(input_var)) / 9,
         lty = 2, col = cols[2], lwd = 2)
text(labels = expression('proportion outliers\n', bold(gained)*  ' by stratifying'), 
     pos = 4, col = cols[2], x = input_var[nid] + diff(range(input_var)) / 10, y = d$gained[nid] + strheight(".") * 2, xpd = NA)

text(labels = expression('proportion outliers\n', bold(changed)*  ' by stratifying'), 
     pos = 4, col = cols[1], x = input_var[nid] + diff(range(input_var)) / 10, y = d$net[nid] - strheight(".") * 2, xpd = NA)
segments(y0 = d$net[nid] + lwh / sqrt(2), y1 = d$net[nid] + lwh / sqrt(2) - strheight(".") * 2,
         x0 = input_var[nid] + diff(range(input_var)) / 100, x1 = input_var[nid] + diff(range(input_var)) / 9,
         lty = 2, col = cols[1], lwd = 2)

text(labels = expression('proportion outliers\n', bold(lost)*  ' by stratifying'), 
     pos = 4, col = cols[3], x = input_var[nid] + diff(range(input_var)) / 10, y = d$lost[nid], xpd = NA)
segments(y0 = d$lost[nid], y1 = d$lost[nid],
         x0 = input_var[nid] + diff(range(input_var)) / 100, x1 = input_var[nid] + diff(range(input_var)) / 9,
         lty = 2, col = cols[3], lwd = 2)


#fake histogram at the top
nx2 = 10
histd <- hist(empirical_obs, breaks = seq(0, max(input_var), by = floor(max(input_var) / nx2 / 10^floor(log10(max(input_var) / nx2))) * 10^floor(log10(max(input_var) / nx2))), 
              plot = F)
histd$relative_density <- histd$density / max(histd$density)
rect(xleft = histd$breaks[1:nx2],
     xright = histd$breaks[2:(nx2+1)],
     ybottom = max(d) + diff(range(d)) / 10,
     ytop = max(d) + diff(range(d)) / 10 + histd$relative_density * diff(range(d)) / 1.25,
     col = "grey90")
ny2 <- 5
yvals2 <- seq(0, max(histd$counts), by = floor(max(histd$counts) / ny2 / 10^floor(log10(max(histd$counts) / ny2))) * 10^floor(log10(max(histd$counts) / ny2)))
yvals2 <- c(yvals2, yvals2[ny2+1] + diff(yvals2[1:2]))
yvals2_loc <- seq(max(d) + diff(range(d)) / 10, 
                  max(d) + diff(range(d)) / 10 + max(yvals2) / max(histd$counts) * diff(range(d)) / 1.25, 
                  length.out = length(yvals2))

segments(y0 = yvals2_loc[1], y1 = max(yvals2_loc), x0 = xvals[1] - diff(range(input_var)) / 30, x1 = xvals[1] - diff(range(input_var)) / 30, lwd = 2)
segments(y0 = yvals2_loc, y1 = yvals2_loc, x0 = xvals[1] - diff(range(input_var)) / 30, x1 = xvals[1] - diff(range(input_var)) / 20, lwd = 2)
text(y = yvals2_loc, x = xvals[1] - diff(range(input_var)) / 20, labels = yvals2, pos = 2)
text(y = mean(range(yvals2_loc)), x = xvals[1] - diff(range(input_var)) / 5, lwd = 2, labels = "# genes in GTEx corr.\nto given Hedges' g", pos = 3, srt = 90, cex = 1.25)

#connect axes
segments(x0 = histd$breaks, x1 = histd$breaks,
         y0 = yvals[1] - diff(range(d)) / 30,
         y1 = yvals2_loc[1], lty = 3, col = "grey50", lwd = 1)

#line for 0?
if(all(par("usr")[3:4] * c(-1,1) > diff(par("usr")[3:4]) / 4 )){
  segments(x0 = xvals[1] - diff(range(input_var)) / 30,
           x1 = xvals[length(xvals)],
           y0 = 0, y1 = 0, lwd = 2, lty = 3)
}

#plot main lines
lines(input_var, d$net, lwd = lwd, col = cols[1])
lines(input_var - lww / sqrt(2), d$gained + lwh / sqrt(2), lwd = lwd, col = cols[2])
lines(input_var, d$lost, lwd = lwd, col = cols[3])
