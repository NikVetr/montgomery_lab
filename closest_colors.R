library(farver)

#### specify parameters ####
color_space     <-  c("hsl", "lab", "lch", "oklab", "oklch")[4]
colorblind_safe <- T #weights over deutan, protan, tritan vision
colorblind_weights <- c(none = 6,
                        deutan = 6,
                        protan = 2,
                        tritan = 0.1)
brewer_cols <- T
rainbow_circle_cols <- F
constrain_to_similar_colors <- T
fb_rat <- 0.9 #how much wider should the constraint be than the quantiles?
user_cols <- NA
user_cols <- c(Transcriptomics = "#4477AA",
               Metabolomics = "#6D4B08",
               Proteomics = "#228833",
               Phosphoproteomics = "#F3A02B",
               ATAC = "#882255",
               Methylation = "#D687B5")
n_rainbow = 6
n_optim_runs <- 100
n_cols_to_add <- 1
colorwheel_space <- c("hsl", "lch", "oklch")[1]
plot_cbs <- c("none", "deutan", "protan", "tritan")[1:4]

#### functions ####
polar2cart <- function(t, r){c(r*cos(t), r * sin(t))}

logit <- function(p){log(p/(1-p))}

invlogit <- function(x){exp(x)/(1+exp(x))}

cyl2cart <- function(t,r,z){c(polar2cart(t,r),z)}

cart2polar <- function(x,y){c(atan(y/x), sqrt(x^2 + y^2))}

shrink <- function(x, a = 0.1){c(x[1] + diff(x)/2*a, x[2] - diff(x)/2*a)}

cvd_convert <- function(cols, severity = 1, linear = T, 
                        type = c("none", "deutan", "protan", "tritan")[1]){
  if(type == "deutan"){
    return(colorspace::deutan(cols, severity = severity, linear = linear))
  } else if(type == "protan"){
    return(colorspace::protan(cols, severity = severity, linear = linear))
  } else if(type == "tritan"){
    return(colorspace::tritan(cols, severity = severity, linear = linear))
  } else {
    return(cols)
  }
}

#set bounds across channels
flex_bounds <- function(x, xlim = c(0,1), rat = 1.1){
  range_x <- diff(range(x))
  mean_x <- mean(x)
  new_x <- mean_x + c(-1,1) / 2 * range_x * rat 
  new_x[1] <- max(xlim[1], new_x[1])
  new_x[2] <- min(xlim[2], new_x[2])
  return(new_x)
}

rescale <- function(M) {
  for (cni in cn)
    M[, cni] <- M[, cni] *
      (cs_ranges[[color_space]]$max[cni] -
         cs_ranges[[color_space]]$min[cni]) +
      cs_ranges[[color_space]]$min[cni]
  M
}


#specify objective function
mean_dist_func <- function(par,
                           curr_cols,
                           bounds_sc, 
                           bounds_l,
                           color_space,
                           cs_ranges,
                           dist_alg = "cie2000",
                           colorblind_safe = TRUE,
                           colorblind_weights = c(none = 6,
                                                  deutan = 6,
                                                  protan = 2,
                                                  tritan = 0.1),
                           init_weights = NULL,
                           init_colors = NULL,
                           return_info = F) {
  
  cn <- colnames(curr_cols)
  
  # transform params to [0,1]
  m <- matrix(par, ncol = length(cn), byrow = TRUE)
  nncols <- nrow(m)
  colnames(m) <- cn
  
  #for luminosity
  if ("l" %in% cn){
    if(nrow(m) > 1) {
      m[2:nrow(m), "l"] <- exp(m[2:nrow(m), "l"])
      m[, "l"]          <- cumsum(m[, "l"])
    } 
    m[, "l"] <- plogis(m[, "l"]) * diff(bounds_l) + bounds_l[1]
  }
  
  #for chroma / saturation
  if (any(cn %in% c("s","c"))){
    m[, intersect(c("s","c"), cn)] <-
      plogis(m[, intersect(c("s","c"), cn)]) *
      diff(bounds_sc) + bounds_sc[1]
  }
  
  #for remaining channels
  other_channels <- !(cn %in% c("s","c", "l"))
  if (any(other_channels)){
    m[, other_channels] <- plogis(m[, other_channels])
  }
  
  #then convert back to real units and recover hex codes
  new_hex <- encode_colour(rescale(m), from = color_space)
  curr_hex <- encode_colour(rescale(curr_cols), from = color_space)
  
  #convert to colorblind-space and compute pairwise distances
  if (colorblind_safe){
    cvd_states <- c("deutan","protan","tritan","none")
  } else {
    cvd_states <- "none"
  }  
  
  d <- setNames(numeric(length(cvd_states)), cvd_states)
  for (cvd_state in cvd_states) {
    
    #simulate colorblindness
    nh <- cvd_convert(new_hex, type = cvd_state)
    ch <- cvd_convert(curr_hex, type = cvd_state)
    
    #get back into appropriate colorspace
    cvals <- farver::decode_colour(ch, to = color_space)
    nvals <- farver::decode_colour(nh, to = color_space)
    
    #find distances between colors
    dists <- c(
      #old colors to new colors
      compare_colour(from = cvals, to = nvals, 
                     method = dist_alg, 
                     from_space = color_space),
      #new colors to new colors
      compare_colour(from = nvals, 
                     method = dist_alg, 
                     from_space = color_space)[upper.tri(diag(length(nh)))]
    )
    
    #compute harmonic mean of distances
    hmean_dist <- 1 / mean(1/dists)
    d[[cvd_state]] <- hmean_dist
    
  }
  
  #compute weighted harmonic mean
  wd <- sum(d * colorblind_weights[cvd_states])
  
  #also compute a penalty term so values can't rocket away
  penalty <- sum(par^2)
  
  #are we trying to also stick to some initial colors?
  if(!is.null(init_weights) & !is.null(init_colors)){
    init_cvals <- farver::decode_colour(init_colors, to = color_space)
    dists_to_init <- sapply(1:nncols, function(ci){
      compare_colour(from = init_cvals[ci,], to = nvals[ci,], 
                     method = dist_alg, 
                     from_space = color_space)})
    weighted_dists_to_init <- sum(dists_to_init * init_weights)
    penalty <- penalty + weighted_dists_to_init
  }
  
  #do we want to return the objective, or other relevant metadata?
  if(return_info){
    out <- list(wd = wd,
                new_hex = new_hex)
  } else {
    #since optim minimizes but we want to maximize the distance, multiply by -1
    out <- -wd + penalty
  }
  
  return(out)
  
}

#### generate colors ####
#generate colors in specified colorspace
if(brewer_cols){
  cols <- RColorBrewer::brewer.pal(8, "Set1")
  dat <- farver::decode_colour(cols, to = color_space)
} 
if(rainbow_circle_cols){
  dat <- do.call(rbind, lapply(1:n_rainbow, function(x) 
    farver::decode_colour("red", to = color_space)))
  dat[,"h"] <- seq(0, 360, length.out = n_rainbow + 1)[-1]
  cols <- encode_colour(dat, from = color_space)
}
if(!all(is.na(user_cols))){
  cols <- user_cols
  dat <- farver::decode_colour(cols, to = color_space)
}

#rescale colorspaces to be out of 1
cs_ranges <- list(
  hsl = list(
    min = c(h = 0, s = 0, l = 0),
    max = c(h = 360, s = 100, l = 100)
  ),
  lab = list(
    min = c(l = 0, a = -128, b = -128),
    max = c(l = 100, a = 127, b = 127)
  ),
  lch = list(
    min = c(l = 0, c = 0, h = 0),
    max = c(l = 100, c = 140, h = 360)
  ),
  oklab = list(
    min = c(l = 0, a = -0.5, b = -0.5),
    max = c(l = 1, a = 0.5, b = 0.5)
  ),
  oklch = list(
    min = c(l = 0, c = 0, h = 0),
    max = c(l = 1, c = 0.4, h = 360)
  )
)

cn <- colnames(dat) #channel or column names
for(cni in cn){
  dat[,cni] <- (dat[,cni] - cs_ranges[[color_space]]$min[cni]) / 
    (cs_ranges[[color_space]]$max[cni] - cs_ranges[[color_space]]$min[cni])  
}


if(constrain_to_similar_colors){
  
  #for the saturation or chroma channel, if one exists
  if(any(cn %in% c("s", "c"))){
    bounds_sc <- flex_bounds(quantile((dat[,2]), probs = c(0.05,0.95)), rat = fb_rat)  
  } else {
    bounds_sc <- c(0,1)
  }
  #for the luminosity channel
  bounds_l <- flex_bounds(quantile((dat[,"l"]), probs = c(0.1,0.9)), rat = fb_rat)
  
} else {
  bounds_sc <- c(0,1)
  bounds_l <- c(0,1)
}

#### optimize colors ####

#initialize parameter vector
starts <- lapply(1:n_optim_runs, function(i){
  rnorm(n_cols_to_add * ncol(dat), sd = 1)
})

#do optimization
all_output <- do.call(rbind, parallel::mclapply(1:n_optim_runs, function(i){
  opt <- optimx::optimx(
    par       = starts[[i]],
    fn        = mean_dist_func,
    curr_cols        = dat,
    bounds_sc       = bounds_sc,
    bounds_l        = bounds_l,
    color_space     = color_space,
    colorblind_weights = colorblind_weights,
    cs_ranges       = cs_ranges,
    method          = c("BFGS")[1],
    colorblind_safe = colorblind_safe)
  return(opt)
}, mc.cores = 12))
  
#find optimal color
output <- all_output[which.min(all_output$value),]
raw_pars <- unlist(output[paste0("p", 1:length(starts[[1]]))])
new_cols <- mean_dist_func(par = raw_pars,
                           curr_cols        = dat,
                           bounds_sc       = bounds_sc,
                           bounds_l        = bounds_l,
                           color_space     = color_space,
                           colorblind_weights = colorblind_weights,
                           cs_ranges       = cs_ranges,
                           colorblind_safe = colorblind_safe, 
                           return_info = T)$new_hex


#### plot results ####

png(filename = "~/color_addition.png", width = 750 * length(plot_cbs), height = 1000, pointsize = 30)
par(mfcol = c(2,4), mar = c(0,0,0,0))

for(plot_cb in plot_cbs){

#plot colorspace w/ old and new colors
plot(x = NULL, y = NULL, pch = 15, col = cols, cex = 10, xlim = c(0,5), ylim = c(0, 1),
     frame.plot = F, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
max_col_num <- max(length(cols), length(new_cols))
col_height <- 0.9 / max_col_num
rect(xleft = rep(0.5, length(cols)), xright = rep(1.5, length(cols)), 
     ybottom = cumsum(c(0, rep(col_height, length(cols)-1))),
     ytop = cumsum(rep(col_height, length(cols))), col = cvd_convert(cols = cols, type = plot_cb))
rect(xleft = rep(2, length(new_cols)), xright = rep(3, length(new_cols)), 
     ybottom = cumsum(c(0, rep(col_height, length(new_cols)-1))),
     ytop = cumsum(rep(col_height, length(new_cols))), col = cvd_convert(cols = new_cols, type = plot_cb))

#should we plot color hex in black or white?
curr_tcol <- c(0:1)[1 + (farver::decode_colour(cvd_convert(cols = cols, type = plot_cb), to = "hsl")[,"l"] > 50)]
new_tcol <- c(0:1)[1 + (farver::decode_colour(cvd_convert(cols = new_cols, type = plot_cb), to = "hsl")[,"l"] > 50)]

#color hex values
text(x = 1, y = cumsum(c(0, rep(col_height, length(cols)-1))) + col_height / 2, 
     labels = cols, col = curr_tcol)
text(x = 2.5, y = cumsum(c(0, rep(col_height, length(new_cols)-1))) + col_height / 2, 
     labels = new_cols, col = new_tcol)

#color labels
text(1, length(cols) * col_height, cex = 1.5, "current colors = O", srt = 0, pos = 3)
text(2.5, length(new_cols) * col_height, cex = 1.5, "new colors = □", srt = 0, pos = 3)

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
plot(1,1,col = "white", xlim = c(-1,1), ylim = c(-1,1), frame.plot = F, xlab = "", 
     ylab = "", xaxt = "n", yaxt = "n", asp = 1)
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
    
    # Get center hue/chroma of this slice
    hue_deg   <- (degrees[ts] + dd) / (2 * pi) * 360
    chroma    <- (radii[rs] + rd)
    
    # Build input for encode_colour depending on space
    fixed_lightness <- ifelse(colorwheel_space == "hsl", 0.5, 0.75)
    wheel_col <- switch(colorwheel_space,
                        hsl   = encode_colour(t(c(h = hue_deg,
                                                  s = chroma * cs_ranges[["hsl"]]$max["s"],
                                                  l = fixed_lightness * cs_ranges[["hsl"]]$max["l"])), from = "hsl"),
                        lch   = encode_colour(t(c(l = fixed_lightness * cs_ranges[["lch"]]$max["l"],
                                                  c = chroma * cs_ranges[["lch"]]$max["c"],
                                                  h = hue_deg)), from = "lch"),
                        oklch = encode_colour(t(c(l = fixed_lightness,
                                                  c = chroma,
                                                  h = hue_deg)), from = "oklch"),
                        stop("Unsupported colorwheel_space: ", colorwheel_space)
    )
    wheel_col <- cvd_convert(cols = wheel_col, type = plot_cb)
    polygon(x = coords[,1], y = coords[,2], border = NA, col = wheel_col)
  }
}

# Decode to chosen colorspace
colorwheel_cols     <- farver::decode_colour(cols, to = colorwheel_space)
colorwheel_new_cols <- farver::decode_colour(new_cols, to = colorwheel_space)

# Get coordinate mapping for current and new colors
s_or_c <- intersect(c("s", "c"), names(cs_ranges[[colorwheel_space]]$max))
coords_curr_cols <- t(apply(colorwheel_cols, 1, function(x)
  polar2cart(t = x["h"] / 360 * 2 * pi,
             r = x[s_or_c] / cs_ranges[[colorwheel_space]]$max[s_or_c])
))

coords_new_cols <- t(apply(colorwheel_new_cols, 1, function(x)
  polar2cart(t = x["h"] / 360 * 2 * pi,
             r = x[s_or_c] / cs_ranges[[colorwheel_space]]$max[s_or_c])
))

# Choose point sizes based on lightness (normalize to 0‥1 if needed)
l_range <- cs_ranges[[colorwheel_space]]               # same list you use elsewhere
l_min   <- l_range$min["l"];  l_max <- l_range$max["l"]

l_old_norm <- (colorwheel_cols[,"l"]     - l_min) / (l_max - l_min)
l_new_norm <- (colorwheel_new_cols[,"l"] - l_min) / (l_max - l_min)

cex_old <- 0.5 + 2 * l_old_norm
cex_new <- 0.5 + 2 * l_new_norm

# Plot original and new colors
points(coords_curr_cols[,1], coords_curr_cols[,2],
       pch = 21, bg = cvd_convert(cols = cols, type = plot_cb),     col = "black", cex = cex_old)
points(coords_new_cols[,1],  coords_new_cols[,2],
       pch = 22, bg = cvd_convert(cols = new_cols, type = plot_cb), col = "black", cex = cex_new)

legend("topleft", legend = paste0("colorwheel in\n", colorwheel_space, " colorspace",
                                  ifelse(plot_cb == "none", "", paste0("\n(", plot_cb, "-type)"))), bty = "n", cex = 1.2)

}

dev.off()