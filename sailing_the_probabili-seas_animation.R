#packages
library(foreach)
library(doParallel)
library(parallel)

#convenience functions
addImg <- function( #taken from https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r/56280430
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE, # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing
  angle = 0
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
              interpolate = interpolate, angle = angle)
}

#horizontal values (avoiding 0 & 1 as dbeta can overflow for some values, and Inf densities mess with polygon())
xr <- seq(-0.01, 1.01, length.out = 500)

#plotting params
boat <- png::readPNG("boat.png")
col <- "#4682B4"
n_waves <- 3
n_transitions <- 5
fps <- 60 #frames per second
fpt <- 60 #frames per transition
thin <- 1
frame_indices <- 1:((fpt*(n_transitions-1))-1)
plotting_indices <- seq(1, max(frame_indices), by = thin)
skycols_gradient <- viridisLite::plasma(15)

#pre-compute Beta densities (ie heights of waves)
dens <- lapply(1:n_transitions, function(ti){
  mean_beta <- runif(n_waves,0.15,0.85)
  conc_beta <- 5 + rexp(n_waves, 0.1)
  a <- mean_beta * conc_beta
  b <- (1-mean_beta) * conc_beta
  dens <- t(sapply(1:n_waves, function(i) dbeta(xr, a[i], b[i])))
  dens
})

#create directories
if(!dir.exists("temp")){dir.create("temp")}
if(!dir.exists("temp/frames")){dir.create("temp/frames")}
if(length(list.files("temp/frames/")) > 0){file.remove(paste0("temp/frames/", list.files("temp/frames/")))}

#initiate cluster
if(!exists("cl")){
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

#iterate through and draw all the frames
foreach(pi=1:length(plotting_indices), .packages = c("png")) %dopar% {
# for(pi in 1:length(plotting_indices)){
  fri <- plotting_indices[pi]
  
  cat(paste0(fri, " "))
  
  #get current transition parameters
  prop_t <- fri %% fpt / fpt
  curr_t <- ceiling(fri / fpt) + ifelse(fri %% fpt == 0, 1, 0)
  
  #get our interpolated matrix of densities (equiv to a mixture of Betas)
  den <- dens[[curr_t]] * (1-prop_t) + dens[[curr_t+1]] * (prop_t)
  
  #start plotting
  png(filename = paste0("temp/frames/", paste0(rep(0, 5 - nchar(pi)), collapse = ""), pi, ".png"), 
      width = 1000, height = 500)
  
  plot(NULL, xaxt = "n", yaxt = "n", ylim = range(dens), xlim = range(xr), xlab = "", ylab = "",
       main = "Sailing the Probabili-seas", cex.main = 3, col.main = col)
  
  #color in the sky and sun
  xyrat <- diff(range(par("usr")[1:2])) / par("pin")[1] / diff(range(par("usr")[3:4])) * par("pin")[2]
  maxrad <- diff(range(par("usr")[1:2])) / 2
  maxrad <- sqrt(maxrad^2 + (maxrad * par("pin")[2] / par("pin")[1])^2)
  rseq <- seq(maxrad, maxrad/10, length.out = length(skycols_gradient))
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = skycols_gradient[1], border = skycols_gradient[1])
  for(sk_i in 1:length(skycols_gradient)){
    pi <- 3.14159265358979 #bultin constant is bugging out?
    polygon(x = cos(seq(0, 2*pi, length.out = 100)) * rseq[sk_i] + 0.5, 
            y = sin(seq(0, 2*pi, length.out = 100)) * rseq[sk_i] / xyrat + mean(range(dens)) * (max(frame_indices) - fri) / max(frame_indices), 
            col = skycols_gradient[sk_i], border = skycols_gradient[sk_i])
  }
  
  #add in our boat at the top of the tallest (primary) wave at its current location
  addImg(boat, x = fri / max(frame_indices), y = max(den[,which.min(abs(fri / max(frame_indices) - xr))]) + 0.5, width = 0.1)
  
  #color in the waves
  wavelets <- lapply(1:n_waves, function(i) sin(seq(0,2*pi*10 * (i/2)^0.25, length.out = length(xr)) + i * 3 * pi +  fri / max(frame_indices) * 8 * pi) / 10)
  
  #first block off the waves with a uniform color 
  for(i in 1:n_waves){
    polygon(x = c(par("usr")[1], xr, par("usr")[2], par("usr")[2], rev(xr), par("usr")[1]), 
            y = c(0, den[i,] + wavelets[[i]], 0, rep(par("usr")[3], length(den[i,]) + 2)), 
            col = skycols_gradient[round(length(skycols_gradient) / 4)], 
            border = skycols_gradient[round(length(skycols_gradient) / 4)], lwd = 2)  
  }
  #then add in the color itself
  for(i in 1:n_waves){
    polygon(x = c(par("usr")[1], xr, par("usr")[2], par("usr")[2], rev(xr), par("usr")[1]), 
            y = c(0, den[i,] + wavelets[[i]], 0, rep(par("usr")[3], length(den[i,]) + 2)), 
            col = adjustcolor(col, 0.2), border = adjustcolor(col, 0.75), lwd = 2)  
  }
  
  dev.off()  
  
}

#stich together animation
file.remove(paste0("~/temp/", list.files("~/temp", pattern = "*.mp4")))
system(paste0("cd \'temp\'; ffmpeg -r ", 
              fps / thin," -f image2 -s 1000x500 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p temp.mp4"))

#reverse and append and loop
system(paste0("cd \'temp\'; ", 
              "ffmpeg -i temp.mp4 -vf reverse rev_temp.mp4; ",
              "touch input.txt;",
              "echo \"file temp.mp4\nfile rev_temp.mp4\" > input.txt;",
              "ffmpeg -f concat -i input.txt -codec copy 1t_temp.mp4; ",
              "ffmpeg -stream_loop 1 -i 1t_temp.mp4 -c copy final.mp4"))
