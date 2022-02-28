#convenience functions
addImg <- function( #taken from https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r/56280430
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

#horizontal values (avoiding 0 & 1 as dbeta can overflow for some values, and Inf densities mess with polygon())
xr <- seq(-0.01, 1.01, length.out = 1000)

#plotting params
boat <- png::readPNG("boat.png")
col <- "#4682B4"
n_waves <- 3
n_transitions <- 5
fps <- 60 #frames per second
fpt <- 60 #frames per transition
frame_indices <- 1:((fpt*(n_transitions-1))-1)

#pre-compute Beta densities (ie heights of waves)
dens <- lapply(1:n_transitions, function(ti){
  mean_beta <- runif(n_waves,0.1,0.9)
  conc_beta <- 5 + rexp(n_waves, 0.1)
  a <- mean_beta * conc_beta
  b <- (1-mean_beta) * conc_beta
  dens <- t(sapply(1:n_waves, function(i) dbeta(xr, a[i], b[i]) + (sin(seq(0,2*pi*10+i,length.out = length(xr))) + 1) / (10 * i)))
  dens
})

#create directories
if(!dir.exists("temp")){dir.create("temp")}
if(!dir.exists("temp/frames")){dir.create("temp/frames")}
if(length(list.files("temp/frames/")) > 0){file.remove(paste0("temp/frames/", list.files("temp/frames/")))}

#iterate through and draw all the frames
for(fri in frame_indices){
  
  cat(paste0(fri, " "))
  
  #get current transition parameters
  prop_t <- fri %% fpt / fpt
  curr_t <- ceiling(fri / fpt) + ifelse(fri %% fpt == 0, 1, 0)
  
  #get our interpolated matrix of densities (equiv to a mixture of Betas)
  den <- dens[[curr_t]] * (1-prop_t) + dens[[curr_t+1]] * (prop_t)
  
  #start plotting
  png(filename = paste0("temp/frames/", paste0(rep(0, 5 - nchar(fri)), collapse = ""), fri, ".png"), 
      width = 1000, height = 500)
  
  plot(NULL, xaxt = "n", yaxt = "n", ylim = range(dens), xlim = range(xr), xlab = "", ylab = "",
       main = "Sailing the Probabili-seas", cex.main = 3, col.main = col)
  
  #make the sky blue
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = adjustcolor("lightblue", 0.4))
  
  #plot a sun in the top left corner
  points(x = par("usr")[1] - 0.01, y = par("usr")[4] * 0.97, cex = 20, pch = 19, col = "#FDD400")
  
  #add in our boat at the top of the tallest wave at its current location
  addImg(boat, x = fri / max(frame_indices), y = max(den[,which.min(abs(fri / max(frame_indices) - xr))]) + 0.5, width = 0.1)
  
  #color in the waves
  for(i in 1:n_waves){
    polygon(x = c(xr, rev(xr)), y = c(den[i,], rep(0, length(den[i,]))), col = adjustcolor(col, 0.2), border = col, lwd = 2)  
  }
  
  dev.off()  
  
}

#stich together animation
file.remove(paste0("~/temp/", list.files("~/temp", pattern = "*.mp4")))
system(paste0("cd \'temp\'; ffmpeg -r ", 
              fps," -f image2 -s 1000x500 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p temp.mp4"))

#reverse and append and loop
system(paste0("cd \'temp\'; ", 
              "ffmpeg -i temp.mp4 -vf reverse rev_temp.mp4; ",
              "touch input.txt;",
              "echo \"file temp.mp4\nfile rev_temp.mp4\" > input.txt;",
              "ffmpeg -f concat -i input.txt -codec copy 1t_temp.mp4; ",
              "ffmpeg -stream_loop 1 -i 1t_temp.mp4 -c copy final.mp4"))
