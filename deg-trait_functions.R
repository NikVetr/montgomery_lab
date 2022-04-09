
#### define functions ####
polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
}

polarp <- function (t, r, col = 1, pch = 1, cex = 1, center = c(0,0)) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  points(xy[, 1] + center[1], xy[, 2] + center[2], cex = cex, col = col, pch = pch)
}

polarl <- function (t, r, col = 1, lwd = 1, center = c(0,0), adjx = 1) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  lines(xy[, 1] * adjx + center[1], xy[, 2] + center[2], lwd = lwd, col = col)
}

arc <- function(t1,t2,r1,r2,res=50,lwd = 1,col=1, mindist = T, self_adjust = 2 * pi / 45, random_selfing = T, clockwise_selfing = T, pointy_selfing = F,
                center = c(0,0), adjx = 1){
  if(mindist){
    if(abs(t1-t2) > pi){
      if(t1 > t2){
        t1 <- t1-2*pi
      } else {
        t2 <- t2-2*pi
      }
    }
  }
  if(abs(t1-t2) < 1E-6){
    
    if(random_selfing){
      lor <- sample(c(-1,1), 1)
    } else {
      lor <- ifelse(clockwise_selfing, -1, 1)
    }
    
    if(pointy_selfing){
      ts <- c(seq(t1, t1 + lor * self_adjust, length.out = res/2), seq(t1 + lor * self_adjust, t1, length.out = res/2))
      rs <- seq(r1, r2, length.out = res)
    } else {
      center <- center + pol2cart(cbind(t2, mean(c(r1, r2))))
      if(random_selfing){
        ts <- seq(t1, t1 + sample(c(-pi, pi), 1), length.out = res)
      } else {
        ts <- seq(t1, t1 + ifelse(clockwise_selfing, -pi, pi), length.out = res)
      }
      rs <- seq(max(c(r1,r2)) - mean(c(r1,r2)), mean(c(r1,r2)) - min(c(r1,r2)), length.out = res)
    }
    
  } else {
    ts <- seq(t1, t2, length.out = res)
    rs <- seq(r1, r2, length.out = res)
  }
  polarl(ts, rs, lwd = lwd,col=col, center = center, adjx = adjx)
}

plotMatrix <- function(mobject, size, location, lwd = 2, lwd_inner = 1.5, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])), lwd = lwd_inner)
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])), lwd = lwd_inner)
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}

draw.contour<-function(a,alpha=0.95,plot.dens=FALSE, line.width=2, line.type=1, limits=NULL, density.res=300,spline.smooth=-1,...){
  ##a is a list or matrix of x and y coordinates (e.g., a=list("x"=rnorm(100),"y"=rnorm(100)))
  ## if a is a list or dataframe, the components must be labeled "x" and "y"
  ## if a is a matrix, the first column is assumed to be x, the second y
  ##alpha is the contour level desired
  ##if plot.dens==TRUE, then the joint density of x and y are plotted,
  ##   otherwise the contour is added to the current plot.
  ##density.res controls the resolution of the density plot
  
  ##A key assumption of this function is that very little probability mass lies outside the limits of
  ## the x and y values in "a". This is likely reasonable if the number of observations in a is large.
  
  require(MASS)
  require(ks)
  if(length(line.width)!=length(alpha)){
    line.width <- rep(line.width[1],length(alpha))
  }
  
  if(length(line.type)!=length(alpha)){
    line.type <- rep(line.type[1],length(alpha))
  }
  
  if(is.matrix(a)){
    a=list("x"=a[,1],"y"=a[,2])
  }
  ##generate approximate density values
  if(is.null(limits)){
    limits=c(range(a$x),range(a$y))
  }
  f1<-kde2d(a$x,a$y,n=density.res,lims=limits)
  
  ##plot empirical density
  if(plot.dens) image(f1,...)
  
  if(is.null(dev.list())){
    ##ensure that there is a window in which to draw the contour
    plot(a,type="n",xlim=limits[1:2],ylim=limits[3:4],...)
  }
  
  ##estimate critical contour value
  ## assume that density outside of plot is very small
  
  zdens <- rev(sort(f1$z))
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)])
  for(cont.level in 1:length(alpha)){
    ##This loop allows for multiple contour levels
    crit.val <- zdens[max(which(Czdens<=alpha[cont.level]))]
    
    ##determine coordinates of critical contour
    b.full=contourLines(f1,levels=crit.val)
    for(c in 1:length(b.full)){
      ##This loop is used in case the density is multimodal or if the desired contour
      ##  extends outside the plotting region
      b=list("x"=as.vector(unlist(b.full[[c]][2])),"y"=as.vector(unlist(b.full[[c]][3])))
      
      ##plot desired contour
      line.dat<-xspline(b,shape=spline.smooth,open=TRUE,draw=FALSE)
      lines(line.dat,lty=line.type[cont.level],lwd=line.width[cont.level])
    }
  }
}

minwhich <- function(x){min(which(x))}

maxwhich <- function(x){max(which(x))}

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

overlaps_with_zero <- function(qi){qi[1,] < 0 & qi[2,] > 0}

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

subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}


prop_greater_than_0 <- function(x) mean(x>0)

logit <- function(p) log(p/(1-p))
logit_prime <- function(p) (1/p + 1 / (1-p))
invlogit <- function(x) {exp(x) / (1+exp(x))}

intersect_rec <- function(x) {
  if(length(x) <= 1){
    return(x)
  }else if(length(x) > 2){
    x[[1]] <- intersect(x[[1]], x[[2]])
    intersect_rec(x[-2])
  } else {
    intersect(x[[1]], x[[2]])
  }
}

logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
# squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
# unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
squish_middle_x <- function(x,f) asinh(x*f)
unsquish_middle_x <- function(x,f) sinh(x)/f
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, shrinkX = 0.95, shrinkY = 0.95, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1*shrinkX, y1*shrinkY, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

line <- function(t,r1,r2,lwd = 1,col=1, center = c(0,0)){
  polarl(rep(t, 2), c(r1,r2), lwd = lwd,col=col, center = center)
}
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
# squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
# unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
squish_middle_x <- function(x,f) asinh(x*f)
unsquish_middle_x <- function(x,f) sinh(x)/f
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}

text_indivcolor <- function(labels_mat, colors_mat, xloc, cex, ...){
  text(labels = labels_mat[,1], col = colors_mat[,1], x = xloc, cex = cex, ...)
  for(i in 2:ncol(labels_mat)){
    prior_string_widths <- apply(do.call(cbind, lapply(1:(i-1), function(li) strwidth(labels_mat[,li], units = "user"))), 1, sum) * cex
    text(labels = labels_mat[,i], col = colors_mat[,i], x = xloc + prior_string_widths, cex = cex, ...)
  }
}

ifelse3 <- function(bool, opt1, opt2){sapply(bool, function(bool_i) ifelse(bool_i, opt1, opt2))}
min2 <- function(vec, opt2){sapply(vec, function(veci) min(veci, opt2))}

#functions from transcr vs. prot stuff
text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = cols[char_i], cex = cex, ...)
  }
}

find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 0.95, srt_height_rescaling_ratio = 1.4){
  
  strwidths <- strwidth(txt)
  strheight <- strheight(txt[1]) * srt_height_rescaling_ratio 
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.5)
  }
  # plot(1:120/100, sapply(log(1:120/100), function(cex) compute_space_remaining(data, cex)), type = "l")
  opt_cex <- suppressWarnings(optimx::optimx(par = par, fn = compute_space_remaining, data = data, hess = NULL, 
                                             method = c('Nelder-Mead'), hessian=FALSE, #can't compute hessian bc of sharp jumps when new line is formed? or maybe not?
                                             control = list(maxit = 1E4, trace = 0, kkt=FALSE)))
  # compute_space_remaining(data = data, par = opt_cex$p1)
  return(list(cex = exp(opt_cex$p1), 
              words_on_lines = put_words_on_lines(data = data, par = exp(opt_cex$p1)),
              space_width_min = space_width_min * exp(opt_cex$p1),
              vertical_space = strheight * exp(opt_cex$p1))
  )
  
}

compute_space_remaining <- function(data, par){
  
  #clarify par-cex relationship
  cex <- exp(par)
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  
  #check that no words are wider than a line
  if(any(strwidths_mod > rectwidth_mod)){
    return(Inf)
  }
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi]
  txt_lines <- list()
  txt_lines[[linei]] <- txi
  
  while(txi < length(txt)){
    txi <- txi + 1
    txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
    current_width <- current_width + strwidths_mod[txi] + space_width_min_mod
    if(current_width > rectwidth_mod){
      txt_lines[[linei]] <- txt_lines[[linei]][-length(txt_lines[[linei]])]
      linei <- linei + 1
      txt_lines[[linei]] <- txi
      current_width <- strwidths_mod[txi]
    }
  }
  
  last_line_width_remaining <- rectwidth_mod - sum(strwidths_mod[txt_lines[[linei]]])
  current_height <- linei * strheight_mod
  space_remaining <- (rectheight_mod - current_height) * rectwidth_mod + (last_line_width_remaining * strheight_mod)
  
  if(space_remaining < 0){return(Inf)} else {return(space_remaining)}
  
}

put_words_on_lines <- function(data, par){
  
  #clarify par-cex relationship
  cex <- par
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi]
  txt_lines <- list()
  txt_lines[[linei]] <- txi
  
  while(txi < length(txt)){
    txi <- txi + 1
    txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
    current_width <- current_width + strwidths_mod[txi] + space_width_min_mod
    if(current_width > rectwidth_mod){
      txt_lines[[linei]] <- txt_lines[[linei]][-length(txt_lines[[linei]])]
      linei <- linei + 1
      txt_lines[[linei]] <- txi
      current_width <- strwidths_mod[txi]
    }
  }
  
  return(txt_lines)
  
}

text_wrapped_words <- function(txt, rect_coords, optimal_word_placement_inf, justified = F, str_height_lower_start_ratio = 0.75, 
                               str_width_lefter_start_ratio = 0.01, rect_rescaling_ratio = 0.95, col = "black", multicolor_words = F, cols_list, ...){
  
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(txt) * optimal_word_placement_inf$cex
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
                  string = word_to_write, pos = 4, cols = cols_list[[words_written]])
      } else {
        text(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
             labels = word_to_write, pos = 4, col = col)
      }
      if(justified){
        curr_x <- curr_x + strwidth(word_to_write) * optimal_word_placement_inf$cex + justified_space_between_words[linei]
      } else {
        curr_x <- curr_x + strwidth(word_to_write) * optimal_word_placement_inf$cex + optimal_word_placement_inf$space_width_min  
      }
      
    }
    curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
    curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio - optimal_word_placement_inf$vertical_space * linei
  }
  
}


seq_incl_0 <- function(from, to, length.out){
  incr <- (to - from) / length.out
  testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
  foo <- 1
  while(length(testseq) > length.out){
    incr <- (to - from) / (length.out-foo)
    testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
    foo <- foo + 1
  }
  foo <- -1
  while(length(testseq) < length.out){
    incr <- (to - from) / (length.out-foo)
    testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
    foo <- foo - 1
  }
  return(testseq)
}

dexp_smooth <- function(x, y, r, reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, interpolate_at = NA){
  
  nx <- length(x)
  if(is.na(interpolate_at[1])){
    n <- nx
    w <- t(sapply(1:(n-1), function(i) c(dexp((x[i] - x[0:(i-1)]), rate = r), dexp(0, rate = r), dexp(-(x[i] - x[(i+1):n]), rate = r))))
    w <- rbind(w, dexp(x[n] - x, rate = r))
  } else {
    n <- length(interpolate_at)
    w <- t(sapply(1:(n-1), function(i) c(dexp((interpolate_at[i] - x[x<=interpolate_at[i]]), rate = r), 
                                         dexp((x[x>interpolate_at[i]] - interpolate_at[i]), rate = r))))
    w <- rbind(w, dexp(x[nx] - x, rate = r))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_trunc_tail){
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:(n-1), function(i) c(pexp((x[i] - x[1]), rate = r), pexp(-(x[i] - x[n]), rate = r))))
    } else {
      tw <- t(sapply(1:(n-1), function(i) c(pexp((interpolate_at[i] - x[1]), rate = r), pexp(-(interpolate_at[i] - x[nx]), rate = r))))
    }
    tw[1,] <- tw[2,]
    tw <- rbind(tw, tw[n-1,])
    tw <- 1/tw
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
    
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_n_pts){
    if(is.na(interpolate_at[1])){
      tw <- cbind(0:(n-1), (n-1):0)
      tw <- 1/tw
    } else {
      tw <- sapply(1:n, function(i) sum(interpolate_at[i] >= x))
      tw <- cbind(tw, nx - tw)
    }
    mintw <- apply(tw, 1, min)
    tw[,1] <- tw[,1] / mintw 
    tw[,2] <- tw[,2] / mintw 
    
    tw[1,] <- tw[2,]; tw[n,] <- tw[n-1,]
    
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(fix_endpoints){
    w[1,] <- c(1, rep(0,nx-1))
    w[n,] <- c(rep(0,nx-1), 1)
  }
  
  wy <- c(w %*% t(t(y)))
  wy <- wy / apply(w,1,sum)
  return(wy)
}

NULL


