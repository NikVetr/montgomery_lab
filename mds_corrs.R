#### load packages ####
library(ape)
library(phangorn)
library(smacof)
library(MotrpacBicQC)
library(pheatmap)
library(data.table)
library(RColorBrewer)

#### specify some custom functions ####
polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
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

#### load Nicole's correlation matrix ####
load("~/data/smontgom/nicole_zcor.RData")

#### get colors ####
col_df = data.frame(row.names=colnames(zcor))
col_df$Tissue = gsub('_.*','',rownames(col_df))
col_df$Time = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[3])
col_df$Sex = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[2])

tissue_cols = bic_animal_tissue_code[,tissue_hex_colour]
names(tissue_cols) = bic_animal_tissue_code[,tissue_name_release]

colours = list(Tissue=tissue_cols[unique(col_df$Tissue)], 
               Time=group_cols[unique(col_df$Time)],
               Sex=sex_cols[c('male','female')])

#fix colors
colours_2 <- colours; 
colours_2$Tissue <- colours_2$Tissue[-which(names(colours_2$Tissue) == "t63-testes")]
names(colours_2$Tissue) <- gsub(names(colours_2$Tissue), pattern = "ovaries", replacement = "gonad")
tissues_2 <- sapply(1:length(names(colours_2$Tissue)), function(x) paste0(strsplit(names(colours_2$Tissue)[x], "-")[[1]][-1], collapse = "-"))
colours_2$Sex <- rev(colours_2$Sex)

#### specify plotting parameters ####
color_by_absolute_correlations <- F

lawOfCosines <- F
spearmanDist <- F
FVU <- T
absoluteValue <- F
footrule <- F
oneMinus <- F

if(lawOfCosines){
  distmat <- sqrt(2-2*(zcor))
  extra = "lawOfCosines"
} else if(spearmanDist){
  distmat <- as.matrix(factoextra::get_dist(zcor, method = "spearman"))
  extra = "spearmanDist"
} else if(FVU){
  distmat <- 1 - zcor^2
  extra = "FVU"
  fossil::tri.ineq(distmat)
} else if(absoluteValue){
  distmat <- sqrt(2-2*abs(zcor))
  extra = "lawOfCosines_||"
} else if(oneMinus){
  distmat <- 1-zcor
  extra = "oneMinus"
} else if(footrule){
  d <- as.data.frame(wide_complete)
  d <- apply(d, 2, order)
  d <- apply(d, 2, rev) #now each col is an ordered list
  distmat <- sapply(1:ncol(d), function(col1) sapply(1:ncol(d), function(col2) sum(abs(d[,col1] - d[,col2]))))
  rownames(distmat) <- colnames(distmat) <- colnames(d)
  distmat <- distmat / max(distmat)
  extra = "footrule"
}
MDSpts <- cmdscale(distmat)
stdMDSpts <- mds(distmat, init = MDSpts, itmax = 1E6, eps = 1E-10)$conf
PCoApts <- pcoa(distmat)$vectors[,1:2] #same as MDSpts

if(!lawOfCosines){
  orig_PCoApts <- PCoApts
  PCoApts[,1] <- PCoApts[,1] - min(PCoApts[,1])
  PCoApts[,1] <- PCoApts[,1] / max(PCoApts[,1]) * 0.8 - 0.6
  PCoApts[,2] <- PCoApts[,2] - min(PCoApts[,2])
  PCoApts[,2] <- PCoApts[,2] / max(PCoApts[,2]) * 1.1 - 0.5
  
  orig_stdMDSpts <- stdMDSpts
  stdMDSpts[,1] <- stdMDSpts[,1] - min(stdMDSpts[,1])
  stdMDSpts[,1] <- stdMDSpts[,1] / max(stdMDSpts[,1]) * 1.525 - 0.8
  stdMDSpts[,2] <- stdMDSpts[,2] - min(stdMDSpts[,2])
  stdMDSpts[,2] <- stdMDSpts[,2] / max(stdMDSpts[,2]) * 1.5 - 0.75
  if(absoluteValue){
    stdMDSpts[,1] <- stdMDSpts[,1] * 0.9
  }
  
}

tiss_cols <- unlist(sapply(col_df$Tissue, function(tiss) colours_2$Tissue[which(names(colours_2$Tissue) == tiss)]))
tiss_names <- stringr::str_to_title(sapply(names(colours_2$Tissue), function(x) paste0(strsplit(x, "-")[[1]][-1], collapse = " ")))
sex_cols <- sapply(col_df$Sex, function(sx) colours_2$Sex[which(names(colours_2$Sex) == sx)])
corr_thresh <- 0.3
if(FVU | absoluteValue){
  edges_to_draw <- which(abs(zcor) > corr_thresh, arr.ind = T)
  edges_to_draw <- edges_to_draw[apply(edges_to_draw, 1, diff) != 0,]
  if(color_by_absolute_correlations){
    corrs_of_edges <- abs(zcor[edges_to_draw])
  } else {
    corrs_of_edges <- (zcor[edges_to_draw])
  }
} else {
  edges_to_draw <- which((zcor) > corr_thresh, arr.ind = T)
  edges_to_draw <- edges_to_draw[apply(edges_to_draw, 1, diff) != 0,]
  corrs_of_edges <- zcor[edges_to_draw]
}
colgrad <- viridis::viridis(n = 100)

if(color_by_absolute_correlations | !(FVU | absoluteValue)){
  corr_edge_col <- colgrad[ceiling((corrs_of_edges - corr_thresh) * 100 / (1-corr_thresh))]
  corr_edge_weight <- corrs_of_edges * 10 - (corr_thresh*10-1)
} else {
  corr_edge_col <- colgrad[ceiling((corrs_of_edges + 1) * 50)]
  corr_edge_weight <- abs(corrs_of_edges) * 10 - (corr_thresh*10-1)
}


#### figure 1 ####
grDevices::cairo_pdf(filename = paste0("~/Documents/figure2a_suggestion/PCoA_Viz_", extra,".pdf"), width = 1151 / 72 / 1.3, height = 931 / 72 / 1.3)

par(mar = c(0,0,2,5), xpd = T)
plot(1,1,xlim = c(-0.6,0.7), ylim = c(-0.5, 0.6), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
segments(x0 = PCoApts[edges_to_draw[,1],1], y0 = PCoApts[edges_to_draw[,1],2], 
         x1 = PCoApts[edges_to_draw[,2],1], y1 = PCoApts[edges_to_draw[,2],2], 
         col = corr_edge_col, lwd = corr_edge_weight)
points(PCoApts, pch = 19, cex = 5, col = sex_cols)
points(PCoApts, pch = 19, cex = 4, col = tiss_cols)
text(x = PCoApts[,1], y = PCoApts[,2], labels = col_df$Time, col = "white", font = 2)

plotMatrix(t(t(rev(tiss_names))), size = c(0.15,0.45), 
           location = c(0.35,0.175), title = F, rownames = F, colnames = F, cex = 0.9)
points(rep(0.52, 18), 0.16 + 1:18 / 39.5, col = colours_2$Tissue, pch = 19, cex = 2)

plotMatrix(t(t(rev(c("Female", "Male")))), size = c(0.075,0.05), 
           location = c(0.55,0.175), title = F, rownames = F, colnames = F, cex = 0.9)
points(rep(0.64275, 2), 0.16 + 1:2 / 39.5, col = colours_2$Sex, pch = 19, cex = 2)
points(rep(0.64275, 2), 0.16 + 1:2 / 39.5, col = "white", pch = 19, cex = 1.5)

xl <- 0.575; yb <- 0.245; xr <- 0.625; yt <- 0.625;
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=colgrad[1:20*5-4]
)
lws <- (c(0.5, 1) * 10 - 4) / 72 / par("pin")[1]
hoffset <- 0.015

if(color_by_absolute_correlations | !(FVU | absoluteValue)){
  text(labels = round(seq(corr_thresh, 1, length.out = 10), 2), y = seq(yb, yt, length.out = 10), x = xr - 0.0075, pos = 4, las=2, cex=0.9)
  text(labels = latex2exp::TeX(paste0(ifelse(FVU | absoluteValue, "|", ""), "$\\rho$", ifelse(FVU | absoluteValue, "|", ""))), y = yt - 0.01, x = (xl + xr) / 2, pos = 3, cex = 2)
  polygon(x = c(xl - hoffset - lws[2]/2,
                xl - hoffset - lws[1]/2,
                xl - hoffset + lws[1]/2,
                xl - hoffset + lws[2]/2), 
          y = c(yt, yb,yb,yt), col = "black")
} else {
  text(labels = round(seq(-1, 1, length.out = 10), 2), y = seq(yb, yt, length.out = 10), x = xr - 0.0075, pos = 4, las=2, cex=0.9)
  text(labels = latex2exp::TeX(paste0("$\\rho$")), y = yt - 0.01, x = (xl + xr) / 2, pos = 3, cex = 2)
  polygon(x = c(xl - hoffset - lws[1]/2, 
                xl - hoffset - lws[2]/2,
                xl - hoffset + lws[2]/2, 
                xl - hoffset + lws[1]/2),
          y = c((yt + yb) / 2, yb,yb,(yt + yb) / 2), col = "black")
  polygon(x = c(xl - hoffset - lws[2]/2,
                xl - hoffset - lws[1]/2,
                xl - hoffset + lws[1]/2,
                xl - hoffset + lws[2]/2), 
          y = c(yt, (yt + yb) / 2,(yt + yb) / 2,yt), col = "black")
}

# text("First Two PCoA Axes", x = 0.5, y = -0.34, cex = 3, font = 2)
if(lawOfCosines){
  text(labels = latex2exp::TeX("$d_{ij} \\approx \\sqrt{2(1-\\rho_{ij})$"), x = 0.275, y = 0.6125)
} else if(FVU){
  text(labels = latex2exp::TeX("$d_{ij} \\approx 1-\\rho_{ij}^{2}$"), x = 0.275, y = 0.6125)
} else if(absoluteValue){
  text(labels = latex2exp::TeX("$d_{ij} \\approx \\sqrt{2(1-|\\rho_{ij}|)$"), x = 0.275, y = 0.6125)
} else if(oneMinus){
  text(labels = latex2exp::TeX("$1-\\rho_{ij}$"), x = 0.275, y = 0.6125)
}

dev.off()

#### figure 2 ####

grDevices::cairo_pdf(filename = paste0("~/Documents/figure2a_suggestion/SMACOF_Viz_", extra,".pdf"), width = 1151 / 72 / 1.3, height = 931 / 72 / 1.3)

par(mar = c(2,2,2,10), xpd = T)
plot(1,1,xlim = c(-0.75,0.85), ylim = c(-0.75,0.75), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
segments(x0 = stdMDSpts[edges_to_draw[,1],1], y0 = stdMDSpts[edges_to_draw[,1],2], 
         x1 = stdMDSpts[edges_to_draw[,2],1], y1 = stdMDSpts[edges_to_draw[,2],2], 
         col = corr_edge_col, lwd = corr_edge_weight)
points(stdMDSpts, pch = 19, cex = 5, col = sex_cols)
points(stdMDSpts, pch = 19, cex = 4, col = tiss_cols)
text(x = stdMDSpts[,1], y = stdMDSpts[,2], labels = col_df$Time, col = "white", font = 2)

plotMatrix(t(t(rev(tiss_names))), size = c(0.2,0.6), 
           location = c(0.85,0.2), title = F, rownames = F, colnames = F, cex = 0.9)
points(rep(1.075, 18), 0.1775 + 1:18 / 29.50, col = colours_2$Tissue, pch = 19, cex = 2)

plotMatrix(t(t(rev(c("Female", "Male")))), size = c(0.1,0.065),
           location = c(0.69, 0.8 - 0.6 / 9), title = F, rownames = F, colnames = F, cex = 0.9)
points(rep(0.8135, 2), 0.7125 + 1:2 / 29.5, col = colours_2$Sex, pch = 19, cex = 2)
points(rep(0.8135, 2), 0.7125 + 1:2 / 29.5, col = "white", pch = 19, cex = 1.5)

xl <- 1; yb <- -0.3; xr <- 1.05; yt <- 0.125;
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=colgrad[1:20*5-4]
)
lws <- (c(0.5, 1) * 10 - 4) / 72 / par("pin")[1]
hoffset <- 0.02



if(color_by_absolute_correlations | !(FVU | absoluteValue)){
  text(labels = round(seq(corr_thresh, 1, length.out = 10), 2), y = seq(yb, yt, length.out = 10), x = xr - 0.0075, pos = 4, las=2, cex=0.9)
  text(labels = latex2exp::TeX(paste0(ifelse(FVU | absoluteValue, "|", ""), "$\\rho$", ifelse(FVU | absoluteValue, "|", ""))), y = yt - 0.01, x = (xl + xr) / 2, pos = 3, cex = 2)
  polygon(x = c(xl - hoffset - lws[2]/2,
                xl - hoffset - lws[1]/2,
                xl - hoffset + lws[1]/2,
                xl - hoffset + lws[2]/2), 
          y = c(yt, yb,yb,yt), col = "black")
} else {
  text(labels = round(seq(-1, 1, length.out = 10), 2), y = seq(yb, yt, length.out = 10), x = xr - 0.0075, pos = 4, las=2, cex=0.9)
  text(labels = latex2exp::TeX(paste0("$\\rho$")), y = yt - 0.01, x = (xl + xr) / 2, pos = 3, cex = 2)
  polygon(x = c(xl - hoffset - lws[1]/2, 
                xl - hoffset - lws[2]/2,
                xl - hoffset + lws[2]/2, 
                xl - hoffset + lws[1]/2),
          y = c((yt + yb) / 2, yb,yb,(yt + yb) / 2), col = "black")
  polygon(x = c(xl - hoffset - lws[2]/2,
                xl - hoffset - lws[1]/2,
                xl - hoffset + lws[1]/2,
                xl - hoffset + lws[2]/2), 
          y = c(yt, (yt + yb) / 2,(yt + yb) / 2,yt), col = "black")
}

# text("Normalized\nMDS Space\n(SMACOF)", x = 0.9, y = -0.65, cex = 2.45, font = 2)
if(lawOfCosines){
  text(labels = latex2exp::TeX("$d_{ij} \\approx \\sqrt{2(1-\\rho_{ij})$"), x = 0.75, y = 0.7)
} else if(FVU){
  text(labels = latex2exp::TeX("$d_{ij} \\approx 1-\\rho_{ij}^{2}$"), x = 0.75, y = 0.7)
} else if(absoluteValue){
  text(labels = latex2exp::TeX("$d_{ij} \\approx \\sqrt{2(1-|\\rho_{ij}|)$"), x = 0.75, y = 0.7)
} else if(oneMinus){
  text(labels = latex2exp::TeX("$1-\\rho_{ij}$"), x = 0.75, y = 0.7)
}

dev.off()

#### examine distance preserving function ####
PCoApts <- orig_PCoApts
stdMDSpts <- orig_stdMDSpts

PCoApts_full <- pcoa(distmat)$vectors
distvec <- as.vector(distmat[upper.tri(distmat)])
PCoA_dists <- lapply(1:ncol(PCoApts_full), function(d) as.vector(as.matrix(dist(PCoApts_full[,1:d]))[upper.tri(distmat)]))
distvec_sort <- sort(distvec)
PCoA_dists_ordered <- lapply(1:length(PCoA_dists), function(d) PCoA_dists[[d]][order(distvec)])
ww <- 0.1; intvs <- seq(0,2,ww)
distvec_means <- sapply(2:length(intvs), function(w) mean(distvec_sort[distvec_sort > intvs[w-1] & distvec_sort < intvs[w]]))
PCoA_means <- lapply(1:length(PCoA_dists_ordered), function(do) sapply(2:length(intvs), function(w) mean(PCoA_dists_ordered[[do]][distvec_sort > intvs[w-1] & distvec_sort < intvs[w]])))


ww <- 0.2; slide_incr <- 0.01
intvs <- seq(0,2-ww,slide_incr)
distvec_means <- sapply(1:length(intvs), function(w) mean(distvec_sort[distvec_sort > intvs[w] & distvec_sort < (intvs[w] + ww)]))
stdMDS_means <- sapply(1:length(intvs), function(w) mean(as.matrix(dist(stdMDSpts))[upper.tri(as.matrix(dist(stdMDSpts)))][order(distvec)][distvec_sort > intvs[w] & distvec_sort < (intvs[w] + ww)]))
PCoA_means <- lapply(1:length(PCoA_dists_ordered), function(do) sapply(1:length(intvs), function(w) mean(PCoA_dists_ordered[[do]][distvec_sort > intvs[w] & distvec_sort < (intvs[w] + ww)])))
cols3 <- viridis::viridis(length(PCoA_dists_ordered), alpha = 0.6)

grDevices::cairo_pdf(filename = paste0("~/Documents/figure2a_suggestion/poor_embedding_", extra,".pdf"), width = 400 / 72, height = 400 / 72)
par(mar = c(4,4,2,1))
plot(1,1,xlim = c(0,ceiling(max(distvec_means, na.rm = T) * 10) / 10), 
     ylim = c(0,ceiling(max(unlist(PCoA_means), na.rm = T) * 10) / 10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
box(which = "plot", lwd = 2)
abline(a = 0, b = 1, xpd = F, lwd = 2)
# points(distvec, PCoA_dists[[40]], pch = 19, col = rgb(0,0,0,0.5))
# points(distvec_sort, PCoA_dists_ordered[[40]], pch = 19, col = rgb(0,0,0,0.5))
for(i in 1:length(PCoA_means)){
  # lines(c(distvec_means[!is.na(distvec_means)]), c(PCoA_means[[i]][!is.na(PCoA_means[[i]])]), lwd = 2, col = cols3[i])
  lines(c(distvec_means), c(PCoA_means[[i]]), lwd = 2, col = cols3[i])
}

segments(x0 = min(distvec_means, na.rm = T), y0 = 0, x1 = min(distvec_means, na.rm = T), y1 = min(distvec_means, na.rm = T), lwd = 2, lty = 3, col = "black")
text(x = min(distvec_means, na.rm = T), y = min(distvec_means, na.rm = T) * 2 / 3, labels = latex2exp::TeX("min. d_{ij}"), pos = 2, srt = 90)
lines(c(distvec_means[!is.na(distvec_means)]), c(stdMDS_means[!is.na(stdMDS_means)]), lwd = 3, col = "orange")
axis(1, at = seq(0,2,by=0.2), labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = 0)
mtext(text = seq(0,2,by=0.2), side = 1, at = seq(0,2,by=0.2), cex = 1, line = 0.35)
axis(2, at = seq(0,2,by=0.2), labels =  rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = 0)
mtext(text = seq(0,2,by=0.2), side = 2, at = seq(0,2,by=0.2), cex = 1, line = 0.5, las = 2)
mtext(text = "Euclidean Distance in PCoA Space", side = 2, cex = 1.5, line = 2)
mtext(text = "Distance in Distance Matrix", side = 1, cex = 1.5, line = 2)


xl <- -0.025; 
yb <- ceiling(max(unlist(PCoA_means), na.rm = T) * 10) / 10 * 0.45; 
xr <- ceiling(max(distvec_means, na.rm = T) * 10) / 10 * 0.035; 
yt <- ceiling(max(unlist(PCoA_means), na.rm = T) * 10) / 10 * 0.95;

rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=cols3[ceiling(seq(1, length(cols3), length.out = 20))], 
  border = cols3[ceiling(seq(1, length(cols3), length.out = 20))]
)
rect(xl,yb,xr,yt, border = 1, lwd = 2)
text(labels = ceiling(seq(1, length(cols3), length.out = 5)), y = seq(yb, yt, length.out = 5), x = xr - 0.02, pos = 4, las=2, cex=0.9)
text(labels = latex2exp::TeX("N_d"), y = yt * 0.98, x = (xl + xr) / 2 -.0025, pos = 3, cex = 0.9, font = 2)
title("Distance-Preserving Inadequacy", cex.main = 1.5, line = 0.25)
text("(means in sliding windows of width 0.2)", x = ceiling(max(distvec_means, na.rm = T) * 10) / 20, y = yt * 1.1, pos = 1)
segments(x0 = xl, y0 = yb * 8 / 9, x1 = xr * 2, y1 = yb * 8 / 9, lwd = 3, col = "orange")
text(x = xr * 2, labels = "2D-SMACOF", y = yb * 79 / 90, pos = 4)
text(x = ceiling(max(unlist(PCoA_means), na.rm = T) * 20) / 30, 
     labels = "1-to-1 line", y = ceiling(max(unlist(PCoA_means), na.rm = T) * 20) / 30, pos = 3, srt = 45)

dev.off()