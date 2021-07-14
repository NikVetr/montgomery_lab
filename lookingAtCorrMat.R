library(pheatmap)
library(data.table)
library(RColorBrewer)
library(MotrpacBicQC) # loading this package gives you access to several color palettes
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
if(!exists("time_de")){load('~/data/smontgom/time_de.RData')} # this is a data.table
# convert to z-scores
time2 = time_de[,.(feature_ID, sex, comparison_group, tissue, logFC, p_value)]
time2[,z := qnorm(p_value/2)*sign(logFC)]
head(time2)
# cast
wide = dcast(time2, value.var='z', feature_ID ~ tissue + sex + comparison_group)
wide_complete = wide[complete.cases(wide),]
rownames(wide_complete) = wide_complete$feature_ID
wide_complete$feature_ID = NULL
head(wide_complete)

zcor = cor(wide_complete,method='spearman')
col_df = data.frame(row.names=colnames(zcor))
col_df$Tissue = gsub('_.*','',rownames(col_df))
col_df$Time = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[3])
col_df$Sex = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[2])

tissue_cols = bic_animal_tissue_code[,tissue_hex_colour]
names(tissue_cols) = bic_animal_tissue_code[,tissue_name_release]


colours = list(Tissue=tissue_cols[unique(col_df$Tissue)], 
               Time=group_cols[unique(col_df$Time)],
               Sex=sex_cols[c('male','female')])

cols = colorRampPalette(c("blue", "white", "red"))(100)
if(!exists("out")){
  out=pheatmap(zcor,
             scale='none',
             clustering_distance_cols = as.dist(1-zcor), ##aren't ya getting squared distances here? or isn't the 'usual' method of transforming between these two more like sqrt(2-2*cor)
             clustering_distance_rows = as.dist(1-zcor),
             cluster_cols = T, # TRUE or FALSE; by default, clustering adds a dendrogram
             cluster_rows = T, 
             annotation_col = col_df, # if you want to annotate your samples/columns with color bars
             annotation_row = col_df, # if you want to annotate your samples/columns with color bars
             show_rownames = F,
             show_colnames = F,
             annotation_colors = colours, # custom colors
             annotation_names_col = F,
             annotation_names_row = F,
             border_color = NA, # remove border
             treeheight_row = 0,
             #cutree_cols = 6,
             #cutree_rows = 6,
             color = cols,
             fontsize = 7)
}
# out

hist(zcor[upper.tri(zcor)], breaks = 20)
pairwiseCorrs <- zcor[upper.tri(zcor)]
binterval <- 0.1
n_in_bins <- sapply(seq(-1, 1 - binterval, by = binterval), function(x) sum(pairwiseCorrs > x & pairwiseCorrs < (x+binterval)))
prop_in_bins <- n_in_bins / sum(n_in_bins)

# col_df$Tissue <- gsub(col_df$Tissue, pattern = c("testes", "ovaries")[1], replacement = "gonad")
# col_df$Tissue <- gsub(col_df$Tissue, pattern = c("testes", "ovaries")[2], replacement = "gonad")
# col_df$Tissue <- gsub(col_df$Tissue, pattern = "t63", replacement = "t64")
tss <- unique(col_df$Tissue)
tps <- unique(col_df$Time)
sex <- unique(col_df$Sex)
match_conditions <- cbind(rep(c(0,1), 4), rep(c(0,0,1,1), 2), c(rep(0,4), rep(1,4))) == 1; colnames(match_conditions) = c("Tissue", "Timepoint", "Sex")

#yes I could do 1:(nrow(zcor)-1), function(row) sapply((row+1):ncol(zcor) of faster but this is cheap & easy so redundancy is ok
if(!exists("condition_matrix")){
  condition_matrix <- sapply(1:(nrow(zcor)), function(row) sapply(1:ncol(zcor), function(col) 
    which(sapply(1:nrow(match_conditions), function(cond) all((as.logical(col_df[row,] == col_df[col,]) == match_conditions[cond,]))))
  ))
  diag(condition_matrix)
}

pairwise_corrs_per_condition <- lapply(1:nrow(match_conditions), function(condition) (zcor[condition_matrix == condition]))
n_in_bins_per_condition <- t(sapply(1:length(pairwise_corrs_per_condition), function(condition) 
  sapply(seq(-1, 1 - binterval, by = binterval), function(x) sum(pairwise_corrs_per_condition[[condition]] > x & pairwise_corrs_per_condition[[condition]] < (x+binterval))))) / 2
apply(n_in_bins_per_condition, 2, sum) - n_in_bins #sanity check
prop_in_bins_per_condition <- n_in_bins_per_condition %*% diag(1/n_in_bins)
prop_in_bins_per_condition[is.na(prop_in_bins_per_condition)] <- 0

cols <- brewer.pal(8, "Dark2")
intervals <- seq(-1, 1 - binterval, by = binterval)
par(xpd=TRUE)
plot(1,1,xlim = c(-1,1.5), ylim = c(0,1.25), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
curr_heights <- rep(0,length(intervals))
for(i in 1:8){
  for(j in 1:length(n_in_bins)){
    rect(xleft = intervals[j], ybottom = curr_heights[j], xright = intervals[j] + binterval, 
         ytop = curr_heights[j] + prop_in_bins_per_condition[i,j], col = cols[i])
    curr_heights[j] <- curr_heights[j] + prop_in_bins_per_condition[i,j]
  }
  rect(xleft = 1.21, ybottom = 0.75 - (i-1) * 0.5 / 8, xright = 1.21 + 0.5 / 8, ytop = 0.75 - i * 0.5 / 8, col = cols[i])
}
for(j in 1:length(n_in_bins)){
  rect(xleft = intervals[j], ybottom = 1.05, xright = intervals[j] + binterval, 
       ytop = 1.05 + prop_in_bins[j] / 2, col = "lightgrey")
  text(x = (intervals + c(rep(0,9),-.01,rep(0,10)))[j] + binterval / 2 , y = 1.05 + prop_in_bins[j] / 2 + 0.025, labels = n_in_bins[j])
}
text(labels = latex2exp::TeX(paste0(" = ", sum(sapply(pairwise_corrs_per_condition, length)) - 144, " = 144^2 - 144")), 
     x = 0.95, y = 1.072, pos = 4)
text(labels = sapply(pairwise_corrs_per_condition, length), x = 1.165 + 0.75 / 8, y = 0.75 - 1:8 * 0.5 / 8 + 0.5 / 16, pos = 4)
text(labels = latex2exp::TeX(paste0(sum(sapply(pairwise_corrs_per_condition, length)), " = 144^2")), x = 1.165 + 0.75 / 8, y = 0.75 - 9 * 0.5 / 8 + 0.5 / 16, pos = 4)
segments(x0 = 1.225 + 0.5 / 8, y0 = 0.25, x1 = 1.35 + 0.5 / 8, y1 = 0.25, lwd = 2)

# plotMatrix(match_conditions + 0, size = c(0.2,0.5), location = c(1,0.25), title = F)
match_conditions_equalsigns <- match_conditions
match_conditions_equalsigns[match_conditions] <- "="; match_conditions_equalsigns[!match_conditions] <- "≠"
plotMatrix(match_conditions_equalsigns, size = c(0.2,0.5), location = c(1,0.25), title = F, rownames = F, colnames = F)
text(x = 1 + c((1:3*2-1)/6)*0.2, y = 0.785, labels = c("Ts", "Tm", "Sx"))
text(x = 1 + c((2*2-1)/6)*0.3, y = 0.85, labels = c("Combination"), cex = 1.5, font = 2)

axis(1, at = seq(-1,1,by=0.2), labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -0.9)
mtext(text = seq(-1,1,by=0.2), side = 1, at = seq(-1,1,by=0.2), cex = 1, line = -0.4)
axis(2, at = seq(0,1,by=0.2), labels =  rep("", 6), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.5)
mtext(text = seq(0,1,by=0.2), side = 2, at = seq(0,1,by=0.2), cex = 1, line = -0.5, las = 2)
title(latex2exp::TeX(paste0("Distribution of Pairwise Spearman's $\\rho$s")), cex.main = 2.5)
mtext(side = 1, text = latex2exp::TeX(paste0("Spearman's $\\rho$")), cex = 1.5, line = 1.25, at = 0)
mtext(side = 2, text = "Bin Proportion", cex = 1.5, line = 1.1, at = 0.5)

#now to get the rest of the figures in at finer granularity
condition_indices <- lapply(1:8, function(x) which(condition_matrix == x, arr.ind = T))

#there's also lots of clever efficient solutions here but whatever, one-offs can be cheap & dirty
if(!exists("already_done")){already_done = F}
if(!already_done){
  arrays <- list()
  arrays_bins <- list()
  for(i in 2:7){
    ncats <- sum(match_conditions[i,])
    a <- array(data = 3, dim = c(c(length(tss), length(tps), length(sex))[match_conditions[i,]], 1E4),
               dimnames = list(tss, tps, sex, paste0("rho_", 1:1E4))[c(match_conditions[i,], T)])
    a_bins <- array(data = 0, dim = c(c(length(tss), length(tps), length(sex))[match_conditions[i,]], length(intervals)),
                    dimnames = list(tss, tps, sex, paste0("bin_", 1:length(intervals)))[c(match_conditions[i,], T)])
    axisnames <- (c(colnames(match_conditions)[match_conditions[i,]], "Rhos"))
    
    if(ncats > 1){
      curr_ind <- a[,,1]; curr_ind[curr_ind == 3] <- 1
    } else if(ncats == 1){
      curr_ind <- a[,1]; curr_ind[curr_ind == 3] <- 1
    }
    
    for(j in 1:nrow(condition_indices[[i]])){
      obs_metadata <- col_df[condition_indices[[i]][j,1],][match_conditions[i,]]
      for(categ in 1:ncats){
        obs_metadata[categ] <- which(as.character(obs_metadata[categ]) == (list(tss, tps, sex)[match_conditions[i,]][[categ]]))
      }
      rho <- zcor[condition_indices[[i]][j,1], condition_indices[[i]][j,2]]
      if(ncats > 1){
        a_bins[as.integer(obs_metadata[1]), as.integer(obs_metadata[2]), max(which(rho > intervals))] <- a_bins[as.integer(obs_metadata[1]), as.integer(obs_metadata[2]), max(which(rho > intervals))] + 1
        a[as.integer(obs_metadata[1]), as.integer(obs_metadata[2]), curr_ind[as.integer(obs_metadata[1]), as.integer(obs_metadata[2])]] <- rho
        curr_ind[as.integer(obs_metadata[1]), as.integer(obs_metadata[2])] <- curr_ind[as.integer(obs_metadata[1]), as.integer(obs_metadata[2])] + 1
      } else if (ncats == 1){
        a_bins[as.integer(obs_metadata), max(which(rho > intervals))] <- a_bins[as.integer(obs_metadata), max(which(rho > intervals))] + 1
        a[as.integer(obs_metadata), curr_ind[as.integer(obs_metadata)]] <- rho
        curr_ind[as.integer(obs_metadata)] <- curr_ind[as.integer(obs_metadata)] + 1
      }
    }
    arrays[[i]] <- a
    arrays_bins[[i]] <- a_bins
  }
  arrays[[1]] <- zcor[condition_matrix == 1]
  arrays_bins[[1]] <- rep(0, length(intervals))
  arrays_bins[[1]][as.integer(names(table(sapply(1:length(arrays[[1]]), function(x) max(which(arrays[[1]][x] > intervals))))))] <- 
    table(sapply(1:length(arrays[[1]]), function(x) max(which(arrays[[1]][x] > intervals))))
  getLastSubArray <- function(a, i){if(length(dim(a)) == 3){return(a[,,i])}else if(length(dim(a)) == 2){return(a[,i])}}
  arrays_total_in_bins <- rbind(arrays_bins[[1]], 
                          t(sapply(2:length(arrays_bins), function(x) sapply(1:length(intervals), function(y) sum(getLastSubArray(arrays_bins[[x]],y))))))
  
  already_done = T
}

####################
#plotting everything
####################

colours_2 <- colours; 
colours_2$Tissue <- colours_2$Tissue[-which(names(colours_2$Tissue) == "t63-testes")]
names(colours_2$Tissue) <- gsub(names(colours_2$Tissue), pattern = "ovaries", replacement = "gonad")
tissues_2 <- sapply(1:length(names(colours_2$Tissue)), function(x) paste0(strsplit(names(colours_2$Tissue)[x], "-")[[1]][-1], collapse = "-"))
colours_2$Sex <- rev(colours_2$Sex)

grDevices::cairo_pdf(filename = "~/Documents/figure2a_suggestion/fig2a_suggestion_image_altOrder.pdf", width = 1090 / 72, height = 1250 / 72)
par(mfrow = c(4,2), mar = c(4,4.5,5,1))
# par(mar = c(4,4.5,5,2))
# lw <- 5; rw <- 4
# layout(matrix(unlist(sapply(1:8, function(x) rep(x, ifelse(x%%2==0,rw,lw)))), byrow = T, nrow = 4))

par(xpd=TRUE)

#plotting first part
plot(1,1,xlim = c(-1,1.5), ylim = c(0,1.25), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
box(which = "figure", lwd = 2, lty = 1.3)
fig_label("a)", cex = 2.5, font = 3, shrinkX = 0.9875, shrinkY = 0.975)
curr_heights <- rep(0,length(intervals))
for(i in 1:8){
  for(j in 1:length(n_in_bins)){
    rect(xleft = intervals[j], ybottom = curr_heights[j], xright = intervals[j] + binterval, 
         ytop = curr_heights[j] + prop_in_bins_per_condition[i,j], col = cols[i])
    curr_heights[j] <- curr_heights[j] + prop_in_bins_per_condition[i,j]
  }
  rect(xleft = 1.21 - 0.025, ybottom = 0.75 - 0.175 - (i-1) * 0.5 / 8, xright = 1.21 + 0.5 / 8 - 0.025, ytop = 0.75 - 0.175 - i * 0.5 / 8, col = cols[i])
}
for(j in 1:length(n_in_bins)){
  rect(xleft = intervals[j], ybottom = 1.05, xright = intervals[j] + binterval, 
       ytop = 1.05 + prop_in_bins[j] / 1.75, col = "lightgrey")
  text(x = (intervals + c(rep(0,9),-.01,rep(0,10)))[j] + binterval / 2 , y = 1.05 + prop_in_bins[j] / 1.75 + 0.025, labels = n_in_bins[j])
}
# text(labels = latex2exp::TeX(paste0(" = ", (sum(sapply(pairwise_corrs_per_condition, length)) - 144) / 2, " = \\frac{144^2 - 144}{2} = $\\left(^{144}_{ 2}\\right)$")), 
#      x = 0.95, y = 1.06, pos = 4)
text(labels = paste0(" = ", (sum(sapply(pairwise_corrs_per_condition, length)) - 144) / 2, " = "), x = 0.95, y = 1.071, pos = 4)
text(labels = latex2exp::TeX(paste0("144^2 - 144")), x = 1.1725, y = 1.1, pos = 4)
text(labels = latex2exp::TeX(paste0("2")), x = 1.2675, y = 1.035, pos = 4)
segments(x0 = 1.205,y0 = 1.07,x1 = 1.4,y1 = 1.07)
text(labels = "=", x = 1.39, y = 1.07, pos = 4)
text(labels = latex2exp::TeX(paste0("$\\left(\\frac{144}{2}\\right)$")), x = 1.425, y = 1.045, pos = 4)
segments(x0 = 1.47,y0 = 1.07,x1 = 1.54,y1 = 1.07, col = "white", lwd = 3)

# text(labels = sapply(pairwise_corrs_per_condition, length) / 2, x = 1.165 + 0.75 / 8 - 0.025, y = 0.75 - 0.175 - 1:8 * 0.5 / 8 + 0.5 / 16, pos = 4)
pcpc <- as.character(c(sapply(pairwise_corrs_per_condition, length)) / 2 ); pcpc[length(pcpc)] <- "(144)"
text(labels = pcpc, x = 1.165 + 0.75 / 8 - 0.025, y = 0.75 - 0.175 - 1:8 * 0.5 / 8 + 0.5 / 16, pos = 4)
# text(labels = latex2exp::TeX(paste0(sum(sapply(pairwise_corrs_per_condition, length)) / 2, " = 144^2")), x = 1.165 + 0.75 / 8 - 0.025, y = 0.75 - 0.175 - 9 * 0.5 / 8 + 0.5 / 16, pos = 4)
text(labels = latex2exp::TeX(paste0(sum(sapply(pairwise_corrs_per_condition, length)) / 2, " = ")), x = 1.165 + 0.75 / 8 - 0.025, y = 0.74 - 0.175 - 9 * 0.5 / 8 + 0.5 / 16, pos = 4)
text(labels = latex2exp::TeX(paste0("$\\left(\\frac{144}{2}\\right)$")), x = 1.4075, y = 0.012, pos = 4)
segments(x0 = 1.4525,y0 = 0.04,x1 = 1.522,y1 = 0.04, col = "white", lwd = 4)
text(labels = "+", x = 1.4475, y = -0.06, pos = 4)
text(labels = "(144)", x = 1.41, y = -0.115, pos = 4)

segments(x0 = 1.225 + 0.5 / 8 - 0.025, y0 = 0.25 - 0.18, x1 = 1.35 + 0.5 / 8 - 0.025, y1 = 0.25 - 0.18, lwd = 1.5)
# segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, lwd = 2)

# plotMatrix(match_conditions + 0, size = c(0.2,0.5), location = c(1,0.25), title = F)
match_conditions_equalsigns <- match_conditions
match_conditions_equalsigns[match_conditions] <- "="; match_conditions_equalsigns[!match_conditions] <- "≠"
plotMatrix(match_conditions_equalsigns, size = c(0.2,0.5), location = c(0.975,0.075), title = F, rownames = F, colnames = F, cex = 1.25)
text(x = 0.975 + c((1:3*2-1)/6)*0.2, y = 0.61, labels = c("Ts", "Tm", "Sx"))
text(x = 1 + c((2*2-1)/6)*0.3, y = 0.675, labels = c("Combination"), cex = 1.5, font = 2)

axis(1, at = seq(-1,1,by=0.2), labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -0.9)
mtext(text = seq(-1,1,by=0.2), side = 1, at = seq(-1,1,by=0.2), cex = 1, line = -0.2)
axis(2, at = seq(0,1,by=0.2), labels =  rep("", 6), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.5)
mtext(text = seq(0,1,by=0.2), side = 2, at = seq(0,1,by=0.2), cex = 1, line = -0.75, las = 2)
title(latex2exp::TeX(paste0("Distribution of Pairwise Spearman's $\\rho$s")), cex.main = 2.5, adj = 0.35)
mtext(side = 1, text = latex2exp::TeX(paste0("Spearman's $\\rho$")), cex = 1.5, line = 2.25, at = 0)
mtext(side = 2, text = "Bin Proportion", cex = 1.5, line = 2, at = 0.5)

#plotting the other 7 figures
for(fig in c(2,4,6,1,3,5,7)){
# for(fig in 1:7){
  # if(fig %% 2 == 1){
  #   par(mar = c(4,4.5,5,0))
  # } else{
  #   par(mar = c(4,4.5,5,2))
  # }
  plot(1,1,xlim = c(-1,1.5), ylim = c(0,1.25), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  fig_label(paste0(letters[1+fig], ")"), cex = 2.5, font = 3, shrinkX = 0.9875, shrinkY = 0.975)
  box(which = "figure", lwd = 2, lty = 1.3)
  curr_heights <- rep(0,length(intervals))
  cols_2 <- colours_2[match_conditions[fig,]]
  if(fig == 1){text(x = 0, y = 0.5, labels = "(no finer-grained breakdown possible)", cex = 1.5)}
  if(length(cols_2) == 1){
    for(i in 1:length(cols_2[[1]])){
      for(j in 1:length(intervals)){
        rect(xleft = intervals[j], ybottom = curr_heights[j], xright = intervals[j] + binterval, 
             ytop = curr_heights[j] + arrays_bins[[fig]][i,j] / arrays_total_in_bins[fig,j], col = cols_2[[1]][i])
        curr_heights[j] <- curr_heights[j] + arrays_bins[[fig]][i,j] / arrays_total_in_bins[fig,j]
      }
      # rect(xleft = 1.21, ybottom = 0.75 - (i-1) * 0.5 / 8, xright = 1.21 + 0.5 / 8, ytop = 0.75 - i * 0.5 / 8, col = cols[i])
    }
  } else if(length(cols_2) == 2){
    shading_angles = seq(from = 45, to = 135, length.out = length(cols_2[[2]]))
    for(i1 in 1:length(cols_2[[1]])){
      for(i2 in 1:length(cols_2[[2]])){
        for(j in 1:length(intervals)){
          rect(xleft = intervals[j], ybottom = curr_heights[j], xright = intervals[j] + binterval, 
               ytop = curr_heights[j] + arrays_bins[[fig]][i1, i2, j] / arrays_total_in_bins[fig,j], col = cols_2[[1]][i1])
          rect(xleft = intervals[j], ybottom = curr_heights[j], xright = intervals[j] + binterval, 
               ytop = curr_heights[j] + arrays_bins[[fig]][i1, i2, j] / arrays_total_in_bins[fig,j], 
               col = cols_2[[2]][i2], border = "transparent", density = 10, angle = shading_angles[i2], lwd = 1.25)
          curr_heights[j] <- curr_heights[j] + arrays_bins[[fig]][i1, i2, j] / arrays_total_in_bins[fig,j]
        }
        # rect(xleft = 1.21, ybottom = 0.75 - (i-1) * 0.5 / 8, xright = 1.21 + 0.5 / 8, ytop = 0.75 - i * 0.5 / 8, col = cols[i])
      }
    }
  }
  
  for(j in 1:length(n_in_bins)){
    rect(xleft = intervals[j], ybottom = 1.05, xright = intervals[j] + binterval, 
         ytop = 1.05 + arrays_total_in_bins[fig,j] / sum(arrays_total_in_bins[fig,]) / max(arrays_total_in_bins[fig,] / sum(arrays_total_in_bins[fig,])) * 0.2, 
         col = grDevices::adjustcolor(cols[fig], alpha.f = 0.75))
    text(x = (intervals + c(rep(0,9),-.01,rep(0,10)))[j] + binterval / 2 , 
         y = 1.05 + arrays_total_in_bins[fig,j] / sum(arrays_total_in_bins[fig,]) / max(arrays_total_in_bins[fig,] / sum(arrays_total_in_bins[fig,])) * 0.2 + 0.025, 
         labels = arrays_total_in_bins[fig,j] / 2)
  }
  
  text(labels = latex2exp::TeX(paste0(" = ", sum(arrays_total_in_bins[fig,]) / 2)),
       x = 0.95, y = 1.07, pos = 4)
  
  #plot legends
  if(fig > 1){
    if(names(cols_2)[1] == "Tissue"){
      plotMatrix(t(t(list(tissues_2, tps, sex)[match_conditions[fig,]][[1]])), size = c(0.3,0.8), 
                 location = c(intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075,0.075), title = F, rownames = F, colnames = F, cex = 0.9)
      text(x = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.175, y = 0.925, labels = c("Tissue"), cex = 1.5, font = 2)
      for(bx in 1:length(cols_2$Tissue)){
        rect(xleft = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.011 + 0.3, 
             ybottom = 0.875 - (bx-1) * 0.8 / length(cols_2$Tissue), 
             xright = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.8 / length(cols_2$Tissue) + 0.011 + 0.3, 
             ytop = 0.875 - bx * 0.8 / length(cols_2$Tissue), col = cols_2$Tissue[bx])
      }
    }
    if(names(cols_2)[1] == "Time"){
      plotMatrix(t(t(list(tissues_2, tps, sex)[match_conditions[fig,]][[1]])), size = c(0.075,0.25), 
                 location = c(intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075,0.075), title = F, rownames = F, colnames = F, cex = 1)
      text(x = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.135, y = 0.375, labels = c("Time"), cex = 1.5, font = 2)
      for(bx in 1:length(cols_2$Time)){
        rect(xleft = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.011 + 0.075, 
             ybottom = 0.325 - (bx-1) * 0.25 / length(cols_2$Time), 
             xright = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.25 / length(cols_2$Time) + 0.011 + 0.075, 
             ytop = 0.325 - bx * 0.25 / length(cols_2$Time), col = cols_2$Time[bx])
      }
    }
    if(names(cols_2)[1] == "Sex"){
      plotMatrix(t(t(list(tissues_2, tps, sex)[match_conditions[fig,]][[1]])), size = c(0.15,0.125), 
                 location = c(intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075,0.075), title = F, rownames = F, colnames = F, cex = 1)
      text(x = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.135, y = 0.25, labels = c("Sex"), cex = 1.5, font = 2)
      for(bx in 1:length(cols_2$Sex)){
        rect(xleft = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.011 + 0.15, 
             ybottom = 0.2 - (bx-1) * 0.125 / length(cols_2$Sex), 
             xright = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.125 / length(cols_2$Sex) + 0.011 + 0.15, 
             ytop = 0.2 - bx * 0.125 / length(cols_2$Sex), col = cols_2$Sex[bx])
      }
    }
    if(length(cols_2) > 1){
      if(names(cols_2)[2] == "Time"){
        horiz_disp <- 0.4
        plotMatrix(t(t(list(tissues_2, tps, sex)[match_conditions[fig,]][[2]])), size = c(0.075,0.25), 
                   location = c(intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + horiz_disp,0.075), title = F, rownames = F, colnames = F, cex = 1)
        text(x = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.135 + horiz_disp, y = 0.375, labels = c("Time"), cex = 1.5, font = 2)
        for(bx in 1:length(cols_2$Time)){
          rect(xleft = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.011 + 0.075 + horiz_disp, 
               ybottom = 0.325 - (bx-1) * 0.25 / length(cols_2$Time), 
               xright = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.25 / length(cols_2$Time) + 0.011 + 0.075 + horiz_disp, 
               ytop = 0.325 - bx * 0.25 / length(cols_2$Time), 
               col = cols_2[[2]][bx], border = "black", density = 30, angle = shading_angles[bx], lwd = 1.25)
        }
      }
      if(names(cols_2)[2] == "Sex"){
        horiz_disp <- ifelse(names(cols_2)[1] == "Tissue", 0.4, 0.1875)
        plotMatrix(t(t(list(tissues_2, tps, sex)[match_conditions[fig,]][[2]])), size = c(0.15,0.125), 
                   location = c(intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + horiz_disp,0.075), title = F, rownames = F, colnames = F, cex = 1)
        text(x = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.135 + horiz_disp, y = 0.25, labels = c("Sex"), cex = 1.5, font = 2)
        for(bx in 1:length(cols_2$Sex)){
          rect(xleft = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.011 + 0.15 + horiz_disp, 
               ybottom = 0.2 - (bx-1) * 0.125 / length(cols_2$Sex), 
               xright = intervals[max(which(arrays_total_in_bins[fig,] != 0)) + 1] + 0.075 + 0.125 / length(cols_2$Sex) + 0.011 + 0.15 + horiz_disp, 
               ytop = 0.2 - bx * 0.125 / length(cols_2$Sex), 
               col = cols_2[[2]][bx], border = "black", density = 30, angle = shading_angles[bx], lwd = 1.25)
        }
      } 
    }
  }
  # text(labels = latex2exp::TeX(paste0(" = ", sum(arrays_total_in_bins[fig,]), " = ", length(cols_2[[1]]), " x ", length(cols_2[[2]]))), 
  #      x = 0.95, y = 1.072, pos = 4)
  # text(labels = sapply(pairwise_corrs_per_condition, length), x = 1.165 + 0.75 / 8, y = 0.75 - 1:8 * 0.5 / 8 + 0.5 / 16, pos = 4)
  # text(labels = latex2exp::TeX(paste0(sum(sapply(pairwise_corrs_per_condition, length)), " = 144^2")), x = 1.165 + 0.75 / 8, y = 0.75 - 9 * 0.5 / 8 + 0.5 / 16, pos = 4)
  # segments(x0 = 1.225 + 0.5 / 8, y0 = 0.25, x1 = 1.35 + 0.5 / 8, y1 = 0.25, lwd = 2)
  # 
  # # plotMatrix(match_conditions + 0, size = c(0.2,0.5), location = c(1,0.25), title = F)
  # match_conditions_equalsigns <- match_conditions
  # match_conditions_equalsigns[match_conditions] <- "="; match_conditions_equalsigns[!match_conditions] <- "≠"
  # plotMatrix(match_conditions_equalsigns, size = c(0.2,0.5), location = c(1,0.25), title = F, rownames = F, colnames = F)
  # text(x = 1 + c((1:3*2-1)/6)*0.2, y = 0.785, labels = c("Ts", "Tm", "Sx"))
  # text(x = 1 + c((2*2-1)/6)*0.3, y = 0.85, labels = c("Combination"), cex = 1.5, font = 2)
  
  axis(1, at = seq(-1,1,by=0.2), labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -0.9)
  mtext(text = seq(-1,1,by=0.2), side = 1, at = seq(-1,1,by=0.2), cex = 1, line = -0.2)
  axis(2, at = seq(0,1,by=0.2), labels =  rep("", 6), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.5)
  mtext(text = seq(0,1,by=0.2), side = 2, at = seq(0,1,by=0.2), cex = 1, line = -0.75, las = 2)
  title(latex2exp::TeX(paste0("Spearman's $\\rho$s")), cex.main = 2.5, adj = 0.375, col.main = cols[fig], line = 3.25)
  title(paste0("\n(", match_conditions_equalsigns[fig,1], " Tissue, ", 
                    match_conditions_equalsigns[fig,2], " Timepoint, ", 
                    match_conditions_equalsigns[fig,3], " Sex)"), cex.main = 1.5, font = 3, adj = 0.35)
  mtext(side = 1, text = latex2exp::TeX(paste0("Spearman's $\\rho$")), cex = 1.5, line = 2.25, at = 0)
  mtext(side = 2, text = "Bin Proportion", cex = 1.5, line = 2, at = 0.5)
  
}

dev.off()

#sanity check
which(zcor == min(zcor[condition_matrix == 4]), arr.ind = T)

#check tissue colors
plot(rep(1,19), 1:19, col = "white"); text(x = rep(1, 19), y = 1:19, labels = names(colours$Tissue), col = colours$Tissue)


#### hive plots ####
# library(HiveR)
# 
# col_df = data.frame(row.names=colnames(zcor))
# col_df$Tissue = gsub('_.*','',rownames(col_df))
# col_df$Time = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[3])
# col_df$Sex = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[2])
# threshold_correlation <- 0.3
# adjMat <- (abs(zcor) > threshold_correlation) + 0
# rownames(adjMat) <- colnames(adjMat) <- col_df$Tissue
# sapply(names(colours$Tissue), function(rowT) sapply(names(colours$Tissue), function(colT) sum(adjMat[rowT, colT])))
# 
# edges_tissue <- as.data.frame(do.call(rbind, sapply(1:(nrow(zcor)-1), function(ri) t(sapply((ri+1):(ncol(zcor)), function(ci) 
#   c(col_df$Tissue[ri], col_df$Tissue[ci], zcor[ri,ci], (sign(zcor[ri,ci]) + 1) / 2 + 1)
#   )))))
# colnames(edges_tissue) <- c("node1", "node2", "weight", "color")
# edges_tissue$weight <- abs(as.numeric(edges_tissue$weight))
# edges_tissue$weight <- edges_tissue$weight / max(edges_tissue$weight) * 2
# source("~/scripts/hive_plots/mod.edge2HPD.R")
# HPD <- mod.edge2HPD(edge_df = edges_tissue[,1:2], edge.weight = edges_tissue$weight, edge.color = edges_tissue$color,
#                     type = "2D", axis.cols = colours$Tissue, 
#                     node.axis = cbind(edges_tissue$node1,match(edges_tissue$node1, names(colours$Tissue))))

#ok... this is dumb and not working, let's just make one from scratch
library(pracma)
polarp <- function (t, r, col = 1, pch = 1, cex = 1, center = c(0,0)) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  points(xy[, 1] + center[1], xy[, 2] + center[2], cex = cex, col = col, pch = pch)
}
polarl <- function (t, r, col = 1, lwd = 1, center = c(0,0)) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  lines(xy[, 1] + center[1], xy[, 2] + center[2], lwd = lwd, col = col)
}
arc <- function(t1,t2,r1,r2,res=1E2,lwd = 1,col=1, mindist = T, self_adjust = 2 * pi / 45, random_selfing = T, clockwise_selfing = T, pointy_selfing = F,
                center = c(0,0)){
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
  polarl(ts, rs, lwd = lwd,col=col, center = center)
}
line <- function(t,r1,r2,lwd = 1,col=1, center = c(0,0)){
  polarl(rep(t, 2), c(r1,r2), lwd = lwd,col=col, center = center)
}


#tissues on spokes
et <- as.data.frame(do.call(rbind, sapply(1:(nrow(zcor)-1), function(ri) t(sapply((ri+1):(ncol(zcor)), function(ci) 
  c(tiss1 = col_df$Tissue[ri], tiss2 = col_df$Tissue[ci], corr = zcor[ri,ci], color = (sign(zcor[ri,ci]) + 1) / 2 + 1,
    sex1 = which(col_df$Sex[ri] == c("male", "female")), sex2 = which(col_df$Sex[ci] == c("male", "female")), 
    time1 = which(col_df$Time[ri] == paste0(c(1,2,4,8), "w")), time2 = which(col_df$Time[ci] == paste0(c(1,2,4,8), "w")))
)))))

axis_length <- 1
et$weight <- abs(as.numeric(et$corr))
et <- et[et$weight > 0.2,]
n_tiss <- length(unique(et$tiss1))
aas <- seq(0, 2*pi - 2*pi / n_tiss, length.out = n_tiss)
rs <-seq(axis_length/8, axis_length, length.out = 8)
et$theta1 <- aas[match(et$tiss1, names(colours$Tissue))]
et$theta2 <- aas[match(et$tiss2, names(colours$Tissue))]
et$r1 <- ((as.numeric(et$sex1)-1)*4 + as.numeric(et$time1 )) / 8
et$r2 <- ((as.numeric(et$sex2)-1)*4 + as.numeric(et$time2 )) / 8
et$color <- c("#2096be", "#e28a4a")[as.numeric(et$color)]

plot(1,1,xlim = c(-1,1), ylim = c(-1,1), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
for(ri in 1:nrow(et)){
  arc(t1 = et[ri,"theta1"], t2 = et[ri,"theta2"], r1 = et[ri,"r1"], r2 = et[ri,"r2"], lwd = et$weight[ri]^1.5 * 5, col = et$color[ri], random_selfing = T)
}
# ri = 10
# arc(t1 = et[ri,"theta1"], t2 = et[ri,"theta2"], r1 = et[ri,"r1"], r2 = et[ri,"r2"], lwd = et$weight[ri]^1.5 * 5, col = et$color[ri])

for(ax in 1:n_tiss){line(t = aas[ax], 0, 1, col = colours$Tissue[ax], lwd = 4)}

# nice_names <- sapply(names(colours$Tissue), function(tissue) substr(paste0(strsplit(tissue, "-")[[1]][-c(1,3)], collapse = " "), start = 1, stop = 10))

nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])
nice_names <- sapply(names(colours$Tissue), function(tissue) nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)])

for(tissue in 1:n_tiss){
  text_loc <- pol2cart(cbind(aas[tissue], 1 + (sin(aas[tissue] - pi) + 1)/30))
  text(x = text_loc[1], y = text_loc[2], labels = nice_names[tissue], srt = aas[tissue] / 2 / pi * 360 - 90, pos = 3, col = colours$Tissue[tissue])
  for(sex in 1:2){for(time in 1:4){
  polarp(t = aas[tissue], r = ((as.numeric(sex)-1)*4 + as.numeric(time)) / 8, col = colours$Sex[sex], pch = 19, cex = 1.4)  
  polarp(t = aas[tissue], r = ((as.numeric(sex)-1)*4 + as.numeric(time)) / 8, col = colours$Time[time], pch = 19, cex = 1)
  }}
}

xl <- -1.15; yb <- 0.5; xr <- -1.15; yt <- 1;
lws <- (c(0.3, 1) * 15 - 4) / 72 / par("pin")[1]
hoffset <- 0
text(labels = round(seq(-1, 1, length.out = 11), 2), y = seq(yb, yt, length.out = 11), x = xr - 0.0075, pos = 4, las=2, cex=0.9)
text(labels = latex2exp::TeX(paste0("$\\rho$")), y = yt - 0.01, x = (xl + xr) / 2, pos = 3, cex = 2)
polygon(x = c(xl - hoffset - lws[1]/2, 
              xl - hoffset - lws[2]/2,
              xl - hoffset + lws[2]/2, 
              xl - hoffset + lws[1]/2),
        y = c((yt + yb) / 2, yb,yb,(yt + yb) / 2), col = "#2096be")
polygon(x = c(xl - hoffset - lws[2]/2,
              xl - hoffset - lws[1]/2,
              xl - hoffset + lws[1]/2,
              xl - hoffset + lws[2]/2), 
        y = c(yt, (yt + yb) / 2,(yt + yb) / 2,yt), col = "#e28a4a")
points(x = c(-1, -1), y = c(0.95, 1), col = colours$Sex, cex = 1.4)
text(x = c(-1, -1), y = c(0.95, 1), labels = c("male", "female"), pos = 4)
points(x = rep(-1,4), y = c(0.75,0.8,0.85,0.9), col = colours$Time, cex = 1.2, pch = 19)
text(x = rep(-1,4), y = c(0.75,0.8,0.85,0.9), labels = paste0(c(1,2,4,8), "w"), pos = 4, cex = 0.9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### let's summon cthulhu ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

colours = list(Tissue=tissue_cols[unique(col_df$Tissue)], 
               Time=group_cols[unique(col_df$Time)],
               Sex=sex_cols[c('male','female')])
colours$Tissue["t56-vastus-lateralis"] <- colorRampPalette(c(as.character(colours$Tissue["t56-vastus-lateralis"]), "black"))(3)[2]
nnmap <- nnmap[nnmap$tissue_name_release != "",]
tiss_ord <- nnmap$tissue_name_release[match(rev(MotrpacBicQC::tissue_order), nnmap$abbreviation)]
colours$Tissue <- colours$Tissue[match(tiss_ord, names(colours$Tissue))]
colours$Tissue <- colours$Tissue[!is.na(colours$Tissue)]

#figure parameters
nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])
nice_names <- sapply(names(colours$Tissue), function(tissue) nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)])

corr_thresh <- 0.3

axis.length <- 1.5
center_rescaler <- 1.25
inner_shifter <- 0.96
no_concentric_arcs <- F
opacity_concentric_arcs <- 1
opacity_nonconcentric_arcs <- 0.25
outer_shifter <- 1.04
line_weight_power <- 2
line_weight_multiplier <- 10
numbers_in_squares <- F
tissue_names_not_colors <- T
tissue_name_cex <- 0.29
across_relationships <- T
adjacent_relationships <- T

# for(tissues_to_include in tissues){

tissues_to_include <- tissues[2]
print(tissues_to_include)
  
grDevices::cairo_pdf(filename = paste0("~/Documents/figure2a_suggestion/network_visualization/", 
                                       ifelse(all(tissues %in% tissues_to_include), "all-tissues", tissues_to_include),
                                       "_corrThresh-", corr_thresh, ".pdf"), 
                                       width = 650 / 72, height = 585 / 72, family="Arial Unicode MS")
par(mar = c(6,9,5,6), xpd = NA)

if(tissue_names_not_colors){
  nice_names[nice_names == "LUNG"] <- "LUNGS"
  nice_names[nice_names == "BAT"] <- "BR-AT"
  inner_shifter <- 0.9
  outer_shifter <- 1.105
}

et <- as.data.frame(do.call(rbind, sapply(1:(nrow(zcor)-1), function(ri) t(sapply((ri+1):(ncol(zcor)), function(ci) 
  c(tiss1 = col_df$Tissue[ri], tiss2 = col_df$Tissue[ci], corr = zcor[ri,ci], color = (sign(zcor[ri,ci]) + 1) / 2 + 1,
    sex1 = which(col_df$Sex[ri] == c("male", "female")), sex2 = which(col_df$Sex[ci] == c("male", "female")), 
    time1 = which(col_df$Time[ri] == paste0(c(1,2,4,8), "w")), time2 = which(col_df$Time[ci] == paste0(c(1,2,4,8), "w")))
)))))

et$weight <- abs(as.numeric(et$corr))
et <- et[et$weight > corr_thresh,]
n_tiss <- length(unique(c(et$tiss1, et$tiss2)))
aas <- list(c(pi/4 + 1E-3, 5*pi/4 - 1E-3), 
            c(3*pi/4 - 1E-3, 7*pi/4 + 1E-3), 
            rev(c(pi/4 - 1E-3, 5*pi/4 + 1E-3)), 
            rev(c(3*pi/4 + 1E-3, 7*pi/4 - 1E-3)))
rs <- 1:(n_tiss)*axis.length/(n_tiss)
  
# seq(axis.length/(n_tiss)/2, axis.length, by = axis.length/n_tiss)
# rs <- rs - rs[1] + diff(rs)[1]/2

et$theta1 <- match(et$sex1, c(1,2))
et$theta2 <- match(et$sex2, c(1,2))
et$r1 <- rs[match(et$tiss1, names(colours$Tissue))]
et$r2 <- rs[match(et$tiss2, names(colours$Tissue))]
et$color <- c("#2096be", "#e28a4a")[as.numeric(et$color)]

et$concentric <- et$tiss1 != et$tiss2
et$opacity <- opacity_nonconcentric_arcs
et$opacity[et$concentric] <- opacity_concentric_arcs

et_wt <- lapply(1:4, function(time) et[et$time1 == time & et$time2 == time & (et$tiss1 %in% tissues_to_include | et$tiss2 %in% tissues_to_include),])
# et_wt <- lapply(et_wt, function(et_sub) et_sub[1:3,])
if(no_concentric_arcs){
  et_wt <- lapply(1:4, function(time) et_wt[[time]][et_wt[[time]]$tiss1 != et_wt[[time]]$tiss2,])
}

plot(1,1,xlim = c(-2,2), ylim = c(-2,2), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
centers <- list(c(-1,1), c(1,1), c(1,-1), c(-1,-1))
centers <- lapply(centers, function(x) x * center_rescaler)
for(time in 1:4){
  if(nrow(et_wt[[time]]) == 0){next()}
  for(ri in 1:nrow(et_wt[[time]])){
    arc(t1 = aas[[time]][et_wt[[time]][ri,"theta1"]], t2 = aas[[time]][et_wt[[time]][ri,"theta2"]], r1 = et_wt[[time]][ri,"r1"], center = centers[[time]] * outer_shifter,
        r2 = et_wt[[time]][ri,"r2"], lwd = et_wt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier, col = adjustcolor(et_wt[[time]]$color[ri], et_wt[[time]]$opacity[ri]), 
        random_selfing = F, clockwise_selfing = (et_wt[[time]][ri,"sex1"] == "1") * c(1,-1,1,-1)[time]*-1 + rev(c(0,1,0,1))[time], self_adjust = 0.3)
  }
  for(sex in 1:2){
    # line(t = aas[[time]][sex], 0, axis.length, col = colours$Sex[sex], lwd = 4, center = centers[[time]])
    # for(tissue in 1:n_tiss){
    #   polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Tissue[tissue], pch = 19, cex = 0.85, center = centers[[time]])
    # }
  }
}

max_outer_radii <- sapply(et_wt, function(etwt) max(etwt$r1 + etwt$r2) / 2)
week_label_locs <- sapply(1:4, function(w) polar2cart(r = max_outer_radii[w], t = (-(0:3)*pi/2 + 3*pi/4)[w]) + centers[[w]])
week_label_locs[week_label_locs == Inf] <- -center_rescaler * 1.2
week_label_locs[week_label_locs == -Inf] <- center_rescaler * 1.2
timelab_nudge <- 0.25

for(i in 1:4){
  text(c("\u2642", "\u2640", "\u2642", "\u2640")[i], col = c(colours$Sex, colours$Sex)[i], 
       x = c(0,1,0,-1)[i]*(2 - sqrt(2) + axis.length)*1.15 - c(2,0,-2.15,0)[i] * 0.06, 
       y = c(1,0,-1,0)[i]*(2 - sqrt(2) + axis.length)*1.275 - c(-1,-1,-1,-1)[i] * 0.06, 
       cex = 4.5, srt = c(45,90,225,270)[i], pos = 1)
  # shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = c(-1,1,1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, 
  #            y = c(1,1,-1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, cex = 4,
  #            srt = c(45,-45,235,135)[i], pos = 1)
  
  #plot weeks
  shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = week_label_locs[1,i] + c(-1,1,1,-1)[i]*timelab_nudge, 
             y = week_label_locs[2,i] + c(1,1,-1,-1)[i]*timelab_nudge, cex = 4,
             srt = c(45,-45,235,135)[i], pos = 1)
  
}

#two / to is where, 1: counter-clockwise, 2: clockwise, 3: opposite
et$tiw_diff <- (as.numeric(et$time2) - as.numeric(et$time1))
et$tiw <- NA
et$tiw[et$tiw_diff == c(-1, 3)[1] | et$tiw_diff == c(-1, 3)[2]] <- 1
et$tiw[et$tiw_diff == c(1, -3)[1] | et$tiw_diff == c(1, -3)[2]] <- 2
et$tiw[et$tiw_diff == c(-2, 2)[1] | et$tiw_diff == c(-2, 2)[2]] <- 3
et$tiw[et$tiw_diff == 0] <- 0

centers_bt <- list(list(c(-2,0), c(0,2)),
                   list(c(0,2), c(2,0)),
                   list(c(2,0), c(0,-2)),
                   list(c(0,-2), c(-2,0)))
centers_bt <- lapply(centers_bt, function(x1) lapply(x1, function(x2) x2 * center_rescaler * inner_shifter))

#specify angles between new axes
aas_bt <- list(list(c(1*pi/4, 7*pi/4), c(5*pi/4, 7*pi/4)),
            list(c(5*pi/4, 7*pi/4), c(3*pi/4,5*pi/4)),
            list(c(3*pi/4,5*pi/4), c(1*pi/4,3*pi/4)),
            list(c(1*pi/4,3*pi/4), c(1*pi/4, 7*pi/4)))

et$theta1 <- ifelse(et$tiw == 1, 2, 1)
et$theta2 <- ifelse(et$theta1 == 1, 2, 1)


rs_bt <- seq(sqrt(2) * center_rescaler * inner_shifter - axis.length, 
             axis.length + sqrt(2) * center_rescaler * inner_shifter, 
             length.out = n_tiss * 2 + 1)[-(n_tiss+1)]
# rs_bt <- c(rs, rs + max(rs) + diff(rs)[1]) - diff(rs)[1] + sqrt(2) * center_rescaler - axis_length
# polarp(r = rs_bt, t = rep(pi/4, length(rs_bt)), center = c(-2.5,-0), cex = , pch = 1, col = "black")



# rs_bt <- seq(sqrt(2) * center_rescaler - max(rs), 
#              sqrt(2) * center_rescaler, 
#              length.out = n_tiss)
# rs_bt <- c(rs_bt, rs_bt + max(rs_bt) - min(rs_bt) + diff(rs_bt)[1])

# rs_bt <- rs_bt - rs_bt[1] + diff(rs_bt)[1]

et$tpairs <- sapply(1:nrow(et), function(ri) paste0(sort(c(et$time1[ri], et$time2[ri])), collapse = ""))
et$close_sex <- NA
et$close_sex[et$tpairs == "12"] <- 1
et$close_sex[et$tpairs == "23"] <- 2
et$close_sex[et$tpairs == "34"] <- 1
et$close_sex[et$tpairs == "14"] <- 2

# et$r1 <- rs_bt[match(et$tiss1, names(colours$Tissue)) + (as.numeric(et$sex1) - 1) * n_tiss]
#gets index in rs_bt for close / far sex
t1_along <- cbind((n_tiss + 1) - match(et$tiss1, names(colours$Tissue)), n_tiss + match(et$tiss1, names(colours$Tissue))) 
et$r1 <- rs_bt[sapply(1:nrow(t1_along), function(ri) t1_along[ri, 2 - as.numeric(et$sex1[ri] == et$close_sex[ri])])]
# et$r2 <- rs_bt[match(et$tiss2, names(colours$Tissue)) + (as.numeric(et$sex2) - 1) * n_tiss]
t2_along <- cbind((n_tiss + 1) - match(et$tiss2, names(colours$Tissue)), n_tiss + match(et$tiss2, names(colours$Tissue)))
et$r2 <- rs_bt[sapply(1:nrow(t2_along), function(ri) t2_along[ri, 2 - as.numeric(et$sex2[ri] == et$close_sex[ri])])]
et$concentric <- et$tiss1 != et$tiss2 | et$sex1 != et$sex2
et$opacity <- opacity_nonconcentric_arcs
et$opacity[et$concentric] <- opacity_concentric_arcs
 

et_bt <- lapply(1:4, function(time) et[et$time1 == time & et$time2 != time & (et$tiss1 %in% tissues_to_include | et$tiss2 %in% tissues_to_include),])
if(no_concentric_arcs){
  et_bt <- lapply(1:4, function(time) et_bt[[time]][et_bt[[time]]$tiss1 != et_bt[[time]]$tiss2 | et_bt[[time]]$sex1 != et_bt[[time]]$sex2,])
}

#hack to fix sex mixup -- TODO find originl bug
# for(time in 1:4){
#   temp <- et_bt[[time]]$sex1
#   et_bt[[time]]$sex1 <- et_bt[[time]]$sex2
#   et_bt[[time]]$sex2 <- temp
# }

# et_bt <- lapply(et_bt, function(et_sub) et_sub[1:10,])
# et_bt <- lapply(1:4, function(et_sub) et_bt[[et_sub]])

if(adjacent_relationships){

  for(time in 1:4){
    for(ri in 1:nrow(et_bt[[time]])){
      if(any(et_bt[[time]][ri,"tiw"] == c(1,2))){
        t1 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta1"]]
        t2 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta2"]]
        r1 = et_bt[[time]][ri,"r1"]
        r2 = et_bt[[time]][ri,"r2"]
        
        #hack to get around 1w and 8w axes switching -- should probs find more principled solution sometime
        if(all(sort(c(et_bt[[time]]$time1[ri], et_bt[[time]]$time2[ri])) == c(1,4))){
          tt <- t1; t1 <- t2; t2 <- tt
          # rt <- r1; r1 <- r2; r2 <- rt
        }
        
        arc(t1 = t1,
            t2 = t2,
            r1 = r1,
            r2 = r2,
            center = centers_bt[[time]][et_bt[[time]][ri,"tiw"]][[1]],
            lwd = et_bt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier,
            col = adjustcolor(et_bt[[time]]$color[ri], et_bt[[time]]$opacity[ri])
            )
      }
    }
  }

}




# time = 1
# ri = 1
# arc(t1 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta1"]], 
#     t2 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta2"]], 
#     r1 = et_bt[[time]][ri,"r1"], 
#     r2 = et_bt[[time]][ri,"r2"],
#     center = centers_bt[[time]][et_bt[[time]][ri,"tiw"]][[1]],
#     lwd = et_bt[[time]]$weight[ri]^1.5 * 5, 
#     col = adjustcolor(et_bt[[time]]$color[ri], 0.85)
# )

if(!tissue_names_not_colors){
for(time in 1:4){
  for(sex in 1:2){
    # line(t = aas[[time]][sex], 0, axis.length, col = colours$Sex[sex], lwd = 4, center = centers[[time]])
    for(tissue in 1:n_tiss){
      polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Sex[sex], pch = 18, cex = 2.65*axis.length, center = centers[[time]])
    }
    for(tissue in 1:n_tiss){
      polarp(t = aas[[time]][sex], r = rs[tissue], col = "white", pch = 18, cex = 1.775*axis.length, center = centers[[time]])
    }
  }
  for(sex in 1:2){
    for(tissue in 1:n_tiss){
      polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Tissue[tissue], pch = 18, cex = 1.775*axis.length, center = centers[[time]])
      
    }
  }
  
  if(numbers_in_squares){
    for(sex in 1:2){
      for(tissue in 1:n_tiss){
        cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
        text(tissue, x = cartesian_coordinates[1], y = cartesian_coordinates[2], col = "white", 
             srt = c(45,-45,-45,45)[time], adj = 0.5, cex = 0.55, font = 2)   
      }
    }
  }
  
  
}
} else {
  for(time in 1:4){
    for(sex in 1:2){
      for(tissue in 1:n_tiss){
        cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
        if((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)){
          name_color = colours$Sex[sex]
        } else {
          name_color = colours$Sex[sex]
        }
        text(nice_names[tissue], x = cartesian_coordinates[1], y = cartesian_coordinates[2], 
             col = adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
                  !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)), 
             srt = c(-45,45,-45,45)[time], adj = 0.5, cex = 1 / strwidth(nice_names[tissue]) * tissue_name_cex, font = 2)   
      }
    }
  }
}
  
tissues <- names(nice_names)
if(across_relationships){
  for(time in 1:4){
    et_bt_across <- et_bt[[time]][et_bt[[time]][,"tiw"] == 3,]
    if(nrow(et_bt_across) == 0){
      next()
    }
    for(ri in 1:nrow(et_bt_across)){
      cartesian_coordinates_1 <- polar2cart(aas[[as.numeric(et_bt_across$time1[ri])]][as.numeric(et_bt_across$sex1[ri])], 
                                            rs[which(tissues == et_bt_across$tiss1[ri])]) + centers[[as.numeric(et_bt_across$time1[ri])]] * inner_shifter
      cartesian_coordinates_2 <- polar2cart(aas[[as.numeric(et_bt_across$time2[ri])]][as.numeric(et_bt_across$sex2[ri])], 
                                            rs[which(tissues == et_bt_across$tiss2[ri])]) + centers[[as.numeric(et_bt_across$time2[ri])]] * inner_shifter
      # cartesian_coordinates_1 <- (cartesian_coordinates_1 - centers[as.numeric(et_bt_across$time1[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time1[ri])][[1]]
      # cartesian_coordinates_2 <- (cartesian_coordinates_2 - centers[as.numeric(et_bt_across$time2[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time2[ri])][[1]]
      # cartesian_coordinates_1 <- cartesian_coordinates_1 * inner_shifter
      # cartesian_coordinates_2 <- cartesian_coordinates_2 * inner_shifter
      segments(x0 = cartesian_coordinates_1[1], y0 = cartesian_coordinates_1[2], x1 = cartesian_coordinates_2[1], y1 = cartesian_coordinates_2[2],
               lwd = et_bt_across$weight[ri]^line_weight_power * line_weight_multiplier,
               col = adjustcolor(et_bt_across$color[ri], et_bt_across$opacity[ri]))
    }
  }
}

xl <- -3.2; yb <- -0.425; xr <- -2.9; yt <- 0.575;
# et_wt[[time]]$weight[ri]^1 * 4
line_weights <- c(0:10/10)^line_weight_power * line_weight_multiplier
lws <- line_weights / 96 / (par("pin")[1]  / 4)
# segments(xl,yb,xl,yt,lwd = 10)
corresponding_heights <- seq(0,abs(yt - yb)/2,length.out = 11)


# line_weight_power <- 2
# line_weight_multiplier <- 10

hoffset <- 0 + c(5:0/5, 1:5/5)^line_weight_power * line_weight_multiplier / 300
voffset_rhos <- -0.1
text(labels = round(seq(-1, 1, length.out = 11), 2), y = seq(yb, yt, length.out = 11) + voffset_rhos, 
     x = xl - 0.0075 + hoffset, pos = 4, las=2, cex=0.9)
text(labels = latex2exp::TeX(paste0("$\\rho$")), y = yt - 0.03 + voffset_rhos, x = (xl), pos = 3, cex = 2)
text(labels = paste0("≥ ", corr_thresh), y = yt + 0.0325 + voffset_rhos, x = xl + 0.225, pos = 3, cex = 1.1, font = 3)
polygon(x = c(xl - lws/2, xl + rev(lws)/2),
        y = c((yb + yt) / 2 + corresponding_heights, yt - corresponding_heights) + voffset_rhos, col = "#e28a4a")
polygon(x = c(xl - lws/2, xl + rev(lws)/2),
        y = c((yb + yt) / 2 - corresponding_heights, yb + corresponding_heights) + voffset_rhos, col = "#2096be")

if(!tissue_names_not_colors){
points(x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), col = colours$Tissue, cex = 2, pch = 15)
if(numbers_in_squares){
    text(1:19, x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), col = "white", 
         adj = 0.5, cex = 0.55, font = 2)   
}
text(x = rep(xl, n_tiss), y = seq(yb - 0.25, yb - 0.25 - (yt - yb) / 11 * n_tiss, length.out = n_tiss), labels = nice_names, pos = 4, cex = 0.7)
}
# addImg(png::readPNG("~/Pictures/deathrats1.png"), 0, 0, width = 1)

dev.off()

#stitch pdfs together
pdftools::pdf_combine(paste0("~/Documents/figure2a_suggestion/network_visualization/", c("all-tissues", tissues), "_corrThresh-", corr_thresh, ".pdf"), 
                      output = paste0("~/Documents/figure2a_suggestion/hive-y_network-visualization.pdf"))



# for debugging
# ct <- 0.3
# fc <- which(abs(zcor) > ct, arr.ind = T)
# mc <- cbind(t(apply(as.matrix(strsplit(rownames(zcor)[fc[,"row"]], "_")), 1, unlist)), 
#             t(apply(as.matrix(strsplit(rownames(zcor)[fc[,"col"]], "_")), 1, unlist)), zcor[fc])
# mc <- as.data.frame(mc[,c(1,4,2,5,3,6,7)])
# colnames(mc) <- c("ts1", "ts2", "s1", "s2", "tm1", "tm2", "corr")
# 
# valid_pairs <- function(pair){
#   pair <- sort(pair)
#   vps <- cbind(c("1w", "2w", "4w", "8w"), c("2w", "4w", "8w", "1w"))
#   vps <- rbind(vps, cbind(vps[,2], vps[,1]))
#   return(any(paste0(pair, collapse = "") == sapply(1:nrow(vps), function(x) paste0(vps[x,], collapse = ""))))
# }
# # mc = mc[apply(mc[,5:6], 1, valid_pairs),]
# mc
# table(unlist(mc[,5:6])) / 2
# mc[intersect(grep(x = mc$ts1, "adre"), grep(x = mc$ts2, "ovar")),]
