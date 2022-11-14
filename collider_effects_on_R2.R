#functions
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}
gradleg <- function(loc, cols, labels, pos = 4, main = "", rasterize = F, border_in = NA, border_out = T, direc = "v", ...){
  
  #get metadata
  n <- length(cols)
  
  #get rect positions
  if(direc == "v"){
    locs <- data.frame(x1 = rep(loc[1], n),
                  x2 = rep(loc[2], n),
                  y1 = seq(loc[3], loc[4], length.out = n+1)[-(n+1)],
                  y2 = seq(loc[3], loc[4], length.out = n+1)[-1]
    )  
  } else {
    locs <- data.frame(x1 = seq(loc[1], loc[2], length.out = n+1)[-(length(n)+1)],
                  x2 = seq(loc[1], loc[2], length.out = n+1)[-1],
                  y1 = rep(loc[3], n),
                  y2 = rep(loc[4], n)
    )
  }
  
  #draw rects
  rect(locs$x1, locs$y1, locs$x2, locs$y2,
       col = cols, border = border_in)
  if(border_out){
    rect(loc[1], loc[3], loc[2], loc[4])
  }
  
  #draw text
  if(direc == "v") text(x = ifelse2(pos == 4, locs$x2, locs$x1), y = (locs$y1 + locs$y2) / 2, pos = pos, labels = labels)
  if(direc == "h") text(y = ifelse2(pos == 1, locs$x1, locs$x2), x = (locs$x1 + locs$x2) / 2, pos = pos, labels = labels)
  
  #draw title
  text(x = mean(loc[1:2]), y = loc[4], labels = main, pos = 3)
  
}

diag2 <- function(x){
  if(length(x) > 1){
    return(diag(x))
  } else {
    if(abs(x%%1)>1E-6){
      return(as.matrix(x))
    } else{
      return(diag(x))
    }
  }
}

xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

R2_from_Cov <- function(p, r_yx, r_xx_prop, long_way = F){
  
  #adjust supplied r_ux if impossible
  if(p == 0){return(0)}
  r_xx <- r_yx^2 + r_xx_prop * (1 - r_yx^2)
  r_ux <- sqrt(r_xx - r_yx^2)
  
  if(long_way){
    #generate covariance matrix
    k <- p + 2
    inds <- 3:k
    R <- diag(k)
    R[1,inds] <- R[inds,1] <- r_yx
    R[2,inds] <- R[inds,2] <- r_ux 
    R[inds,inds] <- diag(p) * (1 - r_xx) + r_xx
    
    A <- R[-inds, -inds]
    B <- R[-inds, inds]
    C <- R[inds, inds]
    Ci <- solve(C) #slow method for general inverses, vs rank 1 update
    
  } else {
    
    A <- diag(2)
    B <- rbind(rep(r_yx, p), rep(r_ux, p))
    Ci <- fastmatrix::sherman.morrison(a = diag2(rep(1 / (1 - r_xx), p)),
                                       b = rep(r_xx, p),
                                       d = rep(1, p), inverted = T)
  }
  
  #solve schur complement
  cond_R <- A - B %*% Ci %*% t(B)
  
  #return R2
  list(R2 = 1 - cond_R[1,1], 
       pr_yu = cond_R[1,2] / sqrt(cond_R[1,1]) / sqrt(cond_R[2,2]),
       semipr_yu = cond_R[1,2] / sqrt(cond_R[1,1]),
       cov_yu = cond_R[1,2])
  
}

# np <- sort(unique(c(0, round(sqrt(2)^(0:20)))))
np <- (2)^(0:10)
np <- setNames(np, np)
r_yxs <- c(0.01, 1:19/20, 0.99)
r_yxs <- setNames(r_yxs, r_yxs)
r_xx_props <- c(0:99/100, 0.999)
r_xx_props <- setNames(r_xx_props, r_xx_props)

out <- parallel::mclapply(r_xx_props, function(xrp){
  abind::abind(lapply(r_yxs, function(yr){
    do.call(rbind, lapply(np, function(ni){
      R2_from_Cov(ni, yr, xrp)
    }))
  }), along = 3)
}, mc.cores = 12)
out <- abind::abind(out, along = 4)
R2s <- out[,"R2",,]
partial_yus <- out[,"pr_yu",,] 
semipartial_yus <- out[,"semipr_yu",,] 
cov_yus <- out[,"cov_yu",,] 

#### plotting just one example ####
par(mfrow = c(1,3), xpd = NA)
cols <- viridis::viridis(length(np))

#plot 1 -- lines
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
     ylim = c(0,1), xlim = c(0,1))

#R2 lines
r_yx <- 0.15
near_r_yx <- order(abs(r_yxs - r_yx))[1:2]
prop_betw <- abs(r_yx - r_yxs[near_r_yx[1]]) / diff(r_yxs[near_r_yx])

for(ni in seq_along(np)){
  interp_R2 <- R2s[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
    R2s[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
  lines(r_xx_props, interp_R2, col = cols[ni])
}

#axes
segments(c(0,0),c(0,0),c(0,1),c(1,0), lwd = 2)
segments(0:5/5,0,0:5/5,-0.015, lwd = 2)
segments(0,0:5/5,-0.015*xyrat(),0:5/5, lwd = 2)
text(0:5/5, -0.015, 0:5/5, pos = 1)
text(-0.015*xyrat(), 0:5/5, 0:5/5, pos = 2)
text(0.5, -0.1, labels = "Proportion Maximum Cov(U, X)", pos = 1)
text(-0.2, 0.5, labels = latex2exp::TeX("Maximum Possible $R^2$ for Y ~ X"), pos = 1, srt = 90)
text(0.5, y = 1.05, paste0("Cov(Y, X) = ", round(r_yx, 2)))

#legend
gradleg(loc = c(1.025, 1.075, 0.5, 1), cols = cols, labels = np, main = latex2exp::TeX("\\textit{n}"), cex = 0.9)

#plot 2 -- residual correlation 
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
     ylim = c(-1,0), xlim = c(0,1))
for(ni in seq_along(np)){
  interp_pryu <- partial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
    partial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
  # interp_pryu <- semipartial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
  #   semipartial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
  # interp_pryu <- cov_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
  #   cov_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
  
  lines(r_xx_props, interp_pryu, col = cols[ni])
}

#axes
segments(x0 = c(0,0), y0 = c(0,0), x1 = c(1,0), y1 = c(0,-1), lwd = 2)
segments(0:5/5,0,0:5/5,0.015, lwd = 2)
segments(0,0:-5/5,-0.015*xyrat(),0:-5/5, lwd = 2)
text(0:5/5, 0.015, 0:5/5, pos = 3)
text(-0.015*xyrat(), 0:-5/5, 0:-5/5, pos = 2)
text(0.5, 0.15, labels = "Proportion Maximum Cov(U, X)", pos = 1)
text(-0.2, -0.5, labels = latex2exp::TeX("$r_{y,u}$ | x"), pos = 1, srt = 90)

#plot 3 -- graph
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
     ylim = c(0,2.5), xlim = c(0,5))
locs <- rbind(c(2,2),
             c(4,2),
             c(1,1),
             c(2,1),
             c(3,1),
             c(5,1),
             c(1,0.5),
             c(2,0.5),
             c(3,0.5),
             c(5,0.5)
)
node_cols <- c("blue", "red", rep("purple", 4), rep("orange", 4))

edges <- rbind(c(1,3),
               c(1,4),
               c(1,5),
               c(1,6),
               c(2,3),
               c(2,4),
               c(2,5),
               c(2,6),
               c(7,3),
               c(8,4),
               c(9,5),
               c(10,6)
)

psl <- c(rep(0.075, nrow(edges) - 4), rep(0.15, 4))
hs <- rep(0.15, nrow(locs))
ws <- rep(0.5, nrow(locs))
node_labels <- c("Y", "U", paste0("$X_", c(1:3, "n"), "$"), paste0("$\\epsilon_", c(1:3, "n"), "$"))

for(i in 1:nrow(locs)){
  rrect(loc = locs[i,], w = ws[i], h = hs[i], border = 1, col = adjustcolor(node_cols[i],0.1),
        lwd = 2, hat_prop = 0)
  text(x = locs[i,1], locs[i,2], labels = latex2exp::TeX(node_labels[i]), font = 2)
}

for(i in 1:nrow(edges)){
  grad_arrow_curve(c(locs[edges[i,1],1], 
                     locs[edges[i,2],1], 
                     locs[edges[i,1],2] - (hs[edges[i,1]]/2 + 0.025) * ifelse(sign(diff(locs[edges[i,],2])) == -1, 1, -1), 
                     locs[edges[i,2],2] + (hs[edges[i,2]]/2 + 0.025) * ifelse(sign(diff(locs[edges[i,],2])) == -1, 1, -1)), 
                   col_alpha = 0.75, direc = "v",
                   prop_shaft_length = psl[i], cols = c(node_cols[edges[i,1]], node_cols[edges[i,2]]),
                   w = 0.2, outline_col = 1, outline_lwd = 1.25)
}

#extra stuff
text(4, 0.75, "...", cex = 3)
text(0.5 - strwidth("{") * 1.5, 1, "n predictors", pos = 2); text(0.5, 1, "{", cex = 2.5)

#### plotting lots simultaneously ####

#packages
library(foreach)
library(doParallel)
library(parallel)

#plotting params
fps <- 60
ns <- 6
thin <- 1
frame_indices <- 1:(fps*ns)
plotting_indices <- seq(1, max(frame_indices), by = thin)
yr_seq <- seq(0.01, 0.99, length.out = fps * ns)

locs <- rbind(c(2,2),
              c(4,2),
              c(1,1),
              c(2,1),
              c(3,1),
              c(5,1),
              c(1,0.5),
              c(2,0.5),
              c(3,0.5),
              c(5,0.5)
)
node_cols <- c("blue", "red", rep("purple", 4), rep("orange", 4))

edges <- rbind(c(1,3),
               c(1,4),
               c(1,5),
               c(1,6),
               c(2,3),
               c(2,4),
               c(2,5),
               c(2,6),
               c(7,3),
               c(8,4),
               c(9,5),
               c(10,6)
)

psl <- c(rep(0.075, nrow(edges) - 4), rep(0.15, 4))
hs <- rep(0.15, nrow(locs))
ws <- rep(0.5, nrow(locs))
node_labels <- c("Y", "U", paste0("$X_", c(1:3, "n"), "$"), paste0("$\\epsilon_", c(1:3, "n"), "$"))

#create directories
main_dir <- "~/Pictures/R2_from_graph/"
frames_dir <- paste0(main_dir, "frames/")
if(!dir.exists(main_dir)){dir.create(main_dir)}
if(!dir.exists(frames_dir)){dir.create(frames_dir)}
if(length(list.files(frames_dir)) > 0){file.remove(paste0(frames_dir, list.files(frames_dir)))}

#initiate cluster
use_multicore <- T
if(use_multicore){
  stopCluster(cl)
  remove(cl)
  if(!exists("cl")){
    cl <- makeCluster(12, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()
}

foreach(fi=1:length(plotting_indices), .packages = c("png")) %dopar% {
# for(fi in seq_along(plotting_indices)[2]){
  
  cat(paste0(fi, " "))
  
  #main change
  r_yx <- yr_seq[plotting_indices[fi]]
  
  #start plotting
  png(filename = paste0(frames_dir, paste0(rep(0, 5 - nchar(fi)), collapse = ""), fi, ".png"), 
      width = 2560, height = 800, type = "cairo", pointsize = 35)
  
  par(mfrow = c(1,3), xpd = NA, mar = c(4,4,3,2))
  cols <- viridis::viridis(length(np))
  
  #plot 1 -- lines for R2s
  plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
       ylim = c(0,1), xlim = c(0,1))
  
  #lines
  near_r_yx <- sort(order(abs(r_yxs - r_yx))[1:2])
  prop_betw <- abs(r_yx - r_yxs[near_r_yx[1]]) / diff(r_yxs[near_r_yx])
  
  for(ni in seq_along(np)){
    interp_R2 <- R2s[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
      R2s[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
    lines(r_xx_props, interp_R2, col = cols[ni], lwd = 2)
  }
  
  #axes
  segments(c(0,0),c(0,0),c(0,1),c(1,0), lwd = 3)
  segments(0:5/5,0,0:5/5,-0.015, lwd = 3)
  segments(0,0:5/5,-0.015*xyrat(),0:5/5, lwd = 3)
  text(0:5/5, -0.015, 0:5/5, pos = 1)
  text(-0.015*xyrat(), 0:5/5, 0:5/5, pos = 2)
  text(0.5, -0.1, labels = latex2exp::TeX("Proportion Maximum Remaining Cor($X_i$, $X_j$) from U"), pos = 1)
  text(-0.175, 0.5, labels = latex2exp::TeX("Maximum Possible $R^2$ for Y ~ X"), pos = 1, srt = 90)
  text(0.2, y = 1.05, latex2exp::TeX(paste0("Cor(Y, $X_i$) = ", round(r_yx, 2), ifelse(nchar(round(r_yx, 2)) == 3, "0", ""))), pos = 4,
       cex = 2)
  
  #legend
  gradleg(loc = c(1.025, 1.075, 0.5, 1), cols = cols, labels = np, main = latex2exp::TeX("\\textit{n}"), cex = 0.9)
  
  #plot 2 -- residual correlation 
  plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
       ylim = c(-1,0), xlim = c(0,1))
  for(ni in seq_along(np)){
    interp_pryu <- partial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[1]]),] * (1-prop_betw) +
      partial_yus[as.character(np[ni]),as.character(r_yxs[near_r_yx[2]]),] * prop_betw
    lines(r_xx_props, interp_pryu, col = cols[ni], lwd = 2)
  }
  
  #axes
  segments(x0 = c(0,0), y0 = c(0,0), x1 = c(1,0), y1 = c(0,-1), lwd = 3)
  segments(0:5/5,0,0:5/5,0.015, lwd = 3)
  segments(0,0:-5/5,-0.015*xyrat(),0:-5/5, lwd = 3)
  text(0:5/5, 0.015, 0:5/5, pos = 3)
  text(-0.015*xyrat(), 0:-5/5, 0:-5/5, pos = 2)
  text(0.5, 0.15, labels = latex2exp::TeX("Proportion Maximum Remaining Cor($X_i$, $X_j$) from U"), pos = 1)
  text(-0.175, -0.5, labels = latex2exp::TeX("$\\rho_{y,u|x}$"), pos = 1, srt = 90)
  
  #plot 3 -- graph
  plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
       ylim = c(0.25,2.25), xlim = c(0.5,5))
  
  for(i in 1:nrow(locs)){
    rrect(loc = locs[i,], w = ws[i], h = hs[i], border = 1, col = adjustcolor(node_cols[i],0.1),
          lwd = 2, hat_prop = 0)
    text(x = locs[i,1], locs[i,2], labels = latex2exp::TeX(node_labels[i]), font = 2)
  }
  
  for(i in 1:nrow(edges)){
    grad_arrow_curve(c(locs[edges[i,1],1], 
                       locs[edges[i,2],1], 
                       locs[edges[i,1],2] - (hs[edges[i,1]]/2 + 0.025) * ifelse(sign(diff(locs[edges[i,],2])) == -1, 1, -1), 
                       locs[edges[i,2],2] + (hs[edges[i,2]]/2 + 0.025) * ifelse(sign(diff(locs[edges[i,],2])) == -1, 1, -1)), 
                     col_alpha = 0.75, direc = "v",
                     prop_shaft_length = psl[i], cols = c(node_cols[edges[i,1]], node_cols[edges[i,2]]),
                     w = 0.2, outline_col = 1, outline_lwd = 2)
  }
  
  #extra stuff
  text(4, 0.75, "...", cex = 3)
  text(0.5 - strwidth("{") * 1.5, 1, "n predictors", pos = 2); text(0.5, 1, "{", cex = 2.5)
  
  dev.off()
  
}


#stich together animation
file.remove(paste0(main_dir, list.files(main_dir, pattern = "*.mp4")))
system(paste0("cd ", main_dir, "; ffmpeg -r ", 
              fps / thin," -f image2 -s 1000x500 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p temp.mp4"))

#reverse and append and loop
system(paste0("cd ", main_dir, "; ", 
              "ffmpeg -i temp.mp4 -vf reverse rev_temp.mp4; ",
              "touch input.txt;",
              "echo \"file temp.mp4\nfile rev_temp.mp4\" > input.txt;",
              "ffmpeg -f concat -i input.txt -codec copy 1t_temp.mp4; ",
              "ffmpeg -stream_loop 1 -i 1t_temp.mp4 -c copy final.mp4"))
