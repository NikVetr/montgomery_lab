#### load packages ####
library(shiny)
library(shinyBS)
library(shinyWidgets)
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(MotrpacBicQC) # loading this package gives you access to several color palettes
library(pracma)

#### specify functions #### 
polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
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
arc <- function(t1,t2,r1,r2,res=50,lwd = 1,col=1, mindist = T, self_adjust = 2 * pi / 45, random_selfing = T, clockwise_selfing = T, pointy_selfing = F,
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

#### load data and colors ####
zcor = as.matrix(read.table("data/zcor.tsv"))
col_df = data.frame(row.names=rownames(zcor))
col_df$Tissue = gsub('_.*','',rownames(col_df))
col_df$Time = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[3])
col_df$Sex = sapply(rownames(col_df), function(x) unname(unlist(strsplit(x, '_')))[2])

colours = list(Tissue=tissue_cols[unique(col_df$Tissue)], 
               Time=group_cols[unique(col_df$Time)],
               Sex=sex_cols[c('male','female')])
colours$Tissue["t56-vastus-lateralis"] <- colorRampPalette(c(as.character(colours$Tissue["t56-vastus-lateralis"]), "black"))(3)[2]

nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])
nnmap <- nnmap[nnmap$tissue_name_release != "",]
tiss_ord <- nnmap$tissue_name_release[match(rev(MotrpacBicQC::tissue_order), nnmap$abbreviation)]
colours$Tissue <- colours$Tissue[match(tiss_ord, names(colours$Tissue))]
colours$Tissue <- colours$Tissue[!is.na(colours$Tissue)]
nice_names <- sapply(names(colours$Tissue), function(tissue) nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)])
tissues <- names(nice_names)

#### specify figure parameters ####
corr_thresh <- 0
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

#### ready pre-computation for plotting ####
# for(tissues_to_include in tissues){

tissues_to_include <- tissues#[1:5]
# print(tissues_to_include)

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
if(no_concentric_arcs){
  et_wt <- lapply(1:4, function(time) et_wt[[time]][et_wt[[time]]$tiss1 != et_wt[[time]]$tiss2,])
}

centers <- list(c(-1,1), c(1,1), c(1,-1), c(-1,-1))
centers <- lapply(centers, function(x) x * center_rescaler)

max_outer_radii <- sapply(et_wt, function(etwt) max(etwt$r1 + etwt$r2) / 2)
week_label_locs <- sapply(1:4, function(w) polar2cart(r = max_outer_radii[w], t = (-(0:3)*pi/2 + 3*pi/4)[w]) + centers[[w]])
week_label_locs[week_label_locs == Inf] <- -center_rescaler * 1.2
week_label_locs[week_label_locs == -Inf] <- center_rescaler * 1.2
timelab_nudge <- 0.25

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

et$tpairs <- sapply(1:nrow(et), function(ri) paste0(sort(c(et$time1[ri], et$time2[ri])), collapse = ""))
et$close_sex <- NA
et$close_sex[et$tpairs == "12"] <- 1
et$close_sex[et$tpairs == "23"] <- 2
et$close_sex[et$tpairs == "34"] <- 1
et$close_sex[et$tpairs == "14"] <- 2

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

xl <- -3.2; yb <- -0.425; xr <- -2.9; yt <- 0.575;
line_weights <- c(0:10/10)^line_weight_power * line_weight_multiplier
lws <- line_weights / 96 / (par("pin")[1]  / 4)
corresponding_heights <- seq(0,abs(yt - yb)/2,length.out = 11)

hoffset <- 0 + c(5:0/5, 1:5/5)^line_weight_power * line_weight_multiplier / 300
voffset_rhos <- -0.1

et_bt_full_allcorrs <- et_bt
et_wt_full_allcorrs <- et_wt

#### plot the thing already ####

time1 = Sys.time()
tissues_to_include <- tissues[3]
et_bt <- lapply(et_bt_full_allcorrs, function(et_bt_sub) et_bt_sub[et_bt_sub$tiss1 == tissues_to_include,])
et_wt <- lapply(et_wt_full_allcorrs, function(et_wt_sub) et_wt_sub[et_wt_sub$tiss1 == tissues_to_include,])

time2 = Sys.time()

par(mar = c(6,9,5,6), xpd = NA)
plot(1,1,xlim = c(-2,2), ylim = c(-2,2), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#plot tissue names

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
      
      cartesian_coordinates_mat <- t(sapply(1:n_tiss, function(tissue) polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]))
      name_color = colours$Sex[sex]
      tissue_cols_for_labels <- sapply(1:n_tiss, function(tissue) adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
                                                                                                   !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)))
      tissue_cex_for_labels <- sapply(1:n_tiss, function(tissue) 1 / strwidth(nice_names[tissue]) * tissue_name_cex)
      text(nice_names, x = cartesian_coordinates_mat[,1], y = cartesian_coordinates_mat[,2], 
           col = tissue_cols_for_labels, srt = c(-45,45,-45,45)[time], adj = 0.5, cex = tissue_cex_for_labels, font = 2)   
    }
    
    # for(tissue in 1:n_tiss){
    #   cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
    #   if((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)){
    #     name_color = colours$Sex[sex]
    #   } else {
    #     name_color = colours$Sex[sex]
    #   }
    #   text(nice_names[tissue], x = cartesian_coordinates[1], y = cartesian_coordinates[2], 
    #        col = adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
    #                                               !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)), 
    #        srt = c(-45,45,-45,45)[time], adj = 0.5, cex = 1 / strwidth(nice_names[tissue]) * tissue_name_cex, font = 2)   
    # }
    
  }
}

# my_plot <- recordPlot() 
time3 = Sys.time()
time3-time2 

#plot arcs within timepoint
for(time in 1:4){
  if(nrow(et_wt[[time]]) == 0){next()}
  for(ri in 1:nrow(et_wt[[time]])){
    arc(t1 = aas[[time]][et_wt[[time]][ri,"theta1"]], t2 = aas[[time]][et_wt[[time]][ri,"theta2"]], r1 = et_wt[[time]][ri,"r1"], center = centers[[time]] * outer_shifter,
        r2 = et_wt[[time]][ri,"r2"], lwd = et_wt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier, col = adjustcolor(et_wt[[time]]$color[ri], et_wt[[time]]$opacity[ri]), 
        random_selfing = F, clockwise_selfing = (et_wt[[time]][ri,"sex1"] == "1") * c(1,-1,1,-1)[time]*-1 + rev(c(0,1,0,1))[time], self_adjust = 0.3)
  }
}

for(i in 1:4){
  text(c("\u2642", "\u2640", "\u2642", "\u2640")[i], col = c(colours$Sex, colours$Sex)[i], 
       x = c(0,1,0,-1)[i]*(2 - sqrt(2) + axis.length)*1.15 - c(2,0,-2.15,0)[i] * 0.06, 
       y = c(1,0,-1,0)[i]*(2 - sqrt(2) + axis.length)*1.275 - c(-1,-1,-1,-1)[i] * 0.06, 
       cex = 4.5, srt = c(45,90,225,270)[i], pos = 1)
  
  #plot weeks
  shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = week_label_locs[1,i] + c(-1,1,1,-1)[i]*timelab_nudge, 
             y = week_label_locs[2,i] + c(1,1,-1,-1)[i]*timelab_nudge, cex = 4,
             srt = c(45,-45,235,135)[i], pos = 1)
  
}

time4 = Sys.time()
time4-time3

#plot adjacent relationships
if(T){
  
  for(time in 1:4){
    if(nrow(et_bt[[time]]) != 0){
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
  
}

time5 = Sys.time()


if(T){
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
      segments(x0 = cartesian_coordinates_1[1], y0 = cartesian_coordinates_1[2], x1 = cartesian_coordinates_2[1], y1 = cartesian_coordinates_2[2],
               lwd = et_bt_across$weight[ri]^line_weight_power * line_weight_multiplier,
               col = adjustcolor(et_bt_across$color[ri], et_bt_across$opacity[ri]))
    }
  }
}

time6 = Sys.time()

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

time7 = Sys.time()
cbind(paste0(1:6, "-", 2:7), diff(c(time1, time2, time3, time4, time5, time6, time7)))


#### create shiny app ####

ui <- fluidPage(
  
  tags$head(
    uiOutput("css")
  ), #cursor: url(https://emoji.discord.st/emojis/Reddit.png), auto;}
  
  # column(10, offset=0, titlePanel(div(img(src="motrpac-M.png", width = 50), "oTrPAC Transcriptome Correlation Network"))),
  
titlePanel(div(
              div(img(src="motrpac-M.png", width = 50), 
                style="position:absolute; top:0; left:0;margin-left:0;"), 
              div(("oTrPAC Transcriptome Correlation Network"), 
                style="position:absolute; top:0; left:0;margin-left:50px;"),
              style="position:relative; top:0; left:0;margin-left:10px;")
           ),

  # div(id="container", width = 625, height = 650, style="background: #FFFFFF;",
  #     div(plotOutput("name_background",
  #                    height = 650,
  #                    width = 625),
  #         style="position:absolute; top:0; left:0;margin-left:50px;"),
  #     div(uiOutput("plot.ui",
  #                  height = 650,
  #                  width = 625),
  #         style="position:absolute; top:0; left:0;margin-left:50px;")
  # )
  
  sidebarLayout(
    sidebarPanel(width = 3, 
                 # verbatimTextOutput("info"),
                 h5(strong("Arcs to Plot:")),
                 checkboxInput("adjacent", "Adjacent Timepoints", value = TRUE),
                 checkboxInput("nonadjacent", "Opposing Timepoints", value = TRUE),
                 checkboxInput("within_tiss", "Within Tissues", value = TRUE),
                 checkboxInput("between_tiss", "Between Tissues", value = TRUE),
                 sliderInput("positive_corr", "Positive Correlations:",
                             min = 0, max = 1,
                             value = c(0.3,1)),
                 sliderInput("negative_corr", "Negative Correlations:",
                             min = -1, max = 0,
                             value = c(-1,-0.3)),
                 img(src="MoTrPAC_blue-figures.png", width = 200), br(), br(),
                 p("This a quick demonstration of how a figure might be made interactive. 
                   My implementation likely leaves much to be desired with respect to speed, but I didn't want to optimize too many things prematurely.
                   Many more graphical parameters can be included in this right panel for the user to manipulate, but it seemed unwise to go too crazy.")
                 #,
                 # br(), br(),
                 # p("Hovering over the top right square -- at (1,1) -- changes the color of all squares to green. 
                 #   How can I change the cursor to, say, a pointing hand while that happens?")
                 
                 
    ), position = "right",
    
    # mainPanel(
    #   # imageOutput("bayes"),
    #   plotOutput("name_background", width = "625px", height = "650px"),
    #   uiOutput("plot.ui"),
    #   tableOutput("pointlist")
    # )
  
    mainPanel(
      fluidRow(
        div(id="container", width = 625, height = 650, style="background: #FFFFFF;",
            div(plotOutput("name_background",
                           height = 650,
                           width = 625),
                style="position:absolute; top:0; left:0;margin-left:50px;"),
            conditionalPanel(condition = "output.middleNames",
              div(plotOutput("name_background_solid_names",
                           height = 650,
                           width = 625),
                style="position:absolute; top:0; left:0;margin-left:50px;")),
            div(uiOutput("plot.ui",
                           height = 650,
                           width = 625),
                style="position:absolute; top:0; left:0;margin-left:50px;")
        )
      )
    )
  ))

server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2) # set maximum image size
  
  xy <- reactiveValues(x= numeric(0), y = numeric(0), col = numeric(0), line=numeric(0)) # add new points
  tissue_arcs <- reactiveValues(tissues = rep(T, n_tiss)) # add new points
  hide_solid_names <- reactiveValues(yes = F) # add new points
  
  et_bt_full <- reactive({
    lapply(et_bt_full_allcorrs, function(sub) sub[
                                              ((sub$corr > input$positive_corr[1] & sub$corr < input$positive_corr[2]) |
                                              (sub$corr > input$negative_corr[2] & sub$corr < input$negative_corr[1])) &
                                                ((sub$tiss1 == sub$tiss2) * input$within_tiss | 
                                                 (sub$tiss1 != sub$tiss2) * input$between_tiss),])
  })
  
  et_wt_full <- reactive({
    lapply(et_wt_full_allcorrs, function(sub) sub[
                                              ((sub$corr > input$positive_corr[1] & sub$corr < input$positive_corr[2]) |
                                               (sub$corr > input$negative_corr[2] & sub$corr < input$negative_corr[1])) &
                                                ((sub$tiss1 == sub$tiss2) * input$within_tiss | 
                                                 (sub$tiss1 != sub$tiss2) * input$between_tiss),])
  })
  
  output$css <- renderUI({
    
    if(any(!tissue_arcs$tissues)){
      css_string <- "
        #mainplot {  
        cursor: pointer;}"
    } else {
      css_string <- "
      #mainplot {  
      cursor: default;}"
    }
    tags$style(HTML(css_string))
  })
  
  output$plot.ui <- renderUI({
    plotOutput("mainplot", width = "625px", height = "650px",
               click = "plot_click",
               hover = hoverOpts(
                 "plot_hover",
                 delay = 10,
                 delayType = c("debounce", "throttle"),
                 clip = TRUE,
                 nullOutside = TRUE
               ))
  })
  
  
  # observe hover values and clicks
  observe({
    
    if (is.null(input$plot_hover)){
      return()
    }
    
    #coerce to quadrant one
    q1_coords <- abs(c(input$plot_hover$x, input$plot_hover$y))
    q1_coords[1] <- -q1_coords[1]
    q1_tissue_locs <- rbind(t(sapply(1:n_tiss, function(tissue) polar2cart(aas[[1]][1], rs[tissue]) + centers[[1]])), 
                            t(sapply(1:n_tiss, function(tissue) polar2cart(aas[[1]][2], rs[tissue]) + centers[[1]])))
    q1_dists <- sqrt(apply((matrix(q1_coords, ncol = 2, nrow = n_tiss*2, byrow = T) - q1_tissue_locs)^2, 1, sum))
    
    if(any(q1_dists < 0.125)){
      tissue_to_display <- (which.min(q1_dists)-1) %% n_tiss + 1
      tissue_arcs$tissues <- c(rep(F, tissue_to_display-1), T, rep(F, n_tiss - tissue_to_display))
      hide_solid_names$yes <- T
    } else {
      tissue_arcs$tissues <- rep(T, n_tiss)
      hide_solid_names$yes <- F
    }
    
  })
  
  output$middleNames <- reactive({
    hide_solid_names$yes == F
  })
  outputOptions(output, "middleNames", suspendWhenHidden = FALSE)
  
  
  ##here's a way to incorporate clicks
  # observe({
  #   
  #   if (is.null(input$plot_click)){
  #     return()
  #   }
  #   
  #   isolate({
  #     dists <- sqrt((xy$x - input$plot_click$x)^2 + (xy$y - input$plot_click$y)^2)
  #     pts_to_remove <- which(dists < 0.1)
  #     if(length(pts_to_remove) > 0){
  #       xy$x <- xy$x[-pts_to_remove]
  #       xy$y <- xy$y[-pts_to_remove]
  #       xy$col <- xy$col[-pts_to_remove]
  #     } else{
  #       xy$x <- c(xy$x, input$plot_click$x)
  #       xy$y <- c(xy$y, input$plot_click$y)
  #       xy$col <- c(xy$col, sample(1:10, 1))
  #     }
  #     
  #   })
  # })
  
  # first layer of the figure, the translucent names and other permanent fixtures
  output$name_background <- renderPlot({
    par(xpd = NA, bg = "transparent")
    tissues_to_include <- tissues
    plot(100,100,xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), col = "white", bg = "transparent",
         xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
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
          
          cartesian_coordinates_mat <- t(sapply(1:n_tiss, function(tissue) polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]))
          name_color = colours$Sex[sex]
          tissue_cols_for_labels <- sapply(1:n_tiss, function(tissue) adjustcolor(name_color, 0.5))
          tissue_cex_for_labels <- sapply(1:n_tiss, function(tissue) 1 / strwidth(nice_names[tissue]) * tissue_name_cex)
          text(nice_names, x = cartesian_coordinates_mat[,1], y = cartesian_coordinates_mat[,2], 
               col = tissue_cols_for_labels, srt = c(-45,45,-45,45)[time], adj = 0.5, cex = tissue_cex_for_labels, font = 2)   
        }
        
      }
    }
    
    text(labels = round(seq(-1, 1, length.out = 11), 2), y = seq(yb, yt, length.out = 11) + voffset_rhos, 
         x = xl - 0.0075 + hoffset, pos = 4, las=2, cex=0.9)
    text(labels = latex2exp::TeX(paste0("$\\rho$")), y = yt - 0.03 + voffset_rhos, x = (xl), pos = 3, cex = 2)
    # text(labels = paste0("≥ ", corr_thresh), y = yt + 0.0325 + voffset_rhos, x = xl + 0.225, pos = 3, cex = 1.1, font = 3)
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
    
    for(i in 1:4){
      text(c("\u2642", "\u2640", "\u2642", "\u2640")[i], col = c(colours$Sex, colours$Sex)[i], 
           x = c(0,1,0,-1)[i]*(2 - sqrt(2) + axis.length)*1.15 - c(0.25,0,-0.25,0)[i] * 0.06, 
           y = c(1,0,-1,0)[i]*(2 - sqrt(2) + axis.length)*1.275 - c(0,-1.25,-2,-1.25)[i] * 0.06, 
           cex = 4.5, srt = c(25,90,205,270)[i], pos = 1)
    }
    
  }, execOnResize = T, bg = "transparent")
  
  
  # second layer of the figure, the solid names when all arcs are plotted
  output$name_background_solid_names <- renderPlot({
    par(xpd = NA, bg = "transparent")
    tissues_to_include <- tissues
    plot(100,100,xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), col = "white", bg = "transparent",
         xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
  
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
          
          cartesian_coordinates_mat <- t(sapply(1:n_tiss, function(tissue) polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]))
          name_color = colours$Sex[sex]
          tissue_cols_for_labels <- sapply(1:n_tiss, function(tissue) adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
                                                                                                       !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)))
          tissue_cex_for_labels <- sapply(1:n_tiss, function(tissue) 1 / strwidth(nice_names[tissue]) * tissue_name_cex)
          text(nice_names, x = cartesian_coordinates_mat[,1], y = cartesian_coordinates_mat[,2], 
               col = tissue_cols_for_labels, srt = c(-45,45,-45,45)[time], adj = 0.5, cex = tissue_cex_for_labels, font = 2)   
        }
        
      }
    }
  
  }, execOnResize = T, bg = "transparent")
  
  
  #third layer of the figure, the arcs themselves
  output$mainplot <- renderPlot({
    par(xpd = NA, bg = "transparent")
    
    tissues_to_include <- tissues[tissue_arcs$tissues]
    
    if(length(tissues_to_include) != 1){
      et_bt <- et_bt_full()
      et_wt <- et_wt_full()
    } else {
      et_bt <- lapply(et_bt_full(), function(et_bt_sub) et_bt_sub[et_bt_sub$tiss1 == tissues_to_include | et_bt_sub$tiss2 == tissues_to_include,])
      et_wt <- lapply(et_wt_full(), function(et_wt_sub) et_wt_sub[et_wt_sub$tiss1 == tissues_to_include | et_wt_sub$tiss2 == tissues_to_include,])
    }
    
    # print(tissues_to_include)
    
    plot(100,100,xlim = c(-2.5,2.5), ylim = c(-2.5,2.5), col = "white", bg = "transparent",
         xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    for(time in 1:4){
      if(nrow(et_wt[[time]]) == 0){next()}
      for(ri in 1:nrow(et_wt[[time]])){
        arc(t1 = aas[[time]][et_wt[[time]][ri,"theta1"]], t2 = aas[[time]][et_wt[[time]][ri,"theta2"]], r1 = et_wt[[time]][ri,"r1"], center = centers[[time]] * outer_shifter,
            r2 = et_wt[[time]][ri,"r2"], lwd = et_wt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier, col = adjustcolor(et_wt[[time]]$color[ri], et_wt[[time]]$opacity[ri]), 
            random_selfing = F, clockwise_selfing = (et_wt[[time]][ri,"sex1"] == "1") * c(1,-1,1,-1)[time]*-1 + rev(c(0,1,0,1))[time], self_adjust = 0.3)
      }
    }

    
    for(i in 1:4){
      #plot weeks
      shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = week_label_locs[1,i] + c(-1,1,1,-1)[i]*timelab_nudge, 
                 y = week_label_locs[2,i] + c(1,1,-1,-1)[i]*timelab_nudge, cex = 4,
                 srt = c(45,-45,235,135)[i], pos = 1)
      
    }
    
    if(input$adjacent){
      
      for(time in 1:4){
        if(nrow(et_bt[[time]]) != 0){
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
      
    }
    
    if(input$nonadjacent){
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
          segments(x0 = cartesian_coordinates_1[1], y0 = cartesian_coordinates_1[2], x1 = cartesian_coordinates_2[1], y1 = cartesian_coordinates_2[2],
                   lwd = et_bt_across$weight[ri]^line_weight_power * line_weight_multiplier,
                   col = adjustcolor(et_bt_across$color[ri], et_bt_across$opacity[ri]))
        }
      }
    }
    
    #focal tissue name
    if(hide_solid_names$yes){
      for(time in 1:4){
        for(sex in 1:2){
          
          cartesian_coordinates_mat <- sapply(which(tissue_arcs$tissues), function(tissue) polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]])
          name_color = colours$Sex[sex]
          tissue_cols_for_labels <- sapply(which(tissue_arcs$tissues), function(tissue)
            adjustcolor(name_color, ifelse(!((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0)))
          tissue_cex_for_labels <- sapply(which(tissue_arcs$tissues), function(tissue) 1 / strwidth(nice_names[tissue]) * tissue_name_cex)
          text(nice_names[which(tissue_arcs$tissues)], x = cartesian_coordinates_mat[1], y = cartesian_coordinates_mat[2], 
               col = tissue_cols_for_labels, srt = c(-45,45,-45,45)[time], adj = 0.5, cex = tissue_cex_for_labels, font = 2)   
        }
      }
    }
    

    # points(xy$x, xy$y, pch = 19, cex = 2, col = xy$col)
    # points(x = c(-1,-1,1,1), y = c(-1,1,-1,1), cex = 3, col = square_color$col, pch = 15)

  }, bg = "transparent", execOnResize = T)
  
  ##some debugging stuff
  # output$pointlist <- renderTable(data.frame(x = xy$x, y = xy$y, col = xy$col))
  # output$info <- renderText({
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("NULL\n")
  #     paste0("x=", round(e$x, 2), " y=", round(e$y, 2), "\n")
  #   }
  #   
  #   paste0(
  #     "click: ", xy_str(input$plot_click),
  #     "hover: ", xy_str(input$plot_hover)
  #   )
  # })
  
}

shinyApp(ui, server)
