# MATRIX GREEN RAIN WITH NUCLEOTIDES

#### introductory stuff ####
set.seed(777)

#packages
library(Biostrings)
library(graphics)
library(foreach)
library(doParallel)
library(parallel)

#funtions
logistic <- function(x) 1 / (1+exp(-x))
inv_logistic <- function(x) log(x / (1-x))
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

#define window size
nchar_vert <- 30 #originally 25
pwidth <- 1920
pwidth_in <- pwidth / 72
pheight <- 1080
pheight_in <- pheight / 72

#define background string behavior
font_nucs <- "Source Code Pro"
font_nucs <- "Futura T Medium"
font_fword <- "Matrix"
letter_color_base <- "#00FF41"
letter_color <- colorRampPalette(c(letter_color_base, "white"))(10)[3]
bright_letter_col <- colorRampPalette(c(letter_color_base, "white"))(10)[7]
bright_letter_col <- "#d5fce7"
prob_bright_first_letter <- 0.5
prob_glitch <- 0.02

#FINAL WORD
fword <- "THE MONTGOMERY LAB"
fword <- "THE EXTRACELLULAR MATRIX"
fword <- "THE MATRIX: CORRELATIONS"
fword_letters <- unique(strsplit(fword, "")[[1]])

#figure out string properties
cwidth <- max(sapply(c("A", "C", "G", "T"), function(character) strwidth(s = character, family = font_nucs, units = "in")))
cheight <- strheight(s = "ACGT\nACGT", family = font_nucs, units = "in") - strheight(s = "ACGT", family = font_nucs, units = "in")
cex <- pheight_in / cheight / nchar_vert
n_horiz_strings <- floor(pwidth_in / cwidth / cex)
glitch_height_ratio <- strheight(s = fword, family = font_fword, units = "in") / 
  strheight(s = paste0("ACGT"), family = font_nucs, units = "in")
glitch_height_ratio <- 0.9 #hmm better to manually code maybe

#get background strings
triplets <- names(GENETIC_CODE)[GENETIC_CODE != "*"]
stop_codons <- names(GENETIC_CODE)[GENETIC_CODE == "*"]
length_poss <- 1:6
nrow_matrix <- ceiling((nchar_vert + 100) / 3)
ncol_matrix <- n_seqs <- n_horiz_strings + 10 #add in buffer for later ad hoc shenanigans
AAs_matrix <- rbind(rep("ATG",ncol_matrix), 
                    matrix(sample(triplets, size = nrow_matrix * ncol_matrix, replace = T), nrow = nrow_matrix, ncol = ncol_matrix), 
                    sample(stop_codons, ncol_matrix, replace = T))
nucs_matrix <- sapply(1:ncol(AAs_matrix), function(gene) strsplit(paste0(AAs_matrix[,gene], collapse = ""), "")[[1]])
nwindow_sec <- runif(n_seqs, 0.25, 0.5)

#hide labmember names in nucs_matrix
current_members <- c("Nathan Abell",
                     "Olivia de Goede",
                     "Marianne DeGorter",
                     "Tiffany Eulalio",
                     "Nicole Ersaro",
                     "Nicole Gay",
                     "Mike Gloudemans",
                     "Pagé Goddard",
                     "Emily Greenwald",
                     "Tanner Jensen",
                     "Stephen Montgomery",
                     "Daniel Nachun",
                     "Kameron Rodrigues",
                     "Jarod Rutledge",
                     "Kevin Smith",
                     "Nikki Teran",
                     "Rachel Ungar",
                     "Nikolai Gates Vetr"
)
current_members <- c(
  "Pearson Correlation", "Covariance Matrix", "Correlation Matrix", 
  "Correlation Coefficient", "Spearman Correlation", "Multivariate Normal", 
  "Principal Component Analysis", "Partial Correlation", "Mutual Information", 
  "Regression to the Mean", "Kendall Tau", "Eigensystem", "Factor Analysis", 
  "Linear Regression", "Canonical Correlation", "Autocorrelation", 
  "Cross-Correlation", "Distance Correlation", "Copula", 
  "Curse of Dimensionality"
)

members <- sapply(current_members, function(name) toupper(strsplit(name, " ")[[1]][1]))
nx_appearing <- 2 #how many times should each member's name appear?
members <- rep(members, nx_appearing)
nc_members <- sapply(members, function(name) nchar(name))
name_columns <- sample(1:n_horiz_strings, size = length(members), replace = F) #everyone gets two columns
name_topletter_rows <- sapply(nc_members, function(minbotrow) sample(1:(nchar_vert-minbotrow), 1))
for(ni in 1:length(members)){
  nucs_matrix[name_topletter_rows[ni]:(name_topletter_rows[ni] + nc_members[ni] - 1), name_columns[ni]] <- rev(strsplit(members[ni], "")[[1]])
}

#OR change to arbitrary
change_to_arbitrary <- F
if(change_to_arbitrary){
  arbitrary_word <- "XENOTRANSPLANTATION "
  codeword <- paste0(rev(strsplit(arbitrary_word, "")[[1]]), collapse = "")
  nucs_matrix <- matrix(rep(strsplit(codeword, "")[[1]], 
                 times = ceiling(nrow(nucs_matrix) * ncol(nucs_matrix) / nchar(codeword)))[1:(nrow(nucs_matrix) * ncol(nucs_matrix))], 
                 nrow = nrow(nucs_matrix), ncol = ncol(nucs_matrix))  
}

# AAs <- sample(length_poss, size = n_seqs, replace = T)
# AAs <- sapply(AAs, function(nAA) c("ATG", sample(triplets, size = nAA, replace = T), sample(stop_codons, 1)))
# nucseqs <- sapply(AAs, function(aaseq) paste0(aaseq, collapse = ""))



#to change size, vary x and y limits
# zoom <- 1

#### test basic figure ####
# png(filename = paste0("~/Documents/matrix_animation/frame.png"), width = pwidth, height = pheight, family = font_nucs)
# par(mfrow = c(1, 1), mar = c(0,0,0,0), xpd = NA)
# plot(1,1, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", xlim = c(1, n_horiz_strings) / zoom, ylim = c(1,nchar_vert) / zoom)
# rect(xleft = -1E5, xright = 1E5, ybottom = -1E5, ytop = 1E5, col = "black")
# 
# for(x in -2:(n_horiz_strings)){
#   cex_multiplier <- runif(1, 0.95, 1)
#   yloc <- sample(1:nchar_vert, 1)
#   stp_length <- sample((length_poss+2)*3,1) 
#   seqtp <- paste0(nucs_matrix[yloc:(yloc-1+stp_length), x + 3], collapse = "")
#   charcols <- (c(sample(c(letter_color_base, bright_letter_col), 1, prob = c(1-prob_bright_first_letter, prob_bright_first_letter)), 
#                  colorRampPalette(c(letter_color_base, "black"))(nchar(seqtp) + 2)[1:(nchar(seqtp)-1)]))
#   for(char in 1:(nchar(seqtp)-1)){
#     glitch_or_no <- as.logical(rbinom(1, 1, prob_glitch))
#     letter_to_plot <- ifelse(glitch_or_no, sample(fword_letters, 1), c(substr(seqtp, char, char)))
#     font_to_use <- ifelse(glitch_or_no, font_fword, font_nucs)
#     text(x, yloc + char, letter_to_plot, cex = cex * cex_multiplier * zoom, family = font_to_use, col = charcols[char], pos = 1, xpd = NA)  
#     
#   }
#   
#   
# }
# 
# dev.off()


#### animate falling rain ####

if(!exists("cl")){
  cl <- makeCluster(16, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

renderVideo <- T
file.remove(paste0("~/Documents/matrix_animation/frames/", list.files(path = "~/Documents/matrix_animation/frames/", pattern = ".png")))
if(!dir.exists("~/Documents/matrix_animation/")){
  dir.create("~/Documents/matrix_animation/")  
}
if(!dir.exists("~/Documents/matrix_animation/frames")){
  dir.create("~/Documents/matrix_animation/frames")  
}

#specify animation-specific parameters
framethin <- 1
nfps <- 60

#specify misc. other properties
glitch_duration_in_sec <- 0.25
zoom_start <- 1
zoom_max <- 3
zoom_max_2 <- 12

#timing of events in phase #1 of animation
#first number is start time, second is duration
# generic_rain <- c(0,6) * nfps #on second thought, let's not explicitly constrain this
rain_start_times <- runif(n_seqs, 0.1, 1)
show_stephen <- c(0.25,1.25) * nfps 
hide_stephen <- c(sum(show_stephen) / nfps + 0.25,3) * nfps 
zoom_in <- c(6, 4) * nfps
reveal_fword <- c((sum(zoom_in) - zoom_in[2]/2) / nfps, 3) * nfps
hide_fword <- c(sum(reveal_fword) / nfps + 2, 2) * nfps
zoom_in_again <- c(sum(reveal_fword) / nfps + 1.25, 5) * nfps
extra_buffer <- c(sum(hide_fword) / nfps, 0.5) * nfps
nsec <- sum(extra_buffer) / nfps
nf <- (nsec * nfps)

# pre-roll glitch behavior
glitch_letters <- array(data = NA, dim = c(nrow(nucs_matrix), ncol(nucs_matrix), nf))
n_glitches <- round(prob_glitch * length(glitch_letters) / (glitch_duration_in_sec * nfps))
glitch_starts <- sapply(c(nrow(nucs_matrix), ncol(nucs_matrix), nf), function(x) sample(x, n_glitches, replace = T))
glitch_vals <- sample(fword_letters, n_glitches, replace = T)

for(i in 1:n_glitches){
  glitch_letters[glitch_starts[i,1], glitch_starts[i,2], glitch_starts[i,3]:min((glitch_starts[i,3] + (glitch_duration_in_sec * nfps)), nf)] <- glitch_vals[i]
}

#pre-roll substring lengths
stp_lengths <- sample((length_poss+2)*3, n_seqs, replace = T) 

#pre-roll character first letter colors
first_lett_cols <- sample(c(letter_color_base, bright_letter_col), n_seqs, prob = c(1-prob_bright_first_letter, prob_bright_first_letter), replace = T)  

#specify f-word behavior
ncfw <- nchar(fword)
fword_xloc <- round(n_horiz_strings / 2 - ncfw / 2)
fword_yloc <- ceiling((nchar_vert+1) / 2)
fword_letter_appearances <- sample(reveal_fword[1]:sum(reveal_fword), ncfw)
fword_letter_disappearances <- sample(hide_fword[1]:sum(hide_fword), ncfw)
disapperance_fade_speed_in_frames <- 0.25 * nfps

#make "GO" the last substring to disappear
if(fword == "THE MONTGOMERY LAB"){
  GO_times <- fword_letter_disappearances[stringr::str_locate(string = fword, "GO")]
  last_time <- max(fword_letter_disappearances)
  fword_letter_disappearances[which.max(fword_letter_disappearances)] <- min(GO_times)
  nsec_diff_GO <- 0.25
  fword_letter_disappearances[stringr::str_locate(string = fword, "GO")] <- last_time + c(-1/4,1/4) * nfps * nsec_diff_GO
  fword_letter_disappearances[-stringr::str_locate(string = fword, "GO")] <- fword_letter_disappearances[-stringr::str_locate(string = fword, "GO")] - nsec_diff_GO * nfps
}

#make sure to zoom in on fword approx. center
fword_center <- c(fword_xloc + ncfw / 2 - 5.25, fword_yloc - 3)

#tv grid effect, add line segments?
segment_sep <- 0.05
segment_alpha <- 0.15
base_xlim <- c(1,n_horiz_strings)
base_ylim <- c(1,nchar_vert)
xlocs_tv <- seq(base_xlim[1], base_xlim[2], by = segment_sep * diff(base_xlim) / diff(base_ylim) / 2)
ylocs_tv <- seq(base_ylim[1], base_ylim[2], by = segment_sep)

grid_alpha_range <- c(0, 0.125)
screen_flicker_duration_in_sec <- 0.5
v_grid_alphas <- t(replicate(n = ceiling(nf / (screen_flicker_duration_in_sec * nfps)), 
                             runif(n = length(xlocs_tv) - 1, min = grid_alpha_range[1], grid_alpha_range[2])))
h_grid_alphas <- t(replicate(n = ceiling(nf / (screen_flicker_duration_in_sec * nfps)),
                             runif(n = length(ylocs_tv) - 1, min = grid_alpha_range[1], grid_alpha_range[2])))

#stephen's face
face <- png::readPNG("~/data/smontgom/matrix_stephen.png")
face <- png::readPNG("~/data/kate_matrix.png")
face <- png::readPNG("~/data/pig_matrix.png")
face <- png::readPNG("~/Downloads/file.png")

## now work backwards to make the 2nd wave of strings hit when the fword letters appear ##
fword_inds <- fword_xloc:(fword_xloc+ncfw-1)
#first stack the deck a bit
nwindow_sec[fword_inds] <- runif(ncfw, 
                                 quantile(nwindow_sec[-fword_inds], 0.35), quantile(nwindow_sec[-fword_inds], 0.95))
rain_start_times[fword_inds] <- runif(ncfw, 
                                      quantile(rain_start_times[-fword_inds], 0.35), quantile(rain_start_times[-fword_inds], 0.95))
stp_lengths[fword_inds] <- sample((3:5)*3, ncfw, replace = T) 

#now ensure put the thumb on the scale even more
fword_hit_frames <- sapply(fword_inds, function(x) which(sapply(1:1E4, function(f1) #warning: assumes it hits w/in 1E4 frames
  (nchar_vert + 1 - floor(nwindow_sec[x] * nchar_vert / nfps * f1) + floor(rain_start_times[x] * nfps) + 
     nchar_vert + stp_lengths[x]) == fword_yloc)
))
fword_hit_frames <- sapply(fword_hit_frames, function(hf) min(hf))
nwindow_sec[fword_inds] <- nwindow_sec[fword_inds] * fword_hit_frames / fword_letter_appearances

#now do the actual plotting
foreach(f1=seq(1, nf, framethin), .packages = c("png")) %dopar% {
  # for(f1 in seq(1, nf, framethin)[100:150]){
  logistic_speedup <- 8
  logistic_speedup_2 <- 20
  # zoom <- zoom_start + min(max((f1 - zoom_in[1]) / zoom_in[2] * (zoom_max - zoom_start), 0), zoom_max - zoom_start) #linear zoom
  # zoom <- zoom_start + logistic(((f1 - zoom_in[1] - zoom_in[2] / 2) / zoom_in[2]) * logistic_speedup) * (zoom_max - zoom_start) #sigmoid zoom
  zoom <- zoom_start + logistic(((f1 - zoom_in[1] - zoom_in[2] / 2) / zoom_in[2]) * logistic_speedup) * (zoom_max - zoom_start) + #double sigmoid zoom
    logistic(((f1 - zoom_in_again[1] - zoom_in_again[2] / 2) / zoom_in_again[2]) * logistic_speedup_2) * (zoom_max_2 - zoom_max)
  
  cat(paste0(f1, " "))
  png(filename = paste0("~/Documents/matrix_animation/frames/frame_", 
                        paste0(rep(0, 5-nchar(((f1 - 1) / framethin) + 1)), collapse = ""), ((f1 - 1) / framethin) + 1,".png"), 
      width = pwidth, height = pheight, family = font_nucs)
  par(mfrow = c(1, 1), mar = c(0,0,0,0), xpd = NA)
  
  #figure out plot limits
  xlim <- (c(0, diff(base_xlim)) - fword_center[1]) * 1 / zoom + fword_center[1] + 3
  ylim <- (c(0, diff(base_ylim)) - fword_center[2]) * 1 / zoom + fword_center[2] + 2
  
  #start plotting
  plot(1,1, col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
  
  #black (or other color) background
  rect(xleft = -1E5, xright = 1E5, ybottom = -1E5, ytop = 1E5, col = "black")
  
  #iterative over columns
  for(x in 1:(n_horiz_strings+4)){
    
    yloc <- nchar_vert + 1 - floor(nwindow_sec[x] * nchar_vert / nfps * f1) + floor(rain_start_times[x] * nfps)
    
    stp_length <- stp_lengths[x]
    
    if((yloc+stp_length-1) > 0){
      seqtp <- (paste0(nucs_matrix[max(1,yloc):(yloc+stp_length-1), x + 3], collapse = ""))
      if(yloc < 1){yloc <- 1}
    } else if((yloc+stp_length*2+nchar_vert) > 0){ #simulate second rainfall
      yloc <- yloc + nchar_vert + stp_length
      seqtp <- (paste0(nucs_matrix[max(1,yloc):(yloc+stp_length-1), x + 3], collapse = ""))
      if(yloc < 1){yloc <- 1}
    } else {
      next()
    }
    
    # charcols <- (c(first_lett_cols[x], colorRampPalette(c(letter_color_base, "black"))(nchar(seqtp))[1:(nchar(seqtp)-1)]))
    charcols <- (c(first_lett_cols[x], colorRampPalette(c(letter_color_base, "black"))(stp_length)[1:(stp_length-1)]))
    charcols <- tail(charcols, n = nchar(seqtp))
    
    #draw characters in column
    for(char in 1:(nchar(seqtp)-1)){
      
      if((yloc + char) < 1){
        next
      }
      
      glitch_or_no <- !is.na(glitch_letters[yloc + char, x, f1])
      letter_to_plot <- ifelse(glitch_or_no, glitch_letters[yloc + char, x, f1], c(substr(seqtp, char, char)))
      font_to_use <- ifelse(glitch_or_no, font_fword, font_nucs)
      cex_multiplier <- ifelse(glitch_or_no, 1/glitch_height_ratio, 1)
      y_disp <- ifelse(glitch_or_no, -0.05, 0)
      text(x, yloc + char + y_disp, letter_to_plot, cex = cex * cex_multiplier * zoom, family = font_to_use, col = charcols[char], pos = 1, xpd = NA)  
      
    }
    
    
  }
  
  #black out space underneath f-word if its letter has appeared, and also plot its letter and make it fade out when the time is right
  for(fc in 0:(ncfw-1)){
    if(f1 >= fword_letter_appearances[fc+1]){
      letter_to_plot <- substr(fword, fc+1, fc+1)
      rect(xleft = fword_xloc + fc - 0.51, xright = fword_xloc + fc + 0.5, ybottom = -1E5, ytop = fword_yloc + 0.325, col = "black")
      text(fword_xloc + fc, fword_yloc - 0.05, letter_to_plot, cex = cex * zoom / glitch_height_ratio, family = font_fword, col = bright_letter_col, pos = 1, xpd = NA)  
    }
    if(f1 >= fword_letter_disappearances[fc+1]){
      rect(xleft = fword_xloc + fc - 0.51, xright = fword_xloc + fc + 0.5, ybottom = -1E5, ytop = 1E5, 
           col = grDevices::adjustcolor("black", alpha.f = (f1 - fword_letter_disappearances[fc+1]) / disapperance_fade_speed_in_frames))
    }
  }
  
  #overlay a very low-alpha grid on everything? to mimic tv screen
  segments(x0 = xlocs_tv, x1 = xlocs_tv, y0 = -1E5, y1 = 1E5, lwd = zoom,
           col = grDevices::adjustcolor(1, segment_alpha))
  segments(x0 = -1E5, x1 = 1E5, y0 = ylocs_tv, y1 = ylocs_tv, lwd = zoom, 
           col = grDevices::adjustcolor(1, segment_alpha))
  
  #add little cells to mimic tv screen too
  v_alphas <- v_grid_alphas[floor((f1-1) / (screen_flicker_duration_in_sec * nfps)) + 1,]
  h_alphas <- h_grid_alphas[floor((f1-1) / (screen_flicker_duration_in_sec * nfps)) + 1,]
  for(ri in 1:length(h_alphas)){
    rect(xleft = -1E5, xright = 1E5, ybottom = ylocs_tv[ri], ytop = ylocs_tv[ri+1], 
         col = grDevices::adjustcolor(1, alpha = h_alphas[ri]), border = NA)
  }
  for(ci in 1:length(v_alphas)){
    rect(xleft = xlocs_tv[ci], xright = xlocs_tv[ci+1], ybottom = -1E5, ytop = 1E5, 
         col = grDevices::adjustcolor("black", alpha = v_alphas[ci]), border = NA)
  }
  
  #add and fade out stephen face
  if(f1 <= sum(hide_stephen)){
  face_to_plot <- face
  face_to_plot[,,4] <- face_to_plot[,,4] * (1-min(1,max(0,(f1-hide_stephen[1]) / hide_stephen[2]))) * min(1,max(0,(f1-show_stephen[1]) / show_stephen[2]))
  # addImg(face_to_plot, x = (n_horiz_strings) / 1.75, y = 11, width = 50)
  addImg(face_to_plot, x = (n_horiz_strings) / 1.9, y = 12, width = 45)
  }
  
  dev.off() 
  
}

stopCluster(cl = cl)
rm(cl)

#render video w/ filters
base_filename <- paste0("2DArray_animation_", pwidth, "x" , pheight)
raw_filename <- paste0("raw_", base_filename, ".mp4")
blur_filename <- paste0("blur_", base_filename, ".mp4")
final_filename <- paste0(base_filename, ".mp4")

if(file.exists(paste0("~/Documents/matrix_animation/", raw_filename))){file.remove(paste0("~/Documents/matrix_animation/", raw_filename))}
if(file.exists(paste0("~/Documents/matrix_animation/", blur_filename))){file.remove(paste0("~/Documents/matrix_animation/", blur_filename))}
if(file.exists(paste0("~/Documents/matrix_animation/", final_filename))){file.remove(paste0("~/Documents/matrix_animation/", final_filename))}

system(paste0("cd ~/Documents/matrix_animation; ffmpeg -r ", nfps / framethin," -f image2 -s ", pwidth, "x" , pheight," -i frames/frame_%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ", raw_filename))
system(paste0("cd ~/Documents/matrix_animation; ffmpeg -i ", raw_filename," -vf gblur=sigma=20:steps=6 -pix_fmt yuv420p ", blur_filename))
system(paste0("cd ~/Documents/matrix_animation; ffmpeg -i ", raw_filename," -i ", blur_filename, " -filter_complex \"blend=lighten\" ", final_filename))

