#define functions
find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 0.95, srt_height_rescaling_ratio = 1.4){
  
  strwidths <- strwidth(txt)
  strheight <- strheight(txt[1]) * srt_height_rescaling_ratio 
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  # compute_space_remaining(data, par)
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
                               str_width_lefter_start_ratio = 0.01, rect_rescaling_ratio = 0.95, col = "black"){
  
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(txt) * optimal_word_placement_inf$cex
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      text(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
           labels = word_to_write, pos = 4, col = col)
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

#specify data
txt <- sample(replace = T, size = 80, x = c("Akr1c19" , "Aldh1a2" , "Aqp1" , "Arhgef25" , "B4galt4" , "C1qtnf6" , "C1qtnf6" , "Dtx3" , 
         "Dtx3" , "Fbln2" , "Fbn1" , "Gatm" , "Gls" , "Gprc5b" , "LOC100910255" , "Lox" , "Lox" , 
         "Ltbp3" , "Ndrg4" , "Pdlim7" , "Sparc" , "Sparc"))
rect_coords <- list(x0=0.275, x1=1.5195, y0 = 2.751, y1 = 0.11547)

#perform analysis
optimal_word_placement_inf <- find_optimal_cex_and_lines(txt = txt, rect_coords = rect_coords)

#check plot
plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
rect(xleft = rect_coords$x0, ybottom = rect_coords$y0, xright = rect_coords$x1, ytop = rect_coords$y1, lwd = 1, col = "black")
text_wrapped_words(txt, rect_coords, optimal_word_placement_inf, justified = T, col = "white")

#individually colored letters
text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = txt_cols[char_i], cex = cex, ...)
  }
}

text_cols_bad <- function(string, cols, x, cex = 1, ...){
  for(char in 1:nchar(string)){
    text(labels = substr(string, char, char), col = cols[char], cex = cex, x = x + strwidth(substr(string, 1, char)) * cex, ...)
  }
}

plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
txt <- "abcdefg"
txt_cols <- 1:nchar(txt)
text_cols_bad(x = 0.5, y = 2, string = txt, cols = txt_cols, cex = 2, pos = 4)
text(x = 0.5, y = 2, srt = 45, labels = "text_cols_bad()", pos = 3)
text_cols(x = 0.5, y = 1, string = txt, cols = txt_cols, cex = 2, pos = 4)
text(x = 0.5, y = 1.1, srt = 45, labels = "text_cols()", pos = 3)
