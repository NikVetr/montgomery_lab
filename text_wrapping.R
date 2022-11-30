#define functions
text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}

text3 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, ...){
  
  #convert text label to expression
  if(all(class(labels) == "character")){
    word_expression <- latex2exp::TeX(labels)  
  } else {
    word_expression <- labels
  }
  
  #find general params
  strw <- strwidth(word_expression, cex = cex)
  strh <- strheight(word_expression, cex = cex)
  base_strh <- strheight(latex2exp::TeX("GIs"), cex = cex)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strh / 2, 0)
  
  #adjust in case of ledding
  ebot <- strheight(latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", labels)), cex = cex) - strh
  etop <- strheight(latex2exp::TeX(gsub("\\^", "a", labels)), cex = cex) - strh
  ebottop <- strheight(latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u|\\^", "a", labels)), cex = cex) - strh
  
  #ugh this was obnoxious to figure out
  adj_ledding <- ifelse(abs(ebottop - (ebot + etop)) > 1E-6, 
                        ebot / 2 - (ebottop - ebot) / 2, 
                        ebot / 2 - etop / 2)
  adj_ledding <- adj_ledding - ifelse(base_strh > strh, (base_strh - strh) / 2, 0)
  
  #print the text itself
  text(x = adj_x, y = adj_y + adj_ledding, labels = word_expression, pos = NULL, cex = cex, col = col, ...)
  
  #draw a box around it if desired
  if(drect){
    rect(xleft = adj_x - strw / 2, 
         xright = adj_x + strw / 2, 
         ybottom = adj_y - strh / 2 + adj_ledding, 
         ytop = adj_y + strh / 2 + adj_ledding, border = col)
  }
}


text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text2(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = cols[char_i], cex = cex, ...)
  }
}

find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 1, str_height_rescaling_ratio = 1, fixed_cex = NA,
                                       newlines = NA){
  
  strwidths <- strwidth(latex2exp::TeX(txt))
  strheight <- (strheight("M\nM") - strheight("M")) * str_height_rescaling_ratio
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  if(all(is.na(newlines))){
    newlines <- rep(F, length(txt))
  }
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight,
               newlines = newlines)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.95)
  }
  # plot(1:120/100, sapply(log(1:120/100), function(cex) compute_space_remaining(data, cex)), type = "l")
  if(is.na(fixed_cex)){
    opt_cex <- suppressWarnings(optimx::optimx(par = par, fn = compute_space_remaining, data = data, hess = NULL, 
                                               method = c('nlminb', 'nlm', 'BFGS', 'Nelder-Mead')[4], hessian=FALSE, #can't compute hessian bc of sharp jumps when new line is formed? or maybe not?
                                               control = list(maxit = 1E4, trace = 0, kkt=FALSE)))
    
    #check if it worked
    compute_space_remaining(data = data, par = opt_cex$p1)
    compute_space_remaining(data = data, par = opt_cex$p1 + 0.0001)
    
    return(list(cex = exp(opt_cex$p1), 
                words_on_lines = put_words_on_lines(data = data, par = opt_cex$p1),
                space_width_min = space_width_min * exp(opt_cex$p1),
                vertical_space = strheight * exp(opt_cex$p1),
                vertical_space_noWS = strheight("M") * exp(opt_cex$p1))
    )
  } else {
    return(list(cex = fixed_cex, 
                words_on_lines = put_words_on_lines(data = data, par = fixed_cex),
                space_width_min = space_width_min * fixed_cex,
                vertical_space = strheight * fixed_cex,
                vertical_space_noWS = strheight("M") * fixed_cex)
    )
  }
  
  
}

compute_space_remaining <- function(data, par){
  put_words_on_lines(data, par, return_space_remaining = T)
}

put_words_on_lines <- function(data, par, return_space_remaining = F){
  
  #clarify par-cex relationship
  cex <- exp(par)
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  newlines <- data$newlines
  
  #check that no words are wider than a line
  if(any(strwidths_mod > rectwidth_mod) & return_space_remaining){
    return(Inf)
  }
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi] + space_width_min_mod
  txt_lines <- list()
  txt_lines[[1]] <- txi
  txt_lines[[2]] <- integer(0)
  linei <- linei + newlines[txi]
  
  while(txi < length(txt)){
    
    txi <- txi + 1
    prop_current_width <- current_width + strwidths_mod[txi]
    add_to_line <- (prop_current_width < rectwidth_mod) | (length(txt_lines[[linei]]) == 0)
    
    if(add_to_line){
      txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
      current_width <- prop_current_width + space_width_min_mod
    } else {
      current_width <- strwidths_mod[txi] + space_width_min_mod
      linei <- linei + 1
      txt_lines[[linei]] <- txi
    }
    
    if(newlines[txi]){
      linei <- linei + 1
      current_width <- 0
      txt_lines[[linei]] <- integer(0)
    }
    
    #print for debugging
    # print(paste0("line: ", linei, ", txt_i: ", txi, ", txt: ", txt[txi], ", nl: ", newlines[txi]))
    
  }
  
  if(return_space_remaining){
    last_line_width_remaining <- rectwidth_mod - sum(strwidths_mod[txt_lines[[linei]]])
    current_height <- linei * strheight_mod
    vspace_remaining <- (rectheight_mod - current_height) * rectwidth_mod
    hspace_remaining <- (last_line_width_remaining * strheight_mod)
    
    if(vspace_remaining < 0 | hspace_remaining < 0){
      return(Inf)
    } else {
      total_space_remaining <- vspace_remaining + hspace_remaining
      return(total_space_remaining)
    }
  } else {
    return(txt_lines)
  }
  
}

text_wrapped_words <- function(txt, rect_coords, optimal_word_placement_inf, justified = F, str_height_lower_start_ratio = 0, 
                               str_width_lefter_start_ratio = 0, rect_rescaling_ratio = 1, col = "black", multicolor_words = F, cols_list,
                               vertically_justified = F,...){
  
  ws_height <- optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS
  cex <- optimal_word_placement_inf$cex
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space_noWS * (1 + str_height_lower_start_ratio) / 2
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(latex2exp::TeX(txt), cex = cex)
  
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  if(vertically_justified){
    space_taken_vertically <- optimal_word_placement_inf$vertical_space_noWS * nlines + 
      (optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS) * (nlines - 1)
    space_left_vertically <- abs(rect_coords$y0 - rect_coords$y1) - space_taken_vertically
    extra_leading <- space_left_vertically / (nlines - 1)
  } else {
    extra_leading <- 0
  }
    
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      
      #adjust for sticky-outy bits
      word_expression <- latex2exp::TeX(word_to_write)
      
      
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = cex,
                  string = word_expression, pos = 4, cols = cols_list[[words_written]])
      } else {
        # text2(x = curr_x, y = curr_y, cex = cex,
        #       labels = word_expression, pos = 4, col = col, drect = T)
        text3(x = curr_x, y = curr_y, cex = cex,
             labels = word_to_write, pos = 4, col = col, drect = T)
        abline(h=curr_y - strheight(latex2exp::TeX("GIs"), cex = cex) / 2, lwd = 0.5)
      }
      if(justified & (linei != nlines)){
        curr_x <- curr_x + strwidth(word_expression) * cex + justified_space_between_words[linei]
      } else {
        curr_x <- curr_x + strwidth(word_expression) * cex + optimal_word_placement_inf$space_width_min  
      }
      
    }
    curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
    curr_y <- curr_y - optimal_word_placement_inf$vertical_space * (1 + str_height_lower_start_ratio) - extra_leading
  }
  
}

string_to_tokens <- function(txt_string){
  
  txt <- strsplit(txt_string, split = " ")[[1]]
  
  #find relevant LaTeX notation & modify tokens accordingly
  math_pairs <- do.call(rbind, lapply(seq_along(txt), function(i){
    mathi = which(strsplit(txt[i], "")[[1]] == "$")
    if(length(mathi) != 0){
      return(data.frame(token = i, math_i = mathi))
    } else {
      return(integer(0))
    }
  }))
  math_pairs <- data.frame(cbind(matrix(math_pairs$token, ncol = 2, byrow = T), 
                                 matrix(math_pairs$math_i, ncol = 2, byrow = T)))
  colnames(math_pairs) <- c("open_i", "close_i", "char_i_open", "char_i_close")
  for(i in 1:nrow(math_pairs)){
    if(math_pairs$open_i[i] == math_pairs$close_i[i]){
      next()
    } else {
      txt[math_pairs$open_i[i]] <- paste0(txt[math_pairs$open_i[i]:math_pairs$close_i[i]], collapse = " ")
      txt <- txt[-((math_pairs$open_i[i]+1):math_pairs$close_i[i])]
    }
  }
  
  bracket_pairs <- list(
    do.call(rbind, lapply(seq_along(txt), function(i){
      openi <- which(strsplit(txt[i], "")[[1]] == "{")
      if(length(openi) != 0){
        slashi <- which(sapply(strsplit(txt[i], "")[[1]], grepl, pattern = '\t'))
        if(length(slashi) == 0){
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = "",
                            style = ""))
        } else {
          styles <- sapply(seq_along(slashi), function(ci) substr(txt[i], slashi[ci], openi[ci]))
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = slashi,
                            style = styles))  
        }
      } else{
        return(integer(0))
      }
    })
    ),
    
    do.call(rbind, lapply(seq_along(txt), function(i){
      closei <- which(strsplit(txt[i], "")[[1]] == "}")
      if(length(closei) != 0){
        return(data.frame(token_close = i, 
                          close_i = closei))
      } else{
        return(integer(0))
      }
    })
    )
  )
  
  bracket_pos <- rbind(data.frame(type = "o", 
             pos = bracket_pairs[[1]]$token_open + 
               bracket_pairs[[1]]$open_i / 
               nchar(txt[bracket_pairs[[1]]$token_open]),
             ind = 1:nrow(bracket_pairs[[1]])),
        data.frame(type = "c", 
                   pos = bracket_pairs[[2]]$token_close + 
                     bracket_pairs[[2]]$close_i / 
                     nchar(txt[bracket_pairs[[2]]$token_close]),
                   ind = 1:nrow(bracket_pairs[[2]]))
  )
  bracket_pos <- bracket_pos[order(bracket_pos$pos), c("type", "ind")]
  bracket_pos$score <- as.numeric(bracket_pos$type == "o") - as.numeric(bracket_pos$type == "c")
  bracket_pos$cumscore <- cumsum(bracket_pos$score)

  bracket_match <- data.frame(do.call(rbind, lapply(1:max(bracket_pos$cumscore), function(mcs){
    cbind(o = bracket_pos$ind[which(bracket_pos$type == "o" & bracket_pos$cumscore == mcs)], 
          c = bracket_pos$ind[which(bracket_pos$type == "c" & bracket_pos$cumscore == (mcs - 1))]) 
  })))
  bracket_match <- bracket_match[order(bracket_match$o),]
  
  # bracket_pairs <- do.call(cbind, bracket_pairs)
  bracket_pairs <- do.call(rbind, apply(bracket_match, 1, function(i) cbind(bracket_pairs[[1]][i[1],], bracket_pairs[[2]][i[2],])))
  
  #adjust tokens to reflect text modifications
  for(i in 1:nrow(bracket_pairs)){
    if(bracket_pairs$token_open[i] == bracket_pairs$token_close[i]){
      next()
    } else {
      inds <- bracket_pairs$token_open[i]:bracket_pairs$token_close[i]
      txt[inds] <- sapply(inds, function(j){
        new_txt <- txt[j]
        if(j != bracket_pairs$token_open[i]){
          new_txt <- paste0(bracket_pairs$style[i], new_txt)
        }
        if(j != bracket_pairs$token_close[i]){
          new_txt <- paste0(new_txt, "}")
        }
        new_txt
      })
    }
  }
  
  txt <- gsub(pattern = "\t", replacement = "\\t", x = txt, fixed = T)
  
  list(tokens = gsub(x = txt, pattern = "\n", replacement = ""),
       has_newline = grepl(txt, pattern = "\n"))
  
}

#create plot
plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
abline(h = 0:150/50, col = "lightgrey", lwd = 0.5)

#specify data
# set.seed(1)
txt <- sample(replace = T, size = 40, x = c("Akr1c1$^2_i$9" , "Aldh1$_a$2" , "Aqp1" , "Arh$^2$ge$_f$25" , "B4galt4" , "C1qtnf6" , "C1qtnf6" , "Dtx3" ,
         "Dtx3" , "Fbln2" , "Fbn1" , "Gatm" , "Gls" , "Gprc5b" , "LOC100910255" , "Lox" , "Lox" ,
         "Ltbp3" , "Ndrg4" , "Pdlim7" , "Sparc" , "Sparc"))

#specify alternate data
txt_string <- "Lorem ipsum $dolor_{50}$ sit\n $amet$, \textbf{consectetur \textit{adipiscing} elit}, $sed_1 do$ \textbf{eiusmod tempor incididunt} ut \textbf{\textit{labore et}} dolore amana aliqua."
tokens <- string_to_tokens(txt_string)
txt <- tokens$tokens

rect_coords <- list(x0=0.275, x1=1.5195, y0 = 2.751, y1 = 0.11547)

#perform analysis
optimal_word_placement_inf <- find_optimal_cex_and_lines(txt = txt, rect_coords = rect_coords, fixed_cex = NA,
                                                         newlines = tokens$has_newline)

#check plot
rect(xleft = rect_coords$x0, ybottom = rect_coords$y0, xright = rect_coords$x1, ytop = rect_coords$y1, lwd = 1, col = "green")
# par(xpd = NA)
text_wrapped_words(txt, rect_coords, optimal_word_placement_inf, justified = F, col = "red", 
                   str_width_lefter_start_ratio = 0, vertically_justified = F) #TODO justification is screwed up :/


plot_col_version <- F
if(plot_col_version){
  plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  rect(xleft = rect_coords$x0, ybottom = rect_coords$y0, xright = rect_coords$x1, ytop = rect_coords$y1, lwd = 1, col = "black")
  cols_list <- lapply(txt, function(word_i) sample((0:3)[-2], nchar(word_i), replace = T))
  text_wrapped_words(txt, rect_coords, optimal_word_placement_inf, justified = T, multicolor_words = T, cols_list = cols_list)
}

#individually colored letters
# text_cols <- function(string, cols, x, y, cex = 1, ...){
#   for(char_i in 1:nchar(string)){
#     txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
#     text(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = txt_cols[char_i], cex = cex, ...)
#   }
# }
# 
# text_cols_bad <- function(string, cols, x, cex = 1, ...){
#   for(char in 1:nchar(string)){
#     text(labels = substr(string, char, char), col = cols[char], cex = cex, x = x + strwidth(substr(string, 1, char)) * cex, ...)
#   }
# }
# 
# plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
# txt <- "abcdefg"
# txt_cols <- 1:nchar(txt)
# text_cols_bad(x = 0.5, y = 2, string = txt, cols = txt_cols, cex = 2, pos = 4)
# text(x = 0.5, y = 2, srt = 45, labels = "text_cols_bad()", pos = 3)
# text_cols(x = 0.5, y = 1, string = txt, cols = txt_cols, cex = 2, pos = 4)
# text(x = 0.5, y = 1.1, srt = 45, labels = "text_cols()", pos = 3)
