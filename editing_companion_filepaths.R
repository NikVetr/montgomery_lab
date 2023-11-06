#useful functions
extract_paste0 <- function(txt){
  
  #identify location of script keyword? 
  #identify location of paste0( string
  string <- paste0(txt, collapse = "")
  paste0_locs <- unlist(gregexpr("paste0", string))
  
  #check we have paste0s in here
  if(all(paste0_locs == -1)){
    return("")
  }
  
  #process last paste0, then swap it in w/ a placeholder and call this function again
  last_paste0 <- tail(paste0_locs, 1)
  open_paren_locs <- unlist(gregexpr("\\(", string))
  first_open_paren <- open_paren_locs[min(which(open_paren_locs > last_paste0))]
  
  #check we have an open paren in here
  if(all(is.na(open_paren_locs))){
    return("")
  }
  
  #now find the corresponding close paren
  #initialize counter to 1, 
  #count next open paren as +1
  #count next closing paren as -1
  #when counter reaches 0, extract contents between parens
  close_paren_locs <- unlist(gregexpr(")", string))
  
  if(all(close_paren_locs == -1) | all(close_paren_locs <= first_open_paren)){
    return("")
  }
  
  paren_pos <- rep(0, tail(close_paren_locs, 1) - first_open_paren)
  paren_pos[(open_paren_locs - first_open_paren + 1)[sign(open_paren_locs - first_open_paren + 1) == 1]] <- 1
  paren_pos[(close_paren_locs - first_open_paren + 1)[sign(close_paren_locs - first_open_paren + 1) == 1]] <- -1
  corresp_close_paren_disp <- which(cumsum(paren_pos) == 0)
  if(length(corresp_close_paren_disp) != 0){
    corresp_close_paren <- first_open_paren + corresp_close_paren_disp - 1  
  } else {
    return("")
  }
  
  #find commas and quotes
  #then find all commas not inside of quotes
  target_substr <- substr(string, first_open_paren, corresp_close_paren)
  quote_char <- "\""
  quote_locs <- unlist(gregexpr(quote_char, target_substr))
  
  #check we have double quotes
  if(all(quote_locs == -1)){
    quote_char <- "\'"
    quote_locs <- unlist(gregexpr(quote_char, target_substr))
  }
  
  if(all(quote_locs == -1)){
    return("")
  }
  
  quote_locs <- t(matrix(quote_locs, nrow = 2))
  # quoted_content <- apply(quote_locs, 1, function(x) substr(target_substr, x[1] + 1, x[2] - 1))
  # comma_locs <- unlist(gregexpr(",", target_substr))
  
  #actually let's keep this simple
  #and just return quoted content separated by *
  # last_paste0_contents <- paste0(quoted_content, collapse = "*")
  
  #scratch that, this runs afoul ifelse statements or more complicated constructions
  #we *do* need something like comma_locs, and to only extract content not in additional functions
  #one way to do that could be to delete all internal parentheticals not preceded by ws, instead of mucking with commas
  internal_open_parens <- setdiff(unlist(gregexpr("\\(", target_substr)), 1)
  internal_close_parens <- setdiff(unlist(gregexpr(")", target_substr)), nchar(target_substr)) 
  
  if(length(internal_open_parens) != 0 | length(internal_close_parens) != 0){
    
    if(length(internal_open_parens) == 0 | length(internal_close_parens) == 0){
      return("")
    } #check we have both open and close parens
    
    #check these aren't in quotes
    internal_open_parens <- internal_open_parens[!sapply(internal_open_parens, function(paren_i){
      lq <- max(which(quote_locs[,1] < paren_i))
      quote_locs[lq,2] > paren_i
    })]
    
    internal_close_parens <- internal_close_parens[!sapply(internal_close_parens, function(paren_i){
      lq <- max(which(quote_locs[,1] < paren_i))
      quote_locs[lq,2] > paren_i
    })]
    
    #now find sequences inside these internal parens and excise them
    internal_paren_vec <- rep(0, nchar(target_substr))
    internal_paren_vec[internal_open_parens] <- 1
    internal_paren_vec[internal_close_parens + 1] <- -1
    chars_outside_parens <- cumsum(internal_paren_vec) == 0
    target_substr <- paste0(strsplit(target_substr, "")[[1]][chars_outside_parens], collapse = "")
    
  }

  #now process new quote locations, and compile output string
  quote_locs <- unlist(gregexpr(quote_char, target_substr))
  quote_locs <- t(matrix(quote_locs, nrow = 2))
  quoted_content <- apply(quote_locs, 1, function(x) substr(target_substr, x[1] + 1, x[2] - 1))
  last_paste0_contents <- paste0(quoted_content, collapse = "*")
  
  #return result or call function again
  if(length(paste0_locs) == 1){
    return(last_paste0_contents)
  } else {
    next_outer_string <- c(substr(string, 1, first_open_paren - 7), 
                           paste0(quote_char, last_paste0_contents, quote_char),
                           substr(string, corresp_close_paren + 1, nchar(string)))
    return(extract_paste0(next_outer_string))
  }
  
}

#get filepaths loaded
script_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/scripts_local-paths/"
scripts <- paste0(script_dir, list.files(script_dir, recursive = T))

root_dirs <- c("~/data/smontgom/",
               "~/repos/",
               "/Volumes/SSD500GB/")

data_entry_lines <- do.call(rbind, lapply(setNames(scripts, scripts), function(script_loc){
  
  #subset to script
  print(script_loc)
  script_lines <- readLines(script_loc, warn = F)
  
  #subset to root directory keyword
  dir_inds <- sapply(root_dirs, function(rdir) grep(pattern = rdir, script_lines))
  out <- do.call(rbind, lapply(1:length(dir_inds), function(dir_i){
    dir_name <- names(dir_inds)[dir_i]
    dir_ind <- dir_inds[[dir_i]]
    
    #find matches
    if(length(dir_ind) != 0){
      do.call(rbind, lapply(dir_ind, function(di){
        target_line <- orig_script_line <- script_lines[di]
        
        #extract target file path from quotes or from paste0
        if(grepl("paste", target_line)){
          
          prop_target_line <- extract_paste0(target_line)
          
          if(prop_target_line == ""){
            
            lines_to_try <- 10
            lines_tried <- 0
            
            while(lines_tried <= lines_to_try & prop_target_line == ""){
              
              #increment counter
              lines_tried <- lines_tried + 1
              
              #check we're not out of bounds
              if((di + lines_tried) > length(script_lines)){
               target_line <- NA
               break
              }
              
              #try new line
              prop_target_line <- extract_paste0(script_lines[di + 0:(lines_tried)])
            }
          }
          
          target_line <- prop_target_line
          
        } else {
          
          start_string <- unlist(gregexpr(dir_name, target_line))
          end_quotes <- sort(c(unlist(gregexpr("\"", target_line)), unlist(gregexpr("\'", target_line))))
          end_quote <- end_quotes[min(which(end_quotes > start_string))]
          target_line <- substr(target_line, start_string, end_quote-1)
        
        }
        
        #return result
        data.frame(script_loc = script_loc, dir_name = dir_name, 
              line_ind = di, 
              orig_script_line = orig_script_line, 
              target_line = target_line,
              commented_out = substr(trimws(orig_script_line), 1, 1) == "#")
        
      }))
      
    } else {
      #if no matches
      return(character(0))
    }
    
  }))
  
  out
  
}))
rownames(data_entry_lines) <- NULL

#make multiple consecutive asterisks single asterisks
data_entry_lines$target_line <- gsub("\\*+", "*", data_entry_lines$target_line)

#trim bit of lingering whitespace
data_entry_lines$target_line <- trimws(data_entry_lines$target_line)

#find bits in shell command
data_entry_lines$in_shell <- sapply(1:nrow(data_entry_lines), function(i){
  script_line <- data_entry_lines$orig_script_line[i]
  start_of_dir <- unlist(gregexpr(pattern = data_entry_lines$dir_name[i], script_line))
  quote_locs <- sort(c(unlist(gregexpr("\"", script_line)), 
                       unlist(gregexpr("\'", script_line))))
  quote_locs <- quote_locs[quote_locs > 0]
  quoted_substr <- substr(script_line, 
                          max(quote_locs[start_of_dir > quote_locs]), 
                          min(quote_locs[start_of_dir < quote_locs]))
  quoted_substr <- gsub("\"", "", quoted_substr)
  quoted_substr <- gsub("\'", "", quoted_substr)
  unlist(gregexpr(pattern = data_entry_lines$dir_name[i], quoted_substr)) != 1
})

#and modify them
data_entry_lines$target_line[data_entry_lines$in_shell] <- 
  sapply(1:length(data_entry_lines$target_line[data_entry_lines$in_shell]), function(i){
    linei <- data_entry_lines$target_line[data_entry_lines$in_shell][i]
    start_of_dir <- unlist(gregexpr(data_entry_lines$dir_name[data_entry_lines$in_shell][i], 
                                    linei))
    poss_ends <- c(unlist(gregexpr(" ", linei)), unlist(gregexpr(";", linei)), nchar(linei))
    poss_ends <- sort(poss_ends[poss_ends > 0])
    end_of_dir <- min(poss_ends[poss_ends > start_of_dir])
    if(end_of_dir != nchar(linei)) end_of_dir <- end_of_dir - 1
    out <- trimws(substr(linei, start_of_dir, end_of_dir))
    if(is.na(file.info(out)$isdir)){
      out <- paste0(out, "*")
    } else if(file.info(out)$isdir){
      out <- paste0(out, "/")
      gsub("//", "/", out)
    }
    out <- gsub("\\*+", "*", out)
  })

#if we end in a wildcard following a /, remove wildcard
wild_ends <- sapply(data_entry_lines$target_line, function(x) #or could have used endsWith, doh
  substr(x, nchar(x) - 1, nchar(x))) == "/*"
data_entry_lines$target_line[wild_ends] <- sapply(data_entry_lines$target_line[wild_ends], function(x) 
  substr(x, 1, nchar(x) - 1))

#also if we have a /*/ too
# wild_ends <- sapply(data_entry_lines$target_line, function(x) #or could have used endsWith, doh
#   substr(x, nchar(x) - 1, nchar(x))) == "/*/"
# data_entry_lines$target_line[wild_ends] <- sapply(data_entry_lines$target_line[wild_ends], function(x) 
#   substr(x, 1, nchar(x) - 2))

#if we're in a folder that exclusively contains a wildcard, just specify the entire folder
data_entry_lines$target_line <- gsub("(/\\*/).*", "", data_entry_lines$target_line)

#check if we're a file or a folder; add a '/' to the end of all folders
data_entry_lines$isdir <- file.info(data_entry_lines$target_line)$isdir
data_entry_lines$isdir[is.na(data_entry_lines$isdir)] <- F
data_entry_lines$target_line[data_entry_lines$isdir] <- 
  gsub("/+", "/", paste0(data_entry_lines$target_line[data_entry_lines$isdir], "/"))


#gsub "thickCluster" bc I file.renamed it internally
data_entry_lines$target_line <- gsub(pattern = "thickCluster", replacement = "thincluster", x = data_entry_lines$target_line)

#gsub "~/repos/gtex-pipeline/" bc I moved the file when I ran out of space
data_entry_lines$target_line <- gsub(pattern = "~/repos/gtex-pipeline/", 
                                     replacement = "/Volumes/SSD500GB/gtex-pipeline/", x = data_entry_lines$target_line)

#gsub "~/repos/MetaXcan/" bc I moved the file when I ran out of space
data_entry_lines$target_line <- gsub(pattern = "~/repos/MetaXcan/", 
                                     replacement = "/Volumes/SSD500GB/MetaXcan/", x = data_entry_lines$target_line)


#final error checking
err_keyw <- "//"
data_entry_lines[grepl(pattern = err_keyw, x = data_entry_lines$target_line),]

#compile external filenames to write to disk
# files <- unique(data_entry_lines$target_line[!data_entry_lines$commented_out])
files <- unique(data_entry_lines$target_line)

#construct a map of these to the original data object
data_entry_lines$target_line_subs <- data_entry_lines$target_line

#subset files to exclude subsets, first of folders
folders <- files[endsWith(files, "/")]
top_level_folders <- sapply(files, function(file_i){
  nonself_folders <- setdiff(folders, file_i)
  nonself_folders <- setdiff(folders, root_dirs) #remove top-level dirs
  containing_folders <- nonself_folders[sapply(nonself_folders,  function(folder_i) startsWith(file_i, folder_i))]
  if(length(containing_folders) == 0){
    return(file_i)
  } else{
    return(containing_folders[which.min(nchar(containing_folders))])
  }
})

files <- files[files %in% top_level_folders]
data_entry_lines$target_line_subs <- top_level_folders[data_entry_lines$target_line_subs]

#then wrt wildcards
wilds <- files[grepl("\\*", files)]
wild_redundancies <- sapply(files, function(file_i){
  nonself_wilds <- setdiff(wilds, file_i)
  containing_wilds <- nonself_wilds[sapply(nonself_wilds,  function(wild_i) grepl(glob2rx(wild_i), file_i))]
  if(length(containing_wilds) == 0){
    return(file_i)
  } else{
    return(containing_wilds[which.min(nchar(containing_wilds))])
  }
})
files <- files[files == wild_redundancies]
data_entry_lines$target_line_subs <- wild_redundancies[data_entry_lines$target_line_subs]

#remove paths already in the correct directory, bc we can move it externally
files <- files[!startsWith(files, "~/repos/MoTrPAC_Complex_Traits/")]
# data_entry_lines <- data_entry_lines[data_entry_lines$target_line %in% files,]
"/Volumes/2TB_External/MoTrPAC_Complex_Traits/" #is the new location for these
data_entry_lines$target_line_subs[
  startsWith(data_entry_lines$target_line_subs, prefix = "~/repos/MoTrPAC_Complex_Traits/")] <- 
  "~/repos/MoTrPAC_Complex_Traits/"

#check that everything is here
all(data_entry_lines$target_line_subs[!startsWith(data_entry_lines$target_line_subs, "~/repos/MoTr")] %in% files)
# data_entry_lines[!(data_entry_lines$target_line_subs %in% files),]$target_line_subs

#finally write to disk
# write.csv(files, file = "~/data/smontgom/pass1b_companion_paper_external_filenames.csv")

#process stragglers
extra_files <- setdiff(files, unlist(read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_filenames.csv")[,2]))
# write.csv(extra_files, file = "~/data/smontgom/pass1b_companion_paper_external_extra_filenames.csv")

#read mapping of files back in and check for validity
files_map <- rbind(read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_filenames_map.csv")[,c("old_path", "new_dir")],
                   read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_extra_filenames_map.csv")[,c("old_path", "new_dir")])
files_map <- rbind(files_map, 
                   cbind(old_path = "~/repos/MoTrPAC_Complex_Traits/", 
                         new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/"),
                   cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_genotypes/GTEx_v8*", 
                         new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"),
                   cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_data/GTEx_Analysis_v8_*_expected_count.gct.gz", 
                         new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"),
                   cbind(old_path = "/Volumes/SSD500GB/gtex-pipeline/expression_data/GTEx_Analysis_v8_*_tpm.gct.gz", 
                         new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"))
all(data_entry_lines$target_line_subs %in% files_map$old_path)
data_entry_lines[!(data_entry_lines$target_line_subs %in% files_map$old_path),]

data_entry_lines$new_dir_name <- setNames(files_map$new_dir, files_map$old_path)[
  data_entry_lines$target_line_subs]

#now preprocess gsub pattern + replacement
data_entry_lines$gsub_from <- data_entry_lines$dir_name
data_entry_lines$gsub_to <- data_entry_lines$new_dir_name
data_entry_lines$gsub_to[startsWith(data_entry_lines$target_line, "~/repos/MoTrPAC_Complex_Traits/")] <- "/Volumes/2TB_External/"
data_entry_lines$new_script_line <- sapply(1:nrow(data_entry_lines), function(i) gsub(pattern = data_entry_lines$gsub_from[i], 
                                                  replacement = data_entry_lines$gsub_to[i], 
                                                  x = data_entry_lines$orig_script_line[i]))
data_entry_lines$new_target <- sapply(1:nrow(data_entry_lines), function(i) gsub(pattern = data_entry_lines$gsub_from[i], 
                                                                                      replacement = data_entry_lines$gsub_to[i], 
                                                                                      x = data_entry_lines$target_line[i]))
file.exists(data_entry_lines$new_target)

wildfile.exists <- function(path, num = F){
  out <- sapply(path, function(path_i){
    parent_folder <- dirname(path_i)
    internal_files <- list.files(parent_folder, full.names = T)
    matches <- grepl(glob2rx(path_i), internal_files)
    if(num){return(sum(matches))}else{return(any(matches))}
  })
  names(out) <- NULL
  out
}

new_file_found <- cbind(wild = wildfile.exists(data_entry_lines$new_target), 
                        reg = file.exists(data_entry_lines$new_target))
data_entry_lines[!apply(new_file_found, 1, any),]

#now substitute in the new lines for the old lines
lineswap_split <- split(data_entry_lines, data_entry_lines$script_loc)
for(script_path in names(lineswap_split)){
  
  print(script_path)
  
  lineswap <- lineswap_split[[script_path]]
  script_lines <- readLines(script_path, warn = F)
  
  if(all(script_lines[lineswap$line_ind] == lineswap$orig_script_line)){
    script_lines[lineswap$line_ind] <- lineswap$new_script_line
    writeLines(script_lines, script_path)
  } else {
    print(paste0("ERROR: lines do not match in <<", script_path, ">>"))
  }
  
}
