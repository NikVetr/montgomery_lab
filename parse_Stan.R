library(dplyr)

stan_code <- readLines("~/repos/egtex-ase/Stan/models/beta-binomial-joint-conc-loc-indiv-pointmix-noncentered.stan")

clean_stan <- function(stan_code) {

  # Regular expressions for detecting types of lines
  declaration_regex <- "^\\s*(int|real|vector|row_vector|matrix|array)"
  lpmf_lpdf_regex <- "_lp(mf|df)\\s*\\("
  target_assignment_regex <- "^\\s*target\\s*(\\+=|=)"
  
  # create a data object to modify
  lines <- stan_code
  
  # Remove comments (text after //)
  lines <- gsub("//.*", "", lines)
  
  #trim whitespace
  lines <- trimws(lines)
  
  # Remove empty lines
  lines <- lines[nchar(lines) > 0]
  
  # Ensure closing brackets are on their own lines
  straggling_close_brackets <- grepl("}", lines) & nchar(lines) != 1
  if(any(straggling_close_brackets)){
    close_bracket_inds <- which(straggling_close_brackets)
    n_lines_inserted <- 0
    for(i in 1:sum(straggling_close_brackets)){
      current_line <- (close_bracket_inds[i] + n_lines_inserted)
      insert_lines <- strsplit(x = lines[current_line], 
                               split = paste0("(?=", "}", ")"), 
                               perl = TRUE)[[1]]
      lines <- append(lines, insert_lines, current_line)
      lines <- lines[-current_line]
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
    }
  }
  
  # Ensure opening brackets are not followed by anything
  # ie only one opening bracket per line
  open_bracket_lines <- grep("\\{", lines)
  straggling_open_bracket_inds <- open_bracket_lines[
    sapply(open_bracket_lines, function(i){
      sum(gregexpr("\\{", lines[i])[[1]] > 0) > 1
    })]
  
  if(length(straggling_open_bracket_inds) > 0){
    n_lines_inserted <- 0
    for(i in 1:length(straggling_open_bracket_inds)){
      current_line <- (straggling_open_bracket_inds[i] + n_lines_inserted)
      insert_lines <- trimws(strsplit(x = lines[current_line], 
                               split = paste0("(?=", "{", ")"), 
                               perl = TRUE)[[1]])
      
      if(insert_lines[1] == "{"){
        insert_lines <- c("", insert_lines)
      }
      
      last_bit_exists <- F
      if(insert_lines[length(insert_lines)] != "{"){
        last_bit <- insert_lines[length(insert_lines)]
        last_bit_exists <- T
      }
      
      insert_lines <- trimws(paste(insert_lines[which(insert_lines == "{")-1],
             insert_lines[which(insert_lines == "{")]
             ))
      
      if(last_bit_exists){
        insert_lines <- c(insert_lines, last_bit)
      }
      
      lines <- append(lines, insert_lines, current_line)
      lines <- lines[-current_line]
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
    }
  }
  
  # get multiline statements on the same line
  # if it's not a "} or ending in a "{" or ";",
  # concatenate lines up to the most recent line after a ";"
  continuation_lines <- !(substr(lines, nchar(lines), nchar(lines)) %in% c("{", ";", "}"))
  if(any(continuation_lines)){
    rlel <- rle(continuation_lines)
    rle_df <- data.frame(val = rlel$values, 
                         start = cumsum(c(0, rlel$lengths)[-length(rlel$lengths)]) + 1,
                         stop = cumsum(rlel$lengths))
    rle_df <- rle_df[rle_df$val,]
    n_lines_removed <- 0
    for(i in 1:nrow(rle_df)){
      line_bounds <- unlist(rle_df[i,2:3] + c(0,1)) - n_lines_removed
      remove_lines <- lines[line_bounds[1]:line_bounds[2]]
      insert_line <- paste0(remove_lines, collapse = " ")
      lines <- append(lines, insert_line, line_bounds[2])
      lines <- lines[-(line_bounds[1]:line_bounds[2])]
      n_lines_removed <- n_lines_removed + length(remove_lines) - 1
    }
  }
  
  #split variable declaration and assignment
  assignment_and_declaration_inds <- which(grepl(declaration_regex, lines) & 
                                             grepl("=", gsub("<[^>]*>", "BOUND_PLACEHOLDER", lines)))
  if(length(assignment_and_declaration_inds) > 0){
    n_lines_inserted <- 0
    for(i in 1:length(assignment_and_declaration_inds)){
      current_line <- assignment_and_declaration_inds[i] + n_lines_inserted
      line <- lines[current_line]
      
      rhs <- paste0(strsplit(gsub("<[^>]*>", "BOUND_PLACEHOLDER", line), "(?==)", perl = T)[[1]][-1], collapse = "")
      lhs <- trimws(substr(line, 1, nchar(line) - nchar(rhs)))
      insert_lines <- c(paste0(lhs, ";"), 
                        paste(tail(strsplit(lhs, " ")[[1]], 1), rhs))
      lines <- append(lines, insert_lines, current_line)
      lines <- lines[-current_line]
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
    }
    
  }
  
  return(lines)
}

stan_code <- clean_stan(stan_code)

parse_stan_blocks <- function(stan_code) {
  
  # Split the code into lines
  lines <- stan_code
  
  # Initialize variables
  current_block <- NULL
  blocks <- list()
  bracket_count <- 0  # Counter for open curly brackets
  block_regex <- "^\\s*(data|transformed data|parameters|transformed parameters|model|functions|generated quantities)\\s*\\{"
  
  # Iterate over the lines
  for (line in lines) {
  
    # Check for block start
    if (grepl(block_regex, line)) {
      current_block <- sub("^\\s*(\\w+(\\s+\\w+)?)\\s*\\{.*", "\\1", line)
      blocks[[current_block]] <- character(0)
      bracket_count <- 1
      
    } else if (!is.null(current_block)) {
      
      # Check for open and close curly brackets
      bracket_count <- bracket_count + 
        sum(gregexpr("\\{", line)[[1]] > 0) - 
        sum(gregexpr("\\}", line)[[1]] > 0)
        
      if (bracket_count == 0) {
        current_block <- NULL
      } else {
        # If we are inside a block, add the line to the current block
        blocks[[current_block]] <-  c(blocks[[current_block]], line)
      }
    }
  }
  
  return(blocks)
}

block_contents <- parse_stan_blocks(stan_code)
block_contents <- block_contents[sapply(block_contents, length) > 0]

parse_declaration <- function(line) {
  # Initialize variables to store extracted information
  declaration_type <- ""
  var_name <- ""
  bounds <- ""
  dimensions <- ""
  distribution <- NA  # Not applicable for declaration lines
  parameters <- NA         # Not applicable for declaration lines
  
  # Regular expressions for different variable types
  scalar_regex <- "^(int|real)\\s*(<[^>]*>)?\\s+(\\w+);"
  vector_regex <- "^(vector|row_vector)\\[(.+)\\]\\s+(\\w+);"
  matrix_regex <- "^matrix\\[(.+)\\]\\s+(\\w+);"
  array_regex <- "^array\\[(.+)\\]\\s+(int|real)\\s*(<[^>]*>)?\\s+(\\w+);"
  
  # Check and extract information for scalar types
  if (grepl(scalar_regex, line)) {
    declaration_type <- sub(scalar_regex, "\\1", line)
    bounds <- sub(scalar_regex, "\\2", line)
    var_name <- sub(scalar_regex, "\\3", line)
    dimensions <- "1"  # Scalar
  }
  # Check and extract information for vector types
  else if (grepl(vector_regex, line)) {
    declaration_type <- sub(vector_regex, "\\1", line)
    dimensions <- sub(vector_regex, "\\2", line)
    var_name <- sub(vector_regex, "\\3", line)
  }
  # Check and extract information for matrix types
  else if (grepl(matrix_regex, line)) {
    declaration_type <- "matrix"
    dimensions <- sub(matrix_regex, "\\1", line)
    var_name <- sub(matrix_regex, "\\2", line)
  }
  # Check and extract information for array types
  else if (grepl(array_regex, line)) {
    declaration_type <- sub(array_regex, "\\2", line)
    dimensions <- sub(array_regex, "\\1", line)
    bounds <- sub(array_regex, "\\3", line)
    var_name <- sub(array_regex, "\\4", line)
  }
  
  return(list(declaration_type = declaration_type, var_name = var_name, bounds = bounds, dimensions = dimensions, distribution = distribution, parameters = parameters))
}

parse_assignment <- function(line) {
  # Initialize variables to store extracted information
  lhs_name <- ""
  lhs_index <- NA
  rhs_expression <- ""
  
  # Regular expression for assignment statements
  # This regex captures the LHS (including index), and the RHS expression
  assignment_regex <- "^(\\w+)(\\[.*\\])?\\s*=\\s*(.+);$"
  
  # Check if the line matches the assignment statement pattern
  if (grepl(assignment_regex, line)) {
    lhs_name <- sub(assignment_regex, "\\1", line)
    lhs_index <- sub(assignment_regex, "\\2", line)
    rhs_expression <- sub(assignment_regex, "\\3", line)
    
    # Removing any additional brackets around the index
    lhs_index <- gsub("\\[|\\]", "", lhs_index)
    
    return(list(lhs_name = lhs_name, lhs_index = lhs_index, 
                rhs_expression = rhs_expression))
    
  }
  
  return(list(lhs_name = lhs_name, lhs_index = lhs_index, rhs_expression = rhs_expression))
}



parse_operation <- function(expression) {
  # Replace Stan-specific operations with R equivalents
  expression <- gsub("\\.\\*", "*", expression)  # Element-wise multiplication
  expression <- gsub("\\./", "/", expression)    # Element-wise division
  
  
  # Handle some of Stan's built-in functions, translating them to R equivalents
  # Example: Stan's normal() becomes R's dnorm()
  # Add similar replacements for other Stan functions as needed
  expression <- gsub("normal\\(", "dnorm(", expression)
  expression <- gsub("bernoulli_logit\\(", "dbinom(", expression)
  
  # ... [Additional replacements for other Stan-specific functions and operations]
  
  return(expression)
}

parse_loop <- function(line) {
  loop_type <- NA
  loop_variable <- NA
  loop_range <- NA
  loop_action <- if (grepl("^\\s*\\}", line)) "close" else "open"
  
  # Handle opening of loops
  if (loop_action == "open") {
    for_loop_regex <- "^for\\s*\\((\\w+)\\s+in\\s+(.+)\\)\\s*\\{"
    while_loop_regex <- "^while\\s*\\((.+)\\)\\s*\\{"
    
    if (grepl(for_loop_regex, line)) {
      loop_type <- "for"
      loop_variable <- sub(for_loop_regex, "\\1", line)
      loop_range <- sub(for_loop_regex, "\\2", line)
    } else if (grepl(while_loop_regex, line)) {
      loop_type <- "while"
      loop_range <- sub(while_loop_regex, "\\1", line)
    }
  }
  
  return(list(loop_type = loop_type, loop_variable = loop_variable, loop_range = loop_range, loop_action = loop_action))
}

parse_sampling <- function(line) {
  # Initialize variables to store extracted information
  variable_name <- ""
  distribution <- ""
  parameters <- ""
  declaration_type <- NA  # Not applicable for sampling lines
  bounds <- NA           # Not applicable for sampling lines
  dimensions <- NA       # Not applicable for sampling lines
  
  # Regular expression for sampling statements
  sampling_regex <- "^(\\w+)\\s*~\\s*(\\w+)\\(([^)]*)\\);"
  
  # Check if the line matches the sampling statement pattern
  if (grepl(sampling_regex, line)) {
    variable_name <- sub(sampling_regex, "\\1", line)
    distribution <- sub(sampling_regex, "\\2", line)
    parameters <- sub(sampling_regex, "\\3", line)
  }
  
  return(list(declaration_type = declaration_type, var_name = variable_name, bounds = bounds, dimensions = dimensions, distribution = distribution, parameters = parameters))
}

parse_lpmf_lpdf <- function(line) {
  # Extract the RHS of the assignment
  rhs <- sub(".*=\\s*(.+);", "\\1", line)
  
  # Isolate the _lpmf or _lpdf function call
  lpmf_lpdf_call <- regmatches(rhs, regexpr("\\w+_lp(mf|df)\\s*\\([^\\)]+\\)", rhs))
  
  # Initialize variables for extracted information
  distribution <- ""
  outcome_variable <- ""
  parameters <- ""
  mixture_probability <- NA  # Default to 0
  
  # Check if there's a mixture probability component
  if (grepl("\\+", rhs)) {
    mixture_probability <- sub("\\s*\\+\\s*\\w+_lp(mf|df).*", "", rhs)
    mixture_probability <- sub("log\\s*\\((.+)\\)", "\\1", mixture_probability)
  }
  
  # Process the _lpmf or _lpdf call
  if (length(lpmf_lpdf_call) > 0) {
    # Extract the distribution
    distribution <- sub("(.+)_lp(mf|df).*", "\\1", lpmf_lpdf_call)
    
    # Extract the outcome variable and parameters
    args <- sub(".*\\(([^|]+)\\s*\\|\\s*(.+)\\)", "\\1,\\2", lpmf_lpdf_call)
    args_split <- strsplit(args, ",")[[1]]
    if (length(args_split) >= 2) {
      outcome_variable <- args_split[1]
      parameters <- paste(args_split[-1], collapse = ", ")
    }
  }
  
  return(list(distribution = distribution, outcome_variable = outcome_variable, 
              parameters = parameters, mixture_probability = mixture_probability))
}



parse_target_assignment <- function(line) {
  
  # Initialize variables for extracted information
  lhs <- "target"
  rhs <- ""
  
  # Regular expression patterns
  target_plus_equal_regex <- "^\\s*target\\s*\\+=\\s*(.+);"
  target_equal_regex <- "^\\s*target\\s*=\\s*target\\s*\\+\\s*(.+);"
  
  # Check for 'target += ...' pattern
  if (grepl(target_plus_equal_regex, line)) {
    rhs <- sub(target_plus_equal_regex, "\\1", line)
  }
  # Check for 'target = target + ...' pattern
  else if (grepl(target_equal_regex, line)) {
    rhs <- sub(target_equal_regex, "\\1", line)
  }
  
  # Further processing if RHS contains lpmf or lpdf
  parsed_info <- if (grepl("_lp(mf|df)\\s*\\(", rhs)) {
    parse_lpmf_lpdf(rhs)
  } else {
    list(rhs_expression = rhs)
  }
  parsed_info["lhs"] <- lhs
  
  return(parsed_info)
}

parse_stan_lines <- function(block_content) {
  
  # create a modifiable lines object
  lines <- block_content
  
  if(length(lines) == 0) return(lines)
  
  # Initialize a data frame to hold line information
  line_info <- data.frame(line = character(), type = character(), loop_depth = integer(), loop_start_line = integer(), stringsAsFactors = FALSE)
  
  # Initialize loop depth counter and stack
  loop_depth <- 0
  open_brackets_stack <- list()
  
  # Regular expressions for detecting types of lines
  declaration_regex <- "^\\s*(int|real|vector|row_vector|matrix|array)"
  lpmf_lpdf_regex <- "_lp(mf|df)\\s*\\("
  target_assignment_regex <- "^\\s*target\\s*(\\+=|=)"
  
  # Iterate over the lines
  for (i in 1:length(lines)) {
    line <- lines[i]
    
    # Initialize variables to store parsed information
    parsed_type <- ""
    parsed_info <- list()
    
    # Determine the type of line and parse accordingly
    if (grepl(declaration_regex, line)) {
      parsed_type <- "declaration"
      parsed_info <- parse_declaration(line)
    } else if (grepl(target_assignment_regex, line)) {
      parsed_type <- "sampling_target"
      parsed_info <- parse_target_assignment(line)
    } else if (grepl(lpmf_lpdf_regex, line)) {
      parsed_type <- "lpmf_lpdf"
      parsed_info <- parse_lpmf_lpdf(line)
    } else if (grepl("=", line)) {
      parsed_type <- "assignment"
      parsed_info <- parse_assignment(line)
    } else if (grepl("~", line)) {
      parsed_type <- "sampling"
      parsed_info <- parse_sampling(line)
    } else if (grepl("^for\\s*\\(|^while\\s*\\(", line)) {
      parsed_type <- "loop"
      parsed_info <- parse_loop(line)
      loop_depth <- loop_depth + 1
      open_brackets_stack[[length(open_brackets_stack) + 1]] <- i
    } else if (grepl("^\\s*\\}", line)) {
      if (length(open_brackets_stack) > 0) {
        start_line <- open_brackets_stack[[length(open_brackets_stack)]]
        parsed_type <- "loop"
        parsed_info <- list(loop_action = "close", 
                            start_line = start_line,
                            loop_type = line_info$loop_type[start_line])
        open_brackets_stack <- open_brackets_stack[-length(open_brackets_stack)]
        loop_depth <- loop_depth - 1
      } else {
        parsed_type <- "operation"
      }
    } else {
      parsed_type <- "other"
    }
    
    # Add parsed information to the data frame
    line_info <- bind_rows(line_info, 
                           c(list(line = line, 
                                  type = parsed_type, 
                                  loop_depth = loop_depth - ifelse(parsed_type == "loop" && parsed_info$loop_action == "open", 1, 0)
                                  ), 
                             parsed_info))
  }
  
  return(line_info)
}

parsed_lines <- lapply(block_contents, parse_stan_lines)
parsed_lines <- lapply(setNames(names(parsed_lines), names(parsed_lines)), 
                       function(block_name){
                         cbind(parsed_lines[[block_name]], block_name = block_name)
                       })

parsed_lines <- bind_rows(parsed_lines)

#### interpreting lines ####

interpret_assignment <- function(line, block_name, dat, samps) {
  # Placeholder for interpret_assignment function
  return("interpret_assignment not yet implemented")
}

interpret_operation <- function(line, block_name, dat, samps) {
  # Placeholder for interpret_operation function
  return("interpret_operation not yet implemented")
}

interpret_loop <- function(line, block_name, dat, samps) {
  # Placeholder for interpret_loop function
  return("interpret_loop not yet implemented")
}

interpret_sampling <- function(line, block_name, dat, samps) {
  # Placeholder for interpret_sampling function
  return("interpret_sampling not yet implemented")
}

interpret_declaration <- function(line, block_name, dat, samps) {
  if (block_name == "parameters") {
    # Read from MCMC samples
    var_name <- line$var_name
    return(paste0(var_name, " <- subset_samps(include = '", var_name, "', samps = samps)"))
  } else if (block_name == "data") {
    # Read from input data
    var_name <- line$var_name
    return(paste0(var_name, " <- dat$", var_name))
  } else {
    # Initialize the variable in R
    var_name <- line$var_name
    return(paste0(var_name, " <- NULL  # Initialize variable"))
  }
}

interpret_lines <- function(line, dat, samps) {
  switch(line$type,
         "declaration" = interpret_declaration(line, line$block_name, dat, samps),
         "assignment" = interpret_assignment(line, line$block_name, dat, samps),
         "operation" = interpret_operation(line, line$block_name, dat, samps),
         "loop" = interpret_loop(line, line$block_name, dat, samps),
         "sampling" = interpret_sampling(line, line$block_name, dat, samps),
         "Unknown type")
}

# Example usage
processed_code <- lapply(parsed_lines, interpret_lines, dat, samps)

