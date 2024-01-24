library(dplyr)
library(cmdstanr)
library(posterior)

stan_code <- readLines("~/repos/egtex-ase/Stan/models/beta-binomial-joint-conc-loc-indiv-pointmix-noncentered.stan") #load stan model spec
load("~/repos/egtex-ase/Stan/output/beta-binomial-joint-conc-loc-indiv-pointmix-noncentered.cmdStanR.fit") #load Stan output
samps <- data.frame(as_draws_df(out$draws()))
#load input files json

subset_samps <- function(var_name, samps) {
  # Regex pattern to match either exact variable name or name followed by indices
  pattern <- paste0("^(", var_name, "(\\.[0-9]+)*\\.|", var_name, ")$")
  matched_cols <- grep(pattern, colnames(samps), value = TRUE)
  return(samps[, matched_cols, drop = F])
}

munge_samps <- function(var_name, df) {
  # Check if the data frame has only one row
  if (nrow(df) == 1) {
    # Process a single row to construct the appropriate data structure
    # Determine the structure of the data (scalar, vector, matrix, or array)
    col_names <- colnames(df)
    col_names_without_var <- gsub(paste0("^", var_name, "\\."), "", col_names)
    if (all(grepl("\\.\\d+\\.", col_names))) {
      # Handle vectors, matrices, and arrays
      # Extract indices from column names
      indices <- lapply(strsplit(col_names_without_var, "\\."), function(x) as.numeric(x))
      
      # Determine if it's a vector, matrix, or array
      max_dims <- sapply(indices, length)
      if (all(max_dims == 1)) {
        # Vector
        return(unlist(df))
      } else if (all(max_dims == 2)) {
        # Matrix
        dim_order <- do.call(rbind, indices)
        mat <- matrix(NA, nrow = max(dim_order[,1]), ncol = max(dim_order[,2]))
        mat[cbind(dim_order[,1], dim_order[,2])] <- unlist(df)
        return(mat)
      } else {
        # Array (higher dimensions)
        array_dims <- apply(do.call(rbind, indices), 2, max)
        arr <- array(NA, dim = array_dims)
        array_indices <- do.call(expand.grid, lapply(array_dims, seq_len))
        arr[as.matrix(array_indices)] <- unlist(df)
        return(arr)
      }
    } else {
      # Scalar or simple vector
      return(unlist(df))
    }
  } else {
    # Recursively apply to each row
    return(lapply(seq_len(nrow(df)), function(i) munge_samps(var_name = var_name, df = df[i, , drop = FALSE])))
  }
}

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

interpret_assignment <- function(line) {
  # Construct the LHS of the assignment
  lhs <- line$lhs_name
  if (nchar(line$lhs_index) > 0) {
    lhs <- paste0(lhs, "[", line$lhs_index, "]")
  }
  
  # Process the RHS of the assignment
  rhs <- interpret_operation(line$rhs_expression)
  
  # Combine LHS and RHS to form the assignment statement
  return(paste0(lhs, " <- ", rhs))
}


interpret_operation <- function(stan_expression) {
  # Replace specific Stan functions with equivalent R functions or code
  # Example: Implementing inv_logit
  stan_expression <- gsub("inv_logit\\(", "plogis(", stan_expression)
  
  # Handle element-wise operations (*, /, +, -) - In R, these are by default element-wise
  # Stan's element-wise operations (e.g., .*) are the same as R's default, so we can remove the '.'
  stan_expression <- gsub("\\.\\*", "*", stan_expression)
  stan_expression <- gsub("\\./", "/", stan_expression)
  stan_expression <- gsub("\\.\\+", "+", stan_expression)
  stan_expression <- gsub("\\.-", "-", stan_expression)
  
  # Ignore certain functions if needed, e.g., to_vector (if it's redundant in R)
  # This can be a no-op if to_vector has no side effects in R context
  stan_expression <- gsub("to_vector\\(", "(", stan_expression)
  
  # Return the translated stan_expression
  return(stan_expression)
}


interpret_loop <- function(line) {
  if (line$loop_action == "open") {
    # Construct the loop header
    loop_variable <- line$loop_variable
    loop_range <- line$loop_range
    loop_header <- paste0("for (", loop_variable, " in ", loop_range, ") {")
    return(loop_header)
  } else if (line$loop_action == "close") {
    # Close the loop
    return("}")
  } else {
    # Handle unexpected cases
    return("# Unexpected loop structure")
  }
}


interpret_sampling <- function(line, dat, samps, sample_index = 1, post_pred_sim = TRUE, sim = TRUE, bound_declarations = NA) {
  # LHS of the sampling statement
  lhs <- line$var_name
  
  # Check if the LHS is in the data block and handle accordingly
  if (lhs %in% names(dat) && !post_pred_sim) {
    if (sim) {
      # Sample from the specified distribution (prior predictive simulation)
      return(paste0(lhs, " <- ", generate_sampling_code(line, sim)))
    } else {
      # Find the density or mass for the parameter in the specified distribution
      return(paste0(lhs, " <- ", generate_density_code(line, dat)))
    }
  } else {
    if (post_pred_sim) {
      # Posterior predictive simulation
      sampled_param <- munge_samps(line$var_name, subset_samps(line$var_name, samps[sample_index, , drop = FALSE]))
      return(paste0(lhs, " <- ", sampled_param))
    } else {
      # Prior predictive simulation
      return(paste0(lhs, " <- ", generate_sampling_code(line, sim)))
    }
  }
}

generate_sampling_code <- function(line) {
  distribution <- line$distribution
  parameters <- line$parameters
  
  # Mapping Stan distributions to R functions
  distribution_map <- list(
    std_normal = "rnorm(n = 1, mean = 0, sd = 1)",
    normal = "rnorm(n = 1, mean = param1, sd = param2)",  # Assuming 'mean, sd' parameterization
    beta = "rbeta(n = 1, shape1 = param1, shape2 = param2)",
    binomial = "rbinom(n = 1, size = param1, prob = param2)",
    gamma = "rgamma(n = 1, shape = param1, rate = param2)",  # Assuming 'shape, rate' parameterization
    student_t = "rt(n = 1, df = param1, ncp = param2)",
    double_exponential = "rlaplace(n = 1, location = param1, scale = param2)",
    exponential = "rexp(n = 1, rate = param1)",
    lognormal = "rlnorm(n = 1, meanlog = param1, sdlog = param2)",
    chi_square = "rchisq(n = 1, df = param1)",
    uniform = "runif(n = 1, min = param1, max = param2)",
    poisson = "rpois(n = 1, lambda = param1)"
    # Add other distributions as needed
  )
  
  # Generate R code for sampling
  r_function <- distribution_map[[distribution]]
  if (is.null(r_function)) {
    return(paste("# No R equivalent for", distribution, "distribution"))
  } else {
    # Replace placeholders with actual parameters
    param_names <- paste0("param", seq_along(strsplit(parameters, ",")[[1]]))
    for (i in seq_along(param_names)) {
      r_function <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], r_function)
    }
    return(r_function)
  }
}


generate_density_code <- function(line, dat) {
  distribution <- line$distribution
  parameters <- line$parameters
  var_name <- line$var_name
  
  # Mapping Stan distributions to R density/mass functions
  density_map <- list(
    normal = "dnorm(x = dat$VAR_NAME, mean = param1, sd = param2)",
    beta = "dbeta(x = dat$VAR_NAME, shape1 = param1, shape2 = param2)",
    binomial = "dbinom(x = dat$VAR_NAME, size = param1, prob = param2)",
    gamma = "dgamma(x = dat$VAR_NAME, shape = param1, rate = param2)",
    student_t = "dt(x = dat$VAR_NAME, df = param1, ncp = param2)",
    double_exponential = "dexp(x = dat$VAR_NAME, rate = param1)",  # Assuming 'rate' parameterization
    exponential = "dexp(x = dat$VAR_NAME, rate = param1)",
    lognormal = "dlnorm(x = dat$VAR_NAME, meanlog = param1, sdlog = param2)",
    chi_square = "dchisq(x = dat$VAR_NAME, df = param1)",
    uniform = "dunif(x = dat$VAR_NAME, min = param1, max = param2)",
    poisson = "dpois(x = dat$VAR_NAME, lambda = param1)"
    # Add other distributions as needed
  )
  
  # Generate R code for density/mass calculation
  r_function <- density_map[[distribution]]
  if (is.null(r_function)) {
    return(paste("# No R equivalent for", distribution, "distribution"))
  } else {
    # Replace placeholder 'VAR_NAME' with the actual variable name from the data
    r_function <- gsub("VAR_NAME", var_name, r_function)
    
    # Replace placeholders with actual parameters
    param_names <- paste0("param", seq_along(strsplit(parameters, ",")[[1]]))
    for (i in seq_along(param_names)) {
      r_function <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], r_function)
    }
    return(r_function)
  }
}



interpret_declaration <- function(line, dat, samps) {
  if (line$block_name == "parameters") {
    # Read from MCMC samples
    var_name <- line$var_name
    return(paste0(var_name, " <- subset_samps(include = '", var_name, "', samps = samps)"))
  } else if (line$block_name == "data") {
    # Read from input data
    var_name <- line$var_name
    return(paste0(var_name, " <- dat$", var_name))
  } else {
    # Initialize the variable in R
    var_name <- line$var_name
    declaration_type <- line$declaration_type
    dimensions <- line$dimensions
    
    # Handle different declaration types
    if (declaration_type == "int" || declaration_type == "real") {
      if (dimensions == "1") {
        # Scalar
        return(paste0(var_name, " <- NA"))
      } else {
        # Vector
        return(paste0(var_name, " <- rep(NA, ", dimensions, ")"))
      }
    } else if (declaration_type == "vector") {
      # Vector
      return(paste0(var_name, " <- rep(NA, ", dimensions, ")"))
    } else if (declaration_type == "matrix") {
      # Matrix
      matrix_dims <- strsplit(dimensions, ",")[[1]]
      return(paste0(var_name, " <- matrix(NA, nrow = ", matrix_dims[1], ", ncol = ", matrix_dims[2], ")"))
    } else if (declaration_type == "array") {
      # Array with more than 2 dimensions
      array_dims <- paste("dim = c(", dimensions, ")", sep = "")
      return(paste0(var_name, " <- array(NA, ", array_dims, ")"))
    } else {
      # Default case for unknown types
      return(paste0("# Unknown type for ", var_name))
    }
  }
}

interpret_lines <- function(line, dat, samps) {
  switch(line$type,
         "declaration" = interpret_declaration(line, dat, samps),
         "assignment" = interpret_assignment(line, dat, samps),
         "operation" = interpret_operation(line, dat, samps),
         "loop" = interpret_loop(line, dat, samps),
         "sampling" = interpret_sampling(line, dat, samps),
         "Unknown type")
}

# Example usage
processed_code <- lapply(parsed_lines, interpret_lines, dat, samps)

line <- parsed_lines[1,]

dcls <- parsed_lines[parsed_lines$type == "declaration",]
sapply(1:nrow(dcls), function(i) cbind(dcls$line[i], dcls$block_name[i], interpret_declaration(dcls[i,], dat, samps)))


asms <- parsed_lines[parsed_lines$type == "assignment",]
sapply(1:nrow(asms), function(i) cbind(asms$line[i], asms$block_name[i], interpret_assignment(asms[i,])))


loops <- parsed_lines[parsed_lines$type == "loop",]
sapply(1:nrow(loops), function(i) cbind(loops$line[i], loops$block_name[i], interpret_loop(loops[i,])))

sss <- parsed_lines[parsed_lines$type == "sampling",]
head(sss[,!apply(apply(sss,2,is.na), 2, all)])
sapply(1:nrow(sss), function(i){print(i); cbind(sss$line[i], sss$block_name[i], 
                                      interpret_sampling(sss[i,], dat, samps, 
                                                         sample_index = 1, 
                                                         post_pred_sim = TRUE, sim = TRUE, bound_declarations = NA))
  })
