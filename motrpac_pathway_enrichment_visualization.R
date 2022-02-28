#load packages
library(data.table)
library(MotrpacBicQC)

#functions
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

#load data
enrich <- read.table("~/data/smontgom/pw-enrich-degs-per-tissue.tsv", sep = "\t", header = T)
sapply(colnames(enrich), function(coli) length(unique(enrich[,coli])))

#snag pathway categories
url <- "https://www.genome.jp/kegg/pathway.html"
html <- rvest::read_html(url)

major_categories <- unlist(lapply(unlist(lapply(as.character(html_elements(html, "h4"))[-(1:2)], function(x) strsplit(x, "\\. ")[[1]][2])), function(x) strsplit(x, "<")[[1]][1]))

categories <- as.character(html_elements(html, "b"))[10:68]
categories <- unlist(lapply(categories, function(x) strsplit(x, ">")[[1]][2]))
categories <- unlist(lapply(categories, function(x) strsplit(x, "<")[[1]][1]))
categories <- data.frame(category_id = unlist(lapply(categories, function(x) strsplit(x, "\\.")[[1]][1])),
                         subcategory_id = unlist(lapply(categories, function(x) strsplit(x, " ")[[1]][1])),
                         subcategory_name = unlist(lapply(categories, function(x) paste0(strsplit(x, " ")[[1]][-1], collapse = " "))))
categories$subcategory_id <- unlist(lapply(categories$subcategory_id, function(x) strsplit(x, "\\.")[[1]][2]))
categories$category_name <- major_categories[as.integer(categories$category_id)]

possible_pathways <- as.character(html_elements(html, "a"))
possible_pathways <- unlist(lapply(possible_pathways, function(x) strsplit(x, ">")[[1]][2]))
possible_pathways <- unlist(lapply(possible_pathways, function(x) strsplit(x, "<")[[1]][1]))
possible_pathways <- possible_pathways[-(2:grep(pattern = "Metabolic pathways", possible_pathways)-1)]

html_text <- unlist(strsplit(as.character(html), "\n"))
pathway_locs <- as.integer(sapply(possible_pathways, function(x) grep(gsub("\\(", "\\\\(", x), html_text)[1]))
category_locs <- as.integer(sapply(categories$name, function(x) grep(x, html_text)[1]))
pathway_categories <- sapply(pathway_locs, function(x) max(which(x > (category_locs))))

pathway_categories <- cbind(pathway = possible_pathways, do.call(rbind, lapply(pathway_categories, function(x) categories[x,])))
head(pathway_categories)
if(!file.exists("~/data/smontgom/motrpac_pathway_categories.txt")){data.table::fwrite(pathway_categories, file = "~/data/smontgom/motrpac_pathway_categories.txt")}
pathway_categories <- data.table::fread(file = "~/data/smontgom/motrpac_pathway_categories.txt")

#subset the enrichment table via significance filter
pval_to_use <- "gost_adj_p_value"
pval_thresh <- 0.1
enrich_sub <- enrich[enrich[,pval_to_use] < pval_thresh,]
enrich_sub <- enrich[enrich$significant,]
n_pathways_per_tissue <- table(enrich_sub[,"tissue"])

#sanity check
nrow(enrich_sub)
unique(sapply(1:nrow(enrich_sub), function(rowi) paste0(enrich_sub[rowi,c("tissue", "term_name")], collapse = ", ")))

#get useful values
tissues <- unique(enrich_sub$tissue)
cols = list(Tissue=tissue_cols[tissues], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue <- c(cols$Tissue, ALL_TISSUES = "grey15")

#get first basic panel info
tissues_per_pathway <- as.data.table(table(enrich_sub$term_name))
enrich_sub$n_tissues_pathway_appears_in <- tissues_per_pathway$N[match(enrich_sub$term_name, tissues_per_pathway$V1)]
num_tissues_across_uniqueness_thresholds <- lapply(tissues, function(tissue) as.data.table(table(enrich_sub$n_tissues_pathway_appears_in[enrich_sub$tissue == tissue])))
names(num_tissues_across_uniqueness_thresholds) <- tissues
cumsum_tissues_across_uniqueness_thresholds <- lapply(tissues, function(tissue) cbind(n_tissues = as.numeric(num_tissues_across_uniqueness_thresholds[[tissue]]$V1),
                                                                                      cumul_sum = cumsum(num_tissues_across_uniqueness_thresholds[[tissue]]$N)))
names(cumsum_tissues_across_uniqueness_thresholds) <- tissues
set_of_allTissues <- as.data.table(table(tissues_per_pathway$N))
set_of_allTissues <- cbind(n_tissues = as.numeric(set_of_allTissues$V1),
                           cumul_sum = cumsum(set_of_allTissues$N))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### find partial correlation jaccard index thing ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


enrich_sub_dt <- as.data.table(enrich_sub)
#approach 1 -- weighted intersects, weighted inverse to n_tissues per pathway
weighted_intersects <- sapply(1:length(tissues), function(ti1) sapply(1:length(tissues), function(ti2)
  sum(1/tissues_per_pathway$N[match(intersect(enrich_sub_dt[tissue == tissues[ti1], term_name], 
      enrich_sub_dt[tissue == tissues[ti2], term_name]), tissues_per_pathway$V1)]) 
))


weighted_intersects_over_unions <- sapply(1:length(tissues), function(ti1) sapply(1:length(tissues), function(ti2)
  sum(1/tissues_per_pathway$N[match(intersect(enrich_sub_dt[tissue == tissues[ti1], term_name], 
                                              enrich_sub_dt[tissue == tissues[ti2], term_name]), tissues_per_pathway$V1)]) / 
    sum(1/tissues_per_pathway$N[match(union(enrich_sub_dt[tissue == tissues[ti1], term_name], 
                                            enrich_sub_dt[tissue == tissues[ti2], term_name]), tissues_per_pathway$V1)])
))
rownames(weighted_intersects_over_unions) <- colnames(weighted_intersects_over_unions) <- tissues

rate = 0.05
exp_dist_cols <- round(cumsum(c(1, dexp(1:100, rate = rate) / min(dexp(1:100, rate = rate)))))
heatcols <- viridis::magma(max(exp_dist_cols))[exp_dist_cols]

plot(1,1,xlim = c(0,length(tissues)), ylim = c(0,length(tissues)), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

for(rowi in 1:length(tissues)){
  text(labels = tissues[rowi], x = length(tissues) + 0.5, y = length(tissues) - rowi + 1, col = cols$Tissue[tissues[rowi]], pos = 4)
  for(colj in 1:length(tissues)){
    points(y = length(tissues) - rowi + 1, x = colj, pch = 15, cex = 4,
           col = heatcols[round(weighted_intersects_over_unions[rowi, colj] * 100) + 1])
    if(rowi == 1){
      text(labels = tissues[colj], x = colj, y = -0.15, col = cols$Tissue[tissues[colj]], pos = 1, srt = 270)
    }
  }
}

plot(1,1,xlim = c(0,length(tissues)), ylim = c(0,length(tissues)), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

for(rowi in 1:length(tissues)){
  text(labels = tissues[rowi], x = length(tissues) + 0.5, y = length(tissues) - rowi + 1, col = cols$Tissue[tissues[rowi]], pos = 4)
  for(colj in 1:length(tissues)){
    points(y = length(tissues) - rowi + 1, x = colj, pch = 15, cex = 3,
           col = heatcols[round(weighted_intersects[rowi, colj] / max(weighted_intersects) * 100) + 1])
    
    if(rowi == 1){
      text(labels = tissues[colj], x = colj, y = -0.15, col = cols$Tissue[tissues[colj]], pos = 1, srt = 270)
    }
    
  }
}

plot(cov2cor(weighted_intersects), weighted_intersects_over_unions)

#four pseudojaccards -- weighted, probability, average subtract out, average similarities from pairs-triplets-tetrads, x-entropy?

mean_3way_partial_jaccard_index <- sapply(1:length(tissues), function(ti1) sapply(1:length(tissues), function(ti2)
  mean(na.rm = T, sapply((1:length(tissues))[-c(ti1,ti2)], function(ti3)
    
    length(setdiff(intersect(enrich_sub_dt[tissue == tissues[ti1], term_name], enrich_sub_dt[tissue == tissues[ti2], term_name]), enrich_sub_dt[tissue == tissues[ti3], term_name])) / 
    length(setdiff(union(enrich_sub_dt[tissue == tissues[ti1], term_name], enrich_sub_dt[tissue == tissues[ti2], term_name]), enrich_sub_dt[tissue == tissues[ti3], term_name]))
    
    
))))
rownames(mean_3way_partial_jaccard_index) <- colnames(mean_3way_partial_jaccard_index) <- tissues

par(mar = c(4,4,4,4))
plot(1,1,xlim = c(0,length(tissues)), ylim = c(0,length(tissues)), xpd = NA,
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sorted_row_inds <- order(cmdscale(1-mean_3way_partial_jaccard_index, k = 1))
for(rowi in 1:length(tissues)){
  text(labels = tissues[rowi], x = length(tissues) + 0.5, y = length(tissues) - rowi + 1, col = cols$Tissue[tissues[rowi]], pos = 4, xpd = NA)
  for(colj in 1:length(tissues)){
    points(y = length(tissues) - rowi + 1, x = colj, pch = 15, cex = 3,
           col = heatcols[round(mean_3way_partial_jaccard_index[sorted_row_inds[rowi], sorted_row_inds[colj]] / max(mean_3way_partial_jaccard_index) * 100) + 1])
    
    if(rowi == 1){
      text(labels = tissues[sorted_row_inds[colj]], x = colj, y = -0.15, col = cols$Tissue[tissues[sorted_row_inds[colj]]], pos = 1, srt = 270, xpd = NA)
    }
    
  }
}

# probability jaccard index
enrich_sub$term_index <- as.integer(as.factor(enrich_sub$term_name))
enrich_sub_list <- lapply(tissues, function(tissue) enrich_sub[enrich_sub$tissue == tissue,c("term_index", "n_tissues_pathway_appears_in")]); names(enrich_sub_list) = tissues
pathways <- sort(unique(unlist(lapply(tissues, function(tissue) enrich_sub_list[[tissue]]$term_index))))
tissue_probabilities <- t(sapply(1:length(tissues), function(ti) rep(0, length(pathways))))
rownames(tissue_probabilities) <- tissues
for(tissue in tissues){
  tissue_probabilities[tissue,enrich_sub_list[[tissue]]$term_index] <- 1/enrich_sub_list[[tissue]]$n_tissues_pathway_appears_in
}
tissue_probability_simplices <- t(sapply(tissues, function(ti) tissue_probabilities[ti,] / sum(tissue_probabilities[ti,]) ))
prob_jaccard <- function(x,y){
  valid_summation_indices <- which(x != 0 & y != 0)
  if(length(valid_summation_indices) == 0){return(0)}
  output <- sum(sapply(valid_summation_indices, function(i) 1/sum(sapply(1:length(x), function(j) max(x[j]/x[i], y[j]/y[i])))))
  if(output == Inf){return(0)}else{return(output)}
}
probability_jaccard_matrix <- sapply(tissues, function(ti1) sapply(tissues, function(ti2) prob_jaccard(x = tissue_probability_simplices[ti1,], y = tissue_probability_simplices[ti2,])))
plot(probability_jaccard_matrix, weighted_intersects_over_unions) #identical! :o as anticipated

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### alternative designs for line graph ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cumulative_less_than <- average_contribution <- optimal_marginal_contribution_start <- optimal_marginal_contribution_end <- greedy_marginal_contribution <- matrix(0, nrow = length(tissues), ncol = length(tissues))
rownames(cumulative_less_than) <- rownames(average_contribution) <- rownames(optimal_marginal_contribution_start) <- rownames(optimal_marginal_contribution_end) <- rownames(greedy_marginal_contribution) <- tissues
average_contribution[,1] <- optimal_marginal_contribution_start[,1] <- optimal_marginal_contribution_end[,1] <- greedy_marginal_contribution[,1] <- sapply(tissues, function(ti) length(unique(enrich_sub_list[[ti]]$term_index)))

#optimal order of tissues, greedy accumulator
enrich_sub$term_index <- as.integer(as.factor(enrich_sub$term_name))
enrich_sub_list <- lapply(tissues, function(tissue) enrich_sub[enrich_sub$tissue == tissue,c("term_index", "n_tissues_pathway_appears_in")])
names(enrich_sub_list) = tissues

tissue_order <- tissues[which.max(sapply(tissues, function(tissue) length(unique(enrich_sub_list[[tissue]]$term_index))))]
for(ti in 1:(length(tissues)-1)){
  tissues_remaining <- setdiff(tissues, tissue_order)
  pathways_already_here <- unique(unlist(lapply(tissue_order, function(ti2) enrich_sub_list[[ti2]]$term_index)))
  marginal_contribution <- sapply(tissues_remaining, function(tissue) length(union(enrich_sub_list[[tissue]]$term_index, pathways_already_here)))
  tissue_order <- c(tissue_order, tissues_remaining[which.max(marginal_contribution)])
}

#optimal order of tissues, exhaustive combinatorics
sum(choose(length(tissues), 1:length(tissues))) #totally doable!
all_possible_combinations <- lapply(1:length(tissues), function(n_tissues) t(combn(x = 1:length(tissues), n_tissues)))
all_possible_sets <- lapply(1:length(tissues), function(n_tissues) lapply(1:nrow(all_possible_combinations[[n_tissues]]),
                                               function(combi) unique(unlist(lapply(1:ncol(all_possible_combinations[[n_tissues]]), 
                                               function(ti) enrich_sub_list[[all_possible_combinations[[n_tissues]][combi,ti]]]$term_index
                                               )))))
n_pathways_in_all_possible_sets <- lapply(1:length(tissues), function(n_tissues) sapply(1:nrow(all_possible_combinations[[n_tissues]]),
                                                             function(combi) length(all_possible_sets[[n_tissues]][[combi]])))
sets_with_max_pathways <- sapply(1:length(tissues), function(n_tissues) 
  tissues[all_possible_combinations[[n_tissues]][which.max(n_pathways_in_all_possible_sets[[n_tissues]]),]])
marginal_optimal_addition <- unlist(lapply(2:length(sets_with_max_pathways), function(nt) setdiff(sets_with_max_pathways[[nt]], sets_with_max_pathways[[nt-1]])))

#find optimal set containing *that* tissue
sets_with_tissue <- lapply(1:length(tissues), function(fti) lapply(1:length(tissues), function(ntiss) 
  which(sapply(1:nrow(all_possible_combinations[[ntiss]]), function(combi) fti %in% all_possible_combinations[[ntiss]][combi,]))))
optimal_sets_to <- t(sapply(1:length(tissues), function(fti) sapply(1:length(tissues), function(ntiss) 
  sets_with_tissue[[fti]][[ntiss]][which.max(n_pathways_in_all_possible_sets[[ntiss]][sets_with_tissue[[fti]][[ntiss]]])]
)))
optimal_sets_to_marginal_contribution <- t(sapply(1:length(tissues), function(fti) sapply(1:length(tissues), function(ntiss) 
  length(setdiff(enrich_sub_list[[fti]]$term_index, unique(unlist(lapply(setdiff(all_possible_combinations[[ntiss]][optimal_sets_to[fti, ntiss],], fti), function(nfti) enrich_sub_list[[nfti]]$term_index))
))))))
optimal_marginal_contribution_end[1:length(tissues), 1:length(tissues)] <- optimal_sets_to_marginal_contribution

#now iterate over this all the tissues performing this procedure to get average, greedy, and optimal marginal quantities (and best possible combination with that tissue included)

for(tissue_to_leave_out in tissues){
  
  print(tissue_to_leave_out)
  subtiss <- setdiff(tissues, tissue_to_leave_out)
  enrich_sub_list_loo <- enrich_sub_list[subtiss]
  
  #greedy
  tissue_order <- subtiss[which.max(sapply(subtiss, function(tissue) length(unique(enrich_sub_list_loo[[tissue]]$term_index))))]
  for(ti in 1:(length(subtiss)-1)){
    tissues_remaining <- setdiff(subtiss, tissue_order)
    pathways_already_here <- unique(unlist(lapply(tissue_order, function(ti2) enrich_sub_list_loo[[ti2]]$term_index)))
    marginal_contribution <- sapply(tissues_remaining, function(tissue) length(union(enrich_sub_list_loo[[tissue]]$term_index, pathways_already_here)))
    tissue_order <- c(tissue_order, tissues_remaining[which.max(marginal_contribution)])
  }
  greedy_marginal_contribution[tissue_to_leave_out,2:length(tissues)] <- sapply(1:length(tissue_order), function(ntiss) 
    length(setdiff(enrich_sub_list[[tissue_to_leave_out]]$term_index ,unique(unlist(lapply(tissue_order[1:ntiss], function(ti) enrich_sub_list_loo[[ti]]$term_index))))))
  
  #all
  all_possible_combinations <- lapply(1:length(subtiss), function(n_tissues) t(combn(x = 1:length(subtiss), n_tissues)))
  all_possible_sets <- lapply(1:length(subtiss), function(n_tissues) lapply(1:nrow(all_possible_combinations[[n_tissues]]),
                                                                            function(combi) unique(unlist(lapply(1:ncol(all_possible_combinations[[n_tissues]]), 
                                                                                                                 function(ti) enrich_sub_list_loo[[all_possible_combinations[[n_tissues]][combi,ti]]]$term_index
                                                                            )))))
  average_contribution[tissue_to_leave_out,2:length(tissues)] <- sapply(1:length(subtiss), function(n_tissues) mean(sapply(1:nrow(all_possible_combinations[[n_tissues]]),
                                                                                           function(combi) length(setdiff(enrich_sub_list[[tissue_to_leave_out]]$term_index, all_possible_sets[[n_tissues]][[combi]])))))
  
  n_pathways_in_all_possible_sets <- lapply(1:length(subtiss), function(n_tissues) sapply(1:nrow(all_possible_combinations[[n_tissues]]),
                                                                                          function(combi) length(all_possible_sets[[n_tissues]][[combi]])))
  
  #optimal start
  sets_with_max_pathways <- sapply(1:length(subtiss), function(n_tissues) 
    subtiss[all_possible_combinations[[n_tissues]][which.max(n_pathways_in_all_possible_sets[[n_tissues]]),]])
  optimal_marginal_contribution_start[tissue_to_leave_out,2:length(tissues)] <- sapply(sets_with_max_pathways, function(seti) length(setdiff(enrich_sub_list[[tissue_to_leave_out]]$term_index, 
                                                               unlist(lapply(seti, function(ti) enrich_sub_list[[ti]]$term_index)))))
  
}

for(tissue in tissues){
  cumsumsub <- as.data.frame(cumsum_tissues_across_uniqueness_thresholds[[tissue]]  )
  cumulative_less_than[tissue,cumsumsub$n_tissues] <- cumsumsub$cumul_sum
  for(i in 2:length(cumulative_less_than[tissue,])){
    if(cumulative_less_than[tissue,i] == 0){
      cumulative_less_than[tissue,i] <- cumulative_less_than[tissue,i-1]
    }
  }
  
}

all_line_options <- list(greedy = (greedy_marginal_contribution), 
                              optimal_end = (optimal_marginal_contribution_end), 
                              optimal_start = (optimal_marginal_contribution_start),
                              average = (average_contribution))
all_line_options_mat <- cbind(greedy = as.vector(greedy_marginal_contribution[,-1]), 
                           optimal_end = as.vector(optimal_marginal_contribution_end[,-1]), 
                           optimal_start = as.vector(optimal_marginal_contribution_start[,-1]),
                           average = as.vector(average_contribution[,-1]))
pairs(all_line_options_mat)
cor(all_line_options_mat)


tissue_specific_corrs <- sapply(all_line_options, function(summary_matrix) sapply(tissues, function(tissue) cor(cumulative_less_than[tissue,], summary_matrix[tissue,])))
for(i in 1:ncol(tissue_specific_corrs)){
  hist(add = ifelse(i==1,F,T), col = adjustcolor(i, 0.5), tissue_specific_corrs[,i], 
       xlim = c(-1,1), ylim = c(0,15), breaks = -10:10/10)
}


#### start plotting ####
#four curve plots -- cumulative x-or-fewer, average marginal, greedy marginal, optimal marginal
#four pseudojaccards -- weighted, probability, average subtract out, average similarities from pairs-triplets-tetrads, x-entropy?

#make frame
grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/marginal_cumulative_contribution.pdf"), 
                     width = 600 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mar = c(3,3,3,5), mfrow = c(1,1), xpd = NA)
plot(1,1,xlim = c(0,max(enrich_sub$n_tissues_pathway_appears_in)), ylim = c(0,max(set_of_allTissues)), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#calculate more useful values
xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(max(enrich_sub$n_tissues_pathway_appears_in) + 0.25, 
                                                          c(n_pathways_per_tissue, ALL_TISSUES = max(set_of_allTissues[,2])) + 
                                                            rnorm(length(n_pathways_per_tissue) + 1, 0, 1E-3)),
                                           rep.fact = 50, adj.max = 0.1)
text(x = max(enrich_sub$n_tissues_pathway_appears_in) / 2 + 0.5, y = max(set_of_allTissues) * 1.05, labels = "Marginal Cumulative Contribution", cex = 2)

#draw horiz axes
segments(x0 = 0, x1 = max(enrich_sub$n_tissues_pathway_appears_in), y0 = 0, y1 = 0, lwd = 2)
segments(x0 = 1:max(enrich_sub$n_tissues_pathway_appears_in), x1 = 1:max(enrich_sub$n_tissues_pathway_appears_in), y0 = 0, y1 = -2, lwd = 2)
text(y = -1, x = 1:max(enrich_sub$n_tissues_pathway_appears_in), labels = 1:max(enrich_sub$n_tissues_pathway_appears_in), pos = 1)
text(labels = "Number of Tissues", y = -10, x = max(enrich_sub$n_tissues_pathway_appears_in) / 2, pos = 1, cex = 2)

#draw vert axis
segments(x0 = 0, x1 = 0, y0 = 0, y1 = max(set_of_allTissues), lwd = 2)
ylocs_vertaxis <- seq(1, max(set_of_allTissues), by = round(diff(seq(1, max(set_of_allTissues), length.out = 10))[1]))
segments(y0 = ylocs_vertaxis, y1 = ylocs_vertaxis, x0 = 0, x1 = -0.15, lwd = 2)
text(x = -0.1, y = ylocs_vertaxis, labels = ylocs_vertaxis, pos = 2)
text(labels = "Number of Pathways in Tissue that Appear in at Most Some Number of Tissues", x = -1.1, y = max(set_of_allTissues) / 2, srt = 90)

for(tissue in c(tissues, "ALL_TISSUES")){
  if(tissue == "ALL_TISSUES"){
    cumsumsub <- set_of_allTissues
  } else {
    cumsumsub <- cumsum_tissues_across_uniqueness_thresholds[[tissue]]  
  }
  cumsumsub <- rbind(c(0,0), cumsumsub)
  if(max(cumsumsub[,"n_tissues"]) < max(enrich_sub$n_tissues_pathway_appears_in)){
    cumsumsub <- rbind(cumsumsub, c(max(enrich_sub$n_tissues_pathway_appears_in), cumsumsub[nrow(cumsumsub), 2]))
  }
  lines(x = cumsumsub[,"n_tissues"], y = cumsumsub[,"cumul_sum"],
        col = cols$Tissue[tissue], lwd = 2)
  text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = 0.75,
       y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], labels = tissue, col = cols$Tissue[tissue], pos = 4)
  segments(x0 = max(cumsumsub[,"n_tissues"]), y0 = max(cumsumsub[,"cumul_sum"]), 
           x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.125, 
           y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
           col = cols$Tissue[tissue], lty = 3)
}

dev.off()

#### plot weighted coloc plot ####
rate = 0.05
exp_dist_cols <- round(cumsum(c(1, dexp(1:100, rate = rate) / min(dexp(1:100, rate = rate)))))
heatcols <- viridis::magma(max(exp_dist_cols))[exp_dist_cols]

for(ji in 1:2){
grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/", c("weighted-probability-jaccard", "mean-3way-jaccard")[ji],".pdf"), 
                     width = 1800 / 72, height = 600 / 72, family="Arial Unicode MS")
layout(cbind(matrix(1,6,6), matrix(c(sapply(list(2:7, 8:13, 14:19), function(n) rep(c(sapply(n, function(x) rep(x,2))),2))), nrow = 6, ncol = 12, byrow = T)))

jaccard_matrix_to_use <- list(weighted_intersects_over_unions, mean_3way_partial_jaccard_index)[[ji]]
par(mar = c(4,3.5,4,5), xpd = NA)
plot(1,1,xlim = c(0,length(tissues)), ylim = c(0,length(tissues)), xpd = NA,
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sorted_row_inds <- order(cmdscale(1-jaccard_matrix_to_use, k = 1))
for(rowi in 1:length(tissues)){
  text(labels = tissues[sorted_row_inds[rowi]], x = length(tissues) + 0.5, y = length(tissues) - rowi + 1, 
       col = cols$Tissue[tissues[sorted_row_inds[rowi]]], pos = 4, xpd = NA, cex = 1.5, font = 2)
  for(colj in 1:length(tissues)){
    points(y = length(tissues) - rowi + 1, x = colj, pch = 15, cex = 7,
           col = heatcols[round(jaccard_matrix_to_use[sorted_row_inds[rowi], sorted_row_inds[colj]] / max(jaccard_matrix_to_use) * 100) + 1])
    
    if(rowi == 1){
      text(labels = tissues[sorted_row_inds[colj]], x = colj-0.225, y = 0.325, col = cols$Tissue[tissues[sorted_row_inds[colj]]], pos = 4, srt = 270, xpd = NA, cex = 1.5, font = 3)
    }
    
  }
}

xl = -0.75; xr = 0; yb = 10; yt = 17.5
rect(xleft = xl, xright = xr, col = heatcols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)], yt= seq(yb, yt, length.out = length(heatcols)+1)[-1])
rect(xleft = xl, xright = xr, ybottom = yb, ytop = yt)
text(labels = 0:10/10, x = xl, pos = 2, y = seq(yb, yt, length.out = 11), cex = 1.25)
text(c("Weighted / Inverse Probability Jaccard Similarity", "Expected Partial 3-Way Jaccard Similarity")[ji], x = length(tissues) / 2 + 0.5, y = length(tissues) + 0.5, pos = 3, cex = 2.5, font = 2)

#plot eigenvector loadings and scree plot
wiou_eigen <- eigen(jaccard_matrix_to_use)

par(mar = c(3.5,0.5,0.5,0.5), xpd = NA)
for(i in 1:length(tissues)){
  plot(100,100,xlim = c(0,length(tissues)), ylim = c(0,max(abs(wiou_eigen$vectors[,i]))), 
       col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  text(paste0("PC ", i, " Loadings"), x = length(tissues) / 2, y = 0.95*max(abs(wiou_eigen$vectors[,i])), cex = 1.25)
  rect(xleft = 1:length(tissues)-1/2, xright = 1:length(tissues)+1/2, 
       ybottom = 0, ytop = abs(wiou_eigen$vectors[,i])[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], 
       col = cols$Tissue[order(abs(wiou_eigen$vectors[,i]), decreasing = T)])
  text(labels = tissues[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], x = 1:length(tissues) - 0.75, 
       y = -max(abs(wiou_eigen$vectors[,i]))*0.01, col = cols$Tissue[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], 
       pos = 4, srt = 270)
}

par(mar = c(5,3,1,2))
plot(cex.lab = 1.5, x = 1:length(tissues), y = wiou_eigen$values, type = "l", xlab = "ordered eigenvalue indices", ylab = "eigenvalue", lwd = 2)
dev.off()
}

pdftools::pdf_combine(paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/", c("weighted-probability-jaccard", "mean-3way-jaccard"),".pdf"), 
                      output = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/jaccard_alternatives.pdf"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### iterate over alternative line plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for(li in 1:length(all_line_options)){

line_matrix <- all_line_options[[li]]
#make frame
grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/", c("greedy", "optimal-ending", "optimal-starting", "average")[li], "_marginal_contribution.pdf"), 
                     width = 600 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mar = c(3,4,3,1), mfrow = c(1,1), xpd = NA)
plot(1,1,xlim = c(-2,17), ylim = c(0,max(line_matrix)), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#calculate more useful values
xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(0.5, line_matrix[,1] + rnorm(17, 0, 1E-3)),
                                           rep.fact = 20, adj.max = 0.1, )
text(x = 17 / 2 - 0.5, y = max(line_matrix) * 1.075, 
     labels = paste0(c("Greedy", "Optimal-Ending", "Optimal-Starting", "Average")[li], " Marginal Contribution"), cex = 2)

#draw horiz axes
segments(x0 = -2, x1 = 17, y0 = 0, y1 = 0, lwd = 2)
segments(x0 = 1:17, x1 = 1:17, y0 = 0, y1 = -2, lwd = 2)
text(y = -1, x = 1:17, labels = 1:17, pos = 1)
text(labels = "Number of Tissues (Ending Count)", y = -7, x = max(enrich_sub$n_tissues_pathway_appears_in) / 2 + 3, pos = 1, cex = 1.75)

#draw vert axis
segments(x0 = -2, x1 = -2, y0 = 0, y1 = max(line_matrix), lwd = 2)
ylocs_vertaxis <- seq(1, max(line_matrix), by = round(diff(seq(1, max(line_matrix), length.out = 10))[1]))
segments(y0 = ylocs_vertaxis, y1 = ylocs_vertaxis, x0 = -2, x1 = -2.15, lwd = 2)
text(x = -2.1, y = ylocs_vertaxis, labels = ylocs_vertaxis, pos = 2)
text(labels = "Number of Pathways Contributed at Given Tissue Count", x = -4, y = max(line_matrix) / 2, srt = 90, cex = 1.5)

for(tissue in tissues){
  
  cumsumsub <- line_matrix[tissue,]
  # cumsumsub <- rbind(c(0,0), cumsumsub)
  lines(x = 1:length(cumsumsub), y = cumsumsub,
        col = cols$Tissue[tissue], lwd = 2)
  text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = 1,
       y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], labels = tissue, col = cols$Tissue[tissue], pos = 2)
  segments(x0 = 1, y0 = cumsumsub[1], 
           x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 0.15, 
           y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
           col = cols$Tissue[tissue], lty = 3)
}

dev.off()
}

pdftools::pdf_combine(c(paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/marginal_cumulative_contribution.pdf"),
                        paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/", c("greedy", "optimal-ending", "optimal-starting", "average"), "_marginal_contribution.pdf")), 
                      output = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/marginal_contribution_alternatives.pdf"))


#~~~~~~~~~~~~~~~#
#### addenda ####
#~~~~~~~~~~~~~~~#

#coerce to matrix of 0s and 1s
str(enrich_sub)
es <- as.data.table(enrich_sub)

pathways <- unique(enrich_sub$term_name)
tissues <- unique(enrich_sub$tissue)
tiss <- tissues[2]
sparse_vector <- function(inds, length){
  x <- rep(0, max(c(length, inds)))
  x[inds] <- 1
  x
}
sparse_mat <- t(sapply(tissues, function(tiss) sparse_vector(match(c(es[tissue == tiss,"term_name"])$term_name, pathways), length(pathways))))
colnames(sparse_mat) <- pathways
sparse_mat_char <- t(sapply(tissues, function(tiss) as.character(sparse_vector(match(c(es[tissue == tiss,"term_name"])$term_name, pathways), length(pathways)))))
colnames(sparse_mat_char) <- pathways

lPCA_cv <- logisticPCA::cv.lpca(sparse_mat, ks = 5, ms = 1:10)
plot(lPCA_cv)
lPCA <- logisticPCA::logisticPCA(sparse_mat, k = nrow(sparse_mat), m = which.min(lPCA_cv))
str(lPCA)
lPCA$PCs

polychoric_corrMat <- as.matrix(psych::polychoric((sparse_mat))$rho)
eigen_pcm <- eigen(polychoric_corrMat)
tetra_propvar <- round((eigen_pcm$values / sum(eigen_pcm$values))*100)
tetra_coords <- sparse_mat %*% eigen_pcm$vectors[,1:5]
MCAout <- FactoMineR::MCA(X = as.data.frame(sparse_mat_char))

panel.1Dnames <- function(x,tissue_names,tissue_colors,...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(-0.05, 1.05, -0.05, 1.05))
  locs <- x
  locs <- locs - min(locs)
  locs <- locs / max(locs)
  text(x = 0.5, y = locs, labels = tissue_names, cex = 1, font = 2, col = tissue_colors)
}

my.text.panel <- function(labels) {
  function(x, y, lbl, ...) {
    if (lbl %in% names(labels)) lbl <- labels[[lbl]]
    text(0.05, 0.5, lbl, srt = 90, cex = 1)
  }
}

var_lab <- paste0(round(MCAout$eig[,2]), "% of Variance")
names(var_lab) <- letters[1:5]

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/MCA_Pairs_Plot.pdf"), 
                     width = 650 / 72, height = 550 / 72, family="Arial Unicode MS")

pairs(MCAout$ind$coord, diag.panel = panel.1Dnames, tissue_names = rownames(MCAout$ind$coord), labels = letters[1:5],
      tissue_colors = cols$Tissue[match(rownames(MCAout$ind$coord), names(cols$Tissue))], 
      text.panel = my.text.panel(var_lab), pch = 19, cex = 1.5, 
      col = adjustcolor(cols$Tissue[match(rownames(MCAout$ind$coord), names(cols$Tissue))], 0.75), oma =  c(2,8,2,2))
xl = 0.0625;xr = 0.075; yb = 0.825; yt = 1
rect(xpd = NA, col = cols$Tissue[-length(cols$Tissue)], border = NA,
     xleft = rep(xl, length(cols$Tissue)-1),
     xright = rep(xr, length(cols$Tissue)-1),
     ybottom = seq(yb, yt, length.out = length(cols$Tissue)-1) - (yt-yb)/(length(cols$Tissue) - 2)/2,
     ytop = seq(yb, yt, length.out = length(cols$Tissue)-1) + (yt-yb)/(length(cols$Tissue) - 2)/2)
text(labels = names(cols$Tissue[-length(cols$Tissue)]), seq(yb, yt, length.out = length(cols$Tissue)-1), pos = 2, x = (xl + xr) / 2, xpd = NA, cex = 0.5)
text(0, 0.5, srt = 90,  xpd = NA, labels = "Multiple Coordinates Analysis (MCA)", cex = 1.2, font = 2)

dev.off()

var_lab <- paste0(tetra_propvar[1:5], "% of Variance")
names(var_lab) <- letters[1:5]
pairs(tetra_coords, diag.panel = panel.1Dnames, tissue_names = rownames(tetra_coords), labels = letters[1:5],
      tissue_colors = cols$Tissue[match(rownames(tetra_coords), names(cols$Tissue))], 
      text.panel = my.text.panel(var_lab), pch = 19, cex = 1.5, 
      col = adjustcolor(cols$Tissue[match(rownames(tetra_coords), names(cols$Tissue))], 0.75), oma =  c(2,8,2,2))
xl = 0.0625;xr = 0.075; yb = 0.825; yt = 1
rect(xpd = NA, col = cols$Tissue[-length(cols$Tissue)], border = NA,
     xleft = rep(xl, length(cols$Tissue)-1),
     xright = rep(xr, length(cols$Tissue)-1),
     ybottom = seq(yb, yt, length.out = length(cols$Tissue)-1) - (yt-yb)/(length(cols$Tissue) - 2)/2,
     ytop = seq(yb, yt, length.out = length(cols$Tissue)-1) + (yt-yb)/(length(cols$Tissue) - 2)/2)
text(labels = names(cols$Tissue[-length(cols$Tissue)]), seq(yb, yt, length.out = length(cols$Tissue)-1), pos = 2, x = (xl + xr) / 2, xpd = NA, cex = 0.5)
text(0, 0.5, srt = 90,  xpd = NA, labels = "Tetrachoric PCA", cex = 1.2, font = 2)


lPCA_coords <- lPCA$PCs[,1:5]
var_lab <- paste0("Axis ", 1:5)
names(var_lab) <- letters[1:5]
pairs(lPCA_coords, diag.panel = panel.1Dnames, tissue_names = rownames(lPCA_coords), labels = letters[1:5],
      tissue_colors = cols$Tissue[match(rownames(lPCA_coords), names(cols$Tissue))], 
      text.panel = my.text.panel(var_lab), pch = 19, cex = 1.5, 
      col = adjustcolor(cols$Tissue[match(rownames(tetra_coords), names(cols$Tissue))], 0.75), oma =  c(2,8,2,2))
xl = 0.0625;xr = 0.075; yb = 0.825; yt = 1
rect(xpd = NA, col = cols$Tissue[-length(cols$Tissue)], border = NA,
     xleft = rep(xl, length(cols$Tissue)-1),
     xright = rep(xr, length(cols$Tissue)-1),
     ybottom = seq(yb, yt, length.out = length(cols$Tissue)-1) - (yt-yb)/(length(cols$Tissue) - 2)/2,
     ytop = seq(yb, yt, length.out = length(cols$Tissue)-1) + (yt-yb)/(length(cols$Tissue) - 2)/2)
text(labels = names(cols$Tissue[-length(cols$Tissue)]), seq(yb, yt, length.out = length(cols$Tissue)-1), pos = 2, x = (xl + xr) / 2, xpd = NA, cex = 0.5)
text(0, 0.5, srt = 90,  xpd = NA, labels = "logistic PCA", cex = 1.2, font = 2)

#look at loadings
MCAout_zeros <- MCAout$var$contrib[grep(rownames(MCAout$var$contrib), pattern = "_0"),]
MCAout_ones <- MCAout$var$contrib[grep(rownames(MCAout$var$contrib), pattern = "_1"),]

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/MCA_Axis1.pdf"), 
                     width = 1000 / 72, height = 400 / 72, family="Arial Unicode MS")
factoextra::fviz_contrib(MCAout, choice = "var", axes = 1, top = 80, xtickslab.rt = 75)
dev.off()

# names(MCAout_ones[,1])[order(MCAout_ones[,1], decreasing = T)]

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/MCA_Axis2.pdf"), 
                     width = 1000 / 72, height = 400 / 72, family="Arial Unicode MS")
factoextra::fviz_contrib(MCAout, choice = "var", axes = 2, top = 80, xtickslab.rt = 76)
dev.off()

# names(MCAout_ones[,2])[order(MCAout_ones[,2], decreasing = T)]

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/MCA_Axis3.pdf"), 
                     width = 1000 / 72, height = 400 / 72, family="Arial Unicode MS")
factoextra::fviz_contrib(MCAout, choice = "var", axes = 3, top = 80, xtickslab.rt = 76)
dev.off()

# names(MCAout_ones[,3])[order(MCAout_ones[,3], decreasing = T)]

pathway_categories$pathway[grep(pathway_categories$pathway, pattern = "Chemical carcinogenesis")] <- "Chemical carcinogenesis"
contrib <- data.frame(c1 = MCAout$var$contrib[,1], c2 = MCAout$var$contrib[,2], c3 = MCAout$var$contrib[,3])
contrib$pathway <- rownames(contrib)
contrib$pathway <- gsub("_1", "", contrib$pathway)
contrib$pathway <- gsub("_0", "", contrib$pathway)
contrib <- cbind(contrib, pathway_categories[match(contrib$pathway, pathway_categories$pathway),c("category_name", "subcategory_name")])
unique_cats <- unique(contrib$category_name)
unique_subcats <- unique(contrib$subcategory_name)

#dimension 1
library(wordcloud)
round(sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c1"]))), T), 2)
head(round(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c1"])), T), 2))

dim1_catfreqs <- sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c1"]))), T)
wordcloud(words = names(dim1_catfreqs), freq = dim1_catfreqs, min.freq = 0.01,
          max.words=200, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))

dim1_subcatfreqs <- sort(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c1"]))), T)
wordcloud(words = names(dim1_subcatfreqs), freq = dim1_subcatfreqs, min.freq = 0.01,
          max.words=200, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))


#dimension 2
round(sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c2"]))), T), 2)
head(round(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c2"])), T), 2))

dim2_catfreqs <- sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c2"]))), T)
wordcloud(words = names(dim2_catfreqs), freq = dim2_catfreqs, min.freq = 0.01,
          max.words=200, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))

dim2_subcatfreqs <- sort(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c2"]))), T)
wordcloud(words = names(dim2_subcatfreqs), freq = dim2_subcatfreqs, min.freq = 0.01,
          max.words=200, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))


#dimension 3
round(sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c3"]))), T), 2)
head(round(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c3"])), T), 2))


dim3_catfreqs <- sort(sort(sapply(unique_cats, function(categ) sum(contrib[contrib$category_name == categ,"c3"]))), T)
wordcloud(words = names(dim3_catfreqs), freq = dim3_catfreqs, min.freq = 0.01,
          max.words=300, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))

dim3_subcatfreqs <- sort(sort(sapply(unique_subcats, function(subcateg) sum(contrib[contrib$subcategory_name == subcateg,"c3"]))), T)
wordcloud(words = names(dim3_subcatfreqs), freq = dim3_subcatfreqs, min.freq = 0.01,
          max.words=300, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))



sort(table(pathway_categories$category_name), T)
head(sort(table(pathway_categories$subcategory_name), T))

wordcloud(words = names(table(pathway_categories$category_name)), freq = table(pathway_categories$category_name), min.freq = 0.01,
          max.words=300, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))
wordcloud(words = names(table(pathway_categories$subcategory_name)), freq = table(pathway_categories$subcategory_name), min.freq = 0.01,
          max.words=300, random.order=FALSE, rot.per=0,
          colors=brewer.pal(8, "Dark2"))


#### make quick draft figure ####


grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/tissue_pathway_enrichment/enrichment_figure_draft.pdf"), 
                     width = 1100 / 72, height = 550 / 72, family="Arial Unicode MS")

layout(matrix(1:2, nrow=1))

#average marginal contribution plot
  li = which(names(all_line_options) == "average")
  line_matrix <- all_line_options[[li]]
  #make frame
  par(mar = c(3,4,3,0), xpd = NA)
  plot(1,1,xlim = c(-2,17), ylim = c(0,max(line_matrix)), 
       col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  
  #label figure
  fig_label(text = "a)", region = "plot", cex = 3, shrinkX = 2, shrinkY = 1.075)
  
  #calculate more useful values
  xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(0.5, line_matrix[,1] + rnorm(17, 0, 1E-3)),
                                             rep.fact = 20, adj.max = 0.1, )
  text(x = 17 / 2 - 0.5, y = max(line_matrix) * 1.075, 
       labels = paste0(c("Greedy", "Optimal-Ending", "Optimal-Starting", "Average")[li], " Marginal Contribution"), cex = 2)
  
  #draw horiz axes
  segments(x0 = -2, x1 = 17, y0 = 0, y1 = 0, lwd = 2)
  segments(x0 = 1:17, x1 = 1:17, y0 = 0, y1 = -2, lwd = 2)
  text(y = -1, x = 1:17, labels = 1:17, pos = 1)
  text(labels = "Number of Tissues (Ending Count)", y = -7, x = max(enrich_sub$n_tissues_pathway_appears_in) / 2 + 3, pos = 1, cex = 1.75)
  
  #draw vert axis
  segments(x0 = -2, x1 = -2, y0 = 0, y1 = max(line_matrix), lwd = 2)
  ylocs_vertaxis <- seq(1, max(line_matrix), by = round(diff(seq(1, max(line_matrix), length.out = 10))[1]))
  segments(y0 = ylocs_vertaxis, y1 = ylocs_vertaxis, x0 = -2, x1 = -2.15, lwd = 2)
  text(x = -2.1, y = ylocs_vertaxis, labels = ylocs_vertaxis, pos = 2)
  text(labels = "Number of Pathways Contributed at Given Tissue Count", x = -4, y = max(line_matrix) / 2, srt = 90, cex = 1.5)
  
  for(tissue in tissues){
    
    cumsumsub <- line_matrix[tissue,]
    # cumsumsub <- rbind(c(0,0), cumsumsub)
    lines(x = 1:length(cumsumsub), y = cumsumsub,
          col = cols$Tissue[tissue], lwd = 2)
    text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = 1,
         y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], labels = tissue, col = cols$Tissue[tissue], pos = 2)
    segments(x0 = 1, y0 = cumsumsub[1], 
             x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] - 0.15, 
             y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
             col = cols$Tissue[tissue], lty = 3)
  }

#MCA pairs plot
  xs = 100
  ys = 100
  nd <- ncol(MCAout$ind$coord)
  nticks <- 5
  var_lab <- paste0(tetra_propvar[1:nd], "% of Variance")
  par(mar = c(1,1,3,1), xpd = NA)
  plot(-1E9,-1E9,xlim = c(0,xs), ylim = c(0,ys), xaxt = "n", yaxt = "n", col = 0, xlab= "", ylab = "", frame.plot = F)
  #label figure
  fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 100, shrinkY = 1.07)
  
  for(ri in 1:nd){
    
    rlim <- seq(0,ys,length.out = nd+1)[c(ri,ri+1)]
    rlim <- rlim + c(diff(rlim)/20,-diff(rlim)/20)# * par("din")[1] / par("din")[2]
    rlim_pts <- rlim + c(diff(rlim)/20,-diff(rlim)/20)
    
    for(ci in 1:nd){
      clim <- seq(0,xs,length.out = nd+1)[c(ci,ci+1)]
      clim <- clim + c(diff(clim)/20,-diff(clim)/20)
      clim_pts <- clim + c(diff(clim)/20,-diff(clim)/20)
      
      #legend
      if(ri == nd & ci == 1){
        xl = clim[1]- xs / 45; xr = clim[1] - xs / 100; yb = rlim[1]; yt = rlim[2]
        rect(xpd = NA, col = cols$Tissue[-length(cols$Tissue)], border = NA,
             xleft = rep(xl, length(cols$Tissue)-1),
             xright = rep(xr, length(cols$Tissue)-1),
             ybottom = seq(yb, yt, length.out = length(cols$Tissue)-1) - (yt-yb)/(length(cols$Tissue) - 2)/2,
             ytop = seq(yb, yt, length.out = length(cols$Tissue)-1) + (yt-yb)/(length(cols$Tissue) - 2)/2)
        text(labels = names(cols$Tissue[-length(cols$Tissue)]), seq(yb, yt, length.out = length(cols$Tissue)-1), pos = 2, x = (xl + xr) / 2, xpd = NA, cex = 0.4)
        text(xs / 2, ys * 1.05, xpd = NA, labels = "Multiple Coordinates Analysis (MCA)", cex = 2, font = 2, pos = 3)
      }
      
      #tick marks
      if(ri == 1 & (ci %% 2) == 1){
        orig_xlim = range(MCAout$ind$coord[,ci])
        tick_vals <- seq(orig_xlim[1], orig_xlim[2], length.out = nticks)
        tick_locs <- (tick_vals - min(tick_vals)) / diff(orig_xlim) * diff(clim_pts) + clim_pts[1]
        ticks_to_label <- (1:nticks)[(1:nticks) %% 2 == 1]
        segments(x0 = tick_locs, x1 = tick_locs, y0 = rlim[1], y1 = rlim[1] - ys / 100 * par("din")[1] / par("din")[2])
        text(x = tick_locs[ticks_to_label], y = rlim[1] - ys / 100 * par("din")[1] / par("din")[2], labels = round(tick_vals[ticks_to_label], 1), 
             pos = 1, cex = 0.75, xpd = NA)
      }
      
      if(ri == nd & (ci %% 2) == 0){
        orig_xlim = range(MCAout$ind$coord[,ci])
        tick_vals <- seq(orig_xlim[1], orig_xlim[2], length.out = nticks)
        tick_locs <- (tick_vals - min(tick_vals)) / diff(orig_xlim) * diff(clim_pts) + clim_pts[1]
        ticks_to_label <- (1:nticks)[(1:nticks) %% 2 == 1]
        segments(x0 = tick_locs, x1 = tick_locs, y0 = rlim[2], y1 = rlim[2] + ys / 100 * par("din")[1] / par("din")[2])
        text(x = tick_locs[ticks_to_label], y = rlim[2] + ys / 100 * par("din")[1] / par("din")[2], labels = round(tick_vals[ticks_to_label], 1), 
             pos = 3, cex = 0.75, xpd = NA)  
      }
      
      if(ci == 1 & (ri %% 2) == 0){
        orig_ylim = range(MCAout$ind$coord[,(nd-ri+1)])
        tick_vals <- seq(orig_ylim[1], orig_ylim[2], length.out = nticks)
        tick_locs <- (tick_vals - min(tick_vals)) / diff(orig_ylim) * diff(rlim_pts) + rlim_pts[1]
        ticks_to_label <- (1:nticks)[(1:nticks) %% 2 == 1]
        segments(y0 = tick_locs, y1 = tick_locs, x0 = clim[1], x1 = clim[1] - xs / 100)
        text(y = tick_locs[ticks_to_label], x = clim[1] - xs / 100, labels = round(tick_vals[ticks_to_label], 1), pos = 2, cex = 0.75, xpd = NA)
      }      
      
      if(ci == nd & (ri %% 2) == 1){
        orig_ylim = range(MCAout$ind$coord[,(nd-ri+1)])
        tick_vals <- seq(orig_ylim[1], orig_ylim[2], length.out = nticks)
        tick_locs <- (tick_vals - min(tick_vals)) / diff(orig_ylim) * diff(rlim_pts) + rlim_pts[1]
        ticks_to_label <- (1:nticks)[(1:nticks) %% 2 == 1]
        segments(y0 = tick_locs, y1 = tick_locs, x0 = clim[2], x1 = clim[2] + xs / 100)
        text(y = tick_locs[ticks_to_label], x = clim[2] + xs / 100, labels = round(tick_vals[ticks_to_label], 1), pos = 4, cex = 0.75, xpd = NA)
      }
      
      #diag column text
      if((nd-ri+1) == ci){
        orig_xlim = range(MCAout$ind$coord[,ci])
        orig_ylim = range(MCAout$ind$coord[,(nd-ri+1)])
        
        text(x = mean(clim_pts), labels = rownames(MCAout$ind$coord), cex = 0.75,
               y = (MCAout$ind$coord[,(nd-ri+1)] - min(MCAout$ind$coord[,(nd-ri+1)])) / diff(orig_ylim) * diff(rlim_pts) + rlim_pts[1], xpd= NA,
               col = adjustcolor(cols$Tissue[match(rownames(MCAout$ind$coord), names(cols$Tissue))], 1))
        rect(xleft = clim[1], xright = clim[2], ybottom = rlim[1], ytop = rlim[2])
        text(labels = var_lab[(nd-ci+1)], x = rlim_pts[1] + diff(rlim_pts) / 50, y = mean(clim), srt = 90, cex = 0.675, xpd = NA)
        
      } else { #off-diag cells
        orig_xlim = range(MCAout$ind$coord[,ci])
        orig_ylim = range(MCAout$ind$coord[,(nd-ri+1)])
        
        points((MCAout$ind$coord[,ci] - min(MCAout$ind$coord[,ci])) / diff(orig_xlim) * diff(clim_pts) + clim_pts[1], 
               (MCAout$ind$coord[,(nd-ri+1)] - min(MCAout$ind$coord[,(nd-ri+1)])) / diff(orig_ylim) * diff(rlim_pts) + rlim_pts[1], xpd= NA,
               pch = 19, col = adjustcolor(cols$Tissue[match(rownames(MCAout$ind$coord), names(cols$Tissue))], 0.75))
        rect(xleft = clim[1], xright = clim[2], ybottom = rlim[1], ytop = rlim[2])
      }
    }
  }
  
dev.off()