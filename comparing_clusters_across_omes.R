#libraries
library(foreach)
library(doParallel)

if(!exists("cl")){
  cl <- makeCluster(15, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

library(MotrpacBicQC)
library(data.table)

#functions
invlogit <- function(x){return(1 / (1+exp(-x)))}
invlogit_deriv <- function(x){return(exp(-x)/((1+exp(-x))^2))}
logit <- function(p){log(p/(1-p))}
text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = cols[char_i], cex = cex, ...)
  }
}
find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 0.95, srt_height_rescaling_ratio = 1.4){
  
  strwidths <- strwidth(txt)
  strheight <- strheight(txt[1]) * srt_height_rescaling_ratio 
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.5)
  }
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
                               str_width_lefter_start_ratio = 0.01, rect_rescaling_ratio = 0.95, col = "black", multicolor_words = F, cols_list, ...){
  
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(txt) * optimal_word_placement_inf$cex
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
             string = word_to_write, pos = 4, cols = cols_list[[words_written]])
      } else {
        text(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
           labels = word_to_write, pos = 4, col = col)
      }
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

#load in data
if(!exists("cluster_membership")){
  load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
}

if(!exists("rna_dea")){
  load('~/data/smontgom/rna_dea_20210114.RData')
  load('~/data/smontgom/dea/transcript_rna_seq_20210126.RData')
  old_rna_dea <- rna_dea
  rna_dea$training_dea <- as.data.table(transcript_rna_seq$training_dea)
  rna_dea$timewise_dea <- as.data.table(transcript_rna_seq$timewise_dea)
  new_rna_dea <- rna_dea
  transcriptome <- as.data.table(transcript_rna_seq$timewise_dea)
}

if(!exists("f2gm")){
  f2gm <- fread("~/data/smontgom/feature-to-gene-map.tsv")
}

if(!exists("selected_zs_finfo")){
  load("~/data/smontgom/merged_zscores_0.1_byfdr_20201202.RData")
  table(selected_zs_finfo$ome)
}

# load("~/data/smontgom/prot_ac_20210126.RData")
# load("~/data/smontgom/prot_ph_20210126.RData")
# str(prot_ac)
# str(prot_ph)

if(!exists("prot_pr")){
  load("~/data/smontgom/prot_pr_20210126.RData")
  str(prot_pr)
  proteome <- prot_pr$timewise_dea
}

##~~ combine proteome and transcriptome ~~##
proteome$gene_id <- f2gm$gene_symbol[match(proteome$feature_ID, f2gm$feature_ID)]
transcriptome$gene_id <- f2gm$gene_symbol[match(transcriptome$feature_ID, f2gm$feature_ID)]
transcriptome <- transcriptome[!is.na(transcriptome$gene_id),]
proteome <- proteome[!is.na(proteome$gene_id),]
transcriptome$gene.sex.tissue.time <- paste0(transcriptome$gene_id, "-", transcriptome$sex, "-", transcriptome$tissue, "-", transcriptome$comparison_group)
proteome$gene.sex.tissue.time <- paste0(proteome$gene_id, "-", proteome$sex, "-", proteome$tissue, "-", proteome$comparison_group)
sum(!is.na(match(transcriptome$gene.sex.tissue.time, proteome$gene.sex.tissue.time)))
transcriptome$logFC_proteome <- proteome$logFC[match(transcriptome$gene.sex.tissue.time, proteome$gene.sex.tissue.time)]
transcriptome$tscore_proteome <- proteome$tscore[match(transcriptome$gene.sex.tissue.time, proteome$gene.sex.tissue.time)]
transcriptome$logFC_se_proteome <- proteome$logFC_se[match(transcriptome$gene.sex.tissue.time, proteome$gene.sex.tissue.time)]
transcriptome <- transcriptome[!is.na(transcriptome$logFC_proteome),]
transcriptome$logFC_diff_prot.minus.trans <- transcriptome$logFC_proteome - transcriptome$logFC
transcriptome$ztscore_diff_prot.minus.trans <- transcriptome$tscore_proteome - transcriptome$zscore

#add extra identified for cases where duplicate gene_id for different feature_id?
unique(transcriptome$gene_id)

#munge into list format
tissues <- sort(unique(transcriptome$tissue))
transcriptome_list <- lapply(tissues, function(tissue_label) transcriptome[transcriptome$tissue == tissue_label,])
names(transcriptome_list) <- tissues


#### plot basic figure ####



#specify colors
nnmap <- as.data.frame(bic_animal_tissue_code[,c(4:5,2,6)])
cols = list(Tissue=tissue_cols[tissues], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
ome_cols <- c(trnscrpt = "#1787b7", prot_col = "#fbaa00", trnscrpt_zscore = "#d03a18", zscore_logFC_corr = "#7851A9", 
              trnscrpt_prot_corr = "#197419", tscore_proteome = "#4B0082", tz_diff = "#8F00FF")
col_opacity <- 0.5

proteins_determine_cluster_assignment <- T
bufferx <- buffery <- 0.1
buffer_hist <- 0.05

plot_trnscript_logFC <- T
plot_protein_logFC <- T
plot_logFC_diff <- T
plot_trnscript_zscore <- F
plot_protein_tscore <- F
plot_TZ_diff <- F
use_zscore_corrs <- F

plot_as_polygons <- T
one_word_per_line <- F
swap_logFC_shrunkLogFC = F
PiC = T
# lower_hist_bin_by <- 0.5
polygon_interval <- 0.20

foreach(cli=1:15, .packages = c("MotrpacBicQC", "data.table")) %dopar% {
# for(cli in 15){
  
grDevices::cairo_pdf(filename = paste0("~/Documents/transcr_v_protein/", ifelse(plot_trnscript_zscore | plot_protein_tscore, "usingT-Z-Scores_", ""),
                                       ifelse(proteins_determine_cluster_assignment, "proteins-determine-cluster-assignment_", ""), 
                                       ifelse(swap_logFC_shrunkLogFC, "shrunkLogFC_", ""), "cluster_", cli,".pdf"), width = 850 / 72, height = 800 / 72, family="Arial Unicode MS")
par(mfrow = c(2,2))

for(tissue in tissues){

cat(paste0("(", cli, ", ", which(tissues == tissue), ") "))
  
if(tissue %in% tissues[1:2]){
  par(mar = c(1,4,5,4))
} else {
  par(mar = c(2,4,4,4))
}
# tissue <- tissues[1]
tissue_abbrev <- nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)]
tissue_nicename <- nnmap$bic_tissue_name[match(tissue, nnmap$tissue_name_release)]
if(length(grep(pattern = "Powder", x = tissue_nicename)) > 0){
  tissue_nicename <- paste0(strsplit(tissue_nicename, split = " ")[[1]][-length(strsplit(tissue_nicename, split = " ")[[1]])], collapse = " ")
}
tissue_col <- nnmap$tissue_hex_colour[match(tissue, nnmap$tissue_name_release)]

# table(transcriptome_list[[tissue]]$comparison_group)

focal_genes <- cluster_membership[cluster_membership$cluster == cli & cluster_membership$tissue == tissue_abbrev & cluster_membership$ome == 
                                    ifelse(proteins_determine_cluster_assignment, "PROT", "TRNSCRPT"), "feature_ID"]
focal_genes <- f2gm$gene_symbol[match(focal_genes, f2gm$feature_ID)]
focal_genes <- focal_genes[!is.na(focal_genes)]
focal_genes <- focal_genes[which(focal_genes %in% transcriptome_list[[tissue]]$gene_id)]
fg_inds <- which(transcriptome_list[[tissue]]$gene_id %in% focal_genes)
transcr_subset <- transcriptome_list[[tissue]][fg_inds,]
transcr_subset <- transcr_subset[order(transcr_subset$gene_id),]

#hack shrunken logFCs in
if(swap_logFC_shrunkLogFC){
  transcr_subset$logFC <- transcr_subset$shrunk_logFC
}


#get rid of invariant genes
logFC_prots <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), logFC_proteome])
logFC_trans <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), logFC])
genes_for_badness_comparison <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), gene_id])
all.eq.numvec <- function(x) all( abs(x - mean(x)) < .Machine$double.eps ^ 0.5 )
bad_inds <- which(apply(logFC_prots, 1, all.eq.numvec) | apply(logFC_trans, 1, all.eq.numvec))
if(length(bad_inds) > 0){
  bad_genes <- genes_for_badness_comparison[bad_inds,1]
  transcr_subset <- transcr_subset[-which(transcr_subset$gene_id %in% bad_genes),]
}


#plotting info
if(plot_as_polygons){
  quantiles_to_use <- c((1 - polygon_interval) / 2, (1/2 + polygon_interval/2))
  
  range_transcr_logFC <- lapply(paste0(c(1,2,4,8), "w"), function(tmpt) sapply(c("female", "male"), function(sx) 
    quantile(transcr_subset$logFC[(transcr_subset$sex == sx & transcr_subset$comparison_group == tmpt)], quantiles_to_use)))
  range_transcr_logFC <- do.call(what = cbind, range_transcr_logFC)
  range_transcr_logFC <- c(min(range_transcr_logFC[1,]), max(range_transcr_logFC[2,]))
  
  range_prot_logFC <- lapply(paste0(c(1,2,4,8), "w"), function(tmpt) sapply(c("female", "male"), function(sx) 
    quantile(transcr_subset$logFC_proteome[(transcr_subset$sex == sx & transcr_subset$comparison_group == tmpt)], quantiles_to_use)))
  range_prot_logFC <- do.call(what = cbind, range_prot_logFC)
  range_prot_logFC <- c(min(range_prot_logFC[1,]), max(range_prot_logFC[2,]))
  logFC_median8w <- sapply(c("female", "male"), function(sx) quantile(transcr_subset$logFC[(transcr_subset$sex == sx & transcr_subset$comparison_group == "8w")], 0.5))
  
  range_prot_logFC_diff <- lapply(paste0(c(1,2,4,8), "w"), function(tmpt) sapply(c("female", "male"), function(sx) 
    quantile(transcr_subset$logFC_diff_prot.minus.trans[(transcr_subset$sex == sx & transcr_subset$comparison_group == tmpt)], quantiles_to_use)))
  range_prot_logFC_diff <- do.call(what = cbind, range_prot_logFC_diff)
  range_prot_logFC_diff <- c(min(range_prot_logFC_diff[1,]), max(range_prot_logFC_diff[2,]))
  
  range_transcr_zscore <- lapply(paste0(c(1,2,4,8), "w"), function(tmpt) sapply(c("female", "male"), function(sx) 
    quantile(transcr_subset$zscore[(transcr_subset$sex == sx & transcr_subset$comparison_group == tmpt)], quantiles_to_use)))
  range_transcr_zscore <- do.call(what = cbind, range_transcr_zscore)
  range_transcr_zscore <- c(min(range_transcr_zscore[1,]), max(range_transcr_zscore[2,]))
  zscore_median8w <- sapply(c("female", "male"), function(sx) quantile(transcr_subset$zscore[(transcr_subset$sex == sx & transcr_subset$comparison_group == "8w")], 0.5))
  
  range_prot_tscore <- lapply(paste0(c(1,2,4,8), "w"), function(tmpt) sapply(c("female", "male"), function(sx) 
    quantile(transcr_subset$tscore_proteome[(transcr_subset$sex == sx & transcr_subset$comparison_group == tmpt)], quantiles_to_use)))
  range_prot_tscore <- do.call(what = cbind, range_prot_tscore)
  range_prot_tscore <- c(min(range_prot_tscore[1,]), max(range_prot_tscore[2,]))
  tscore_median8w <- sapply(c("female", "male"), function(sx) quantile(transcr_subset$tscore_proteome[(transcr_subset$sex == sx & transcr_subset$comparison_group == "8w")], 0.5))
  
  
} else {
  quantiles_to_use <- c(0, 1)
  range_transcr_logFC <- quantile(transcr_subset$logFC, quantiles_to_use)
  range_prot_logFC <- quantile(transcr_subset$logFC_proteome, quantiles_to_use)
  range_transcr_zscore <- quantile(transcr_subset$zscore, quantiles_to_use)
}
range_logFC <- c(min(range_transcr_logFC, range_prot_logFC, range_prot_logFC_diff), max(range_transcr_logFC, range_prot_logFC, range_prot_logFC_diff))

if(plot_trnscript_zscore & plot_protein_tscore){
  range_transcr_zscore <- c(min(range_transcr_zscore, range_prot_tscore), max(range_transcr_zscore, range_prot_tscore))
  range_logFC <- range_transcr_zscore
} else{
  #maybe make it so middle of logFC and zscore coincide? 
  range_transcr_zscore <- max(c(range_transcr_zscore[1] / range_logFC[1], range_transcr_zscore[2] / range_logFC[2], zscore_median8w[1] / logFC_median8w[1])) * range_logFC
}

#start plotting -- make frame
plot(1,1,xlim = c(0,2), ylim = c(0,3), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#timepoint labels
ticks_timepoints <- c(1-1+c(0,1,2,4,8)/8*(1-bufferx)+bufferx/2, 2-1+c(0,1,2,4,8)/8*(1-bufferx)+bufferx/2)
segments(x0 = ticks_timepoints, x1 = ticks_timepoints, lty = 1, lwd = 1.5, y0 = 1, y1 = 0.965)
segments(x0 = ticks_timepoints, x1 = ticks_timepoints, y0 = 1, y1 = 3, lwd = 1, col = "grey50", lty = 3)
text(labels = paste0(c(0,1,2,4,8,0,1,2,4,8), "w"), x = ticks_timepoints, y = 0.99, pos = 1, xpd = NA, cex = 0.75)

#logFC ticks
if(plot_trnscript_logFC){
  ticks_logFC <- seq(range_logFC[1], range_logFC[2], length.out = 8)
  diff_ticks_logFC <- diff(ticks_logFC)[1]
  diff_ticks_logFC <- round(diff_ticks_logFC, digits = ceiling(-log10(diff_ticks_logFC)) + 1) 
  ticks_logFC <- sort(c(seq(0, range_logFC[2], by = diff_ticks_logFC*sign(range_logFC[2])), seq(0, range_logFC[1], by = diff_ticks_logFC*sign(range_logFC[1]))[-1]))
  segments(x0 = 0, x1 = -0.035, lty = 1, lwd = 1.5,
           y0 = (ticks_logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1, 
           y1 = (ticks_logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1)
  text(labels = ticks_logFC, x = -0.01, y = (ticks_logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1, pos = 2, xpd = NA, cex = 0.75)
  text(x = -0.2, y = 2.05, pos = 2, labels = latex2exp::TeX(paste0(ifelse(swap_logFC_shrunkLogFC, "shrunk ", ""), "log_2(FC)")), srt = 90, xpd = NA)
}

#zscore ticks
if(plot_trnscript_zscore){
  
  ticks_zscore <- seq(range_transcr_zscore[1], range_transcr_zscore[2], length.out = 8)
  diff_ticks_zscore <- diff(ticks_zscore)[1]
  diff_ticks_zscore <- round(diff_ticks_zscore, digits = ceiling(-log10(diff_ticks_zscore)) + 1) 
  ticks_zscore <- sort(c(seq(0, range_transcr_zscore[2], by = diff_ticks_zscore*sign(range_transcr_zscore[2])), seq(0, range_transcr_zscore[1], 
                                                                                                                    by = diff_ticks_zscore*sign(range_transcr_zscore[1]))[-1]))
  
  
  segments(x0 = 2, x1 = 2.035, lty = 1, lwd = 1.5,
           y0 = (ticks_zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1, 
           y1 = (ticks_zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1)
  text(labels = ticks_zscore, x = 2.01, y = (ticks_zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1, pos = 4, xpd = NA, cex = 0.75)
  text(x = 2.175, y = 3.0, pos = 4, labels = "T/Z-Score (Signed T/Probit-transformed P-Value)", srt = 270, xpd = NA)
}

#title
shadowtext(1, 3, pos = 3, labels = tissue_nicename, col = tissue_col, xpd = NA, cex = 3)

for(sex in 1:2){
  sex_inds <- which(transcr_subset$sex == c("female", "male")[sex])
    for(timepoint in 1:3){
      timepoint_inds_1 <- which(transcr_subset$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
      timepoint_inds_2 <- which(transcr_subset$comparison_group == paste0(c(1,2,4,8), "w")[timepoint+1])
      data_subset_1 <- transcr_subset[intersect(sex_inds, timepoint_inds_1),]
      data_subset_1_bounds <- sapply(c("logFC", "logFC_proteome", "zscore", "logFC_diff_prot.minus.trans", "tscore_proteome", "ztscore_diff_prot.minus.trans"), function(x)  
        quantile(x = as.data.frame(data_subset_1)[,x], quantiles_to_use))
      data_subset_2 <- transcr_subset[intersect(sex_inds, timepoint_inds_2),]
      data_subset_2_bounds <- sapply(c("logFC", "logFC_proteome", "zscore", "logFC_diff_prot.minus.trans", "tscore_proteome", "ztscore_diff_prot.minus.trans"), function(x)  
        quantile(x = as.data.frame(data_subset_2)[,x], quantiles_to_use))
      
      if(timepoint == 1){
        #plot 0->1 line
        if(plot_trnscript_zscore){
          if(plot_as_polygons){
            polygon(x = c(0*(1-bufferx)+bufferx/2+sex-1, 
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          0*(1-bufferx)+bufferx/2+sex-1),
                    y = c((0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[1,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[2,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                          (0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1),
                    col = adjustcolor(ome_cols[3], col_opacity), border = NA)
          } else {
            segments(x0 = 0*(1-bufferx)+bufferx/2+sex-1, 
                     y0 = (0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1, 
                     x1 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                     y1 = (data_subset_1$zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                     col = adjustcolor(ome_cols[3], col_opacity))    
          }
          
        }
        if(plot_protein_tscore){
          if(plot_as_polygons){
            polygon(x = c(0*(1-bufferx)+bufferx/2+sex-1, 
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          0*(1-bufferx)+bufferx/2+sex-1),
                    y = c((0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[1,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[2,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                          (0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1),
                    col = adjustcolor(ome_cols["tscore_proteome"], col_opacity), border = NA)
          } else {
            segments(x0 = 0*(1-bufferx)+bufferx/2+sex-1, 
                     y0 = (0 - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2*(1-buffery)+buffery/2 + 1, 
                     x1 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                     y1 = (data_subset_1$tscore_proteome - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                     col = adjustcolor(ome_cols["tscore_proteome"], col_opacity))    
          }
          
        }
        if(plot_trnscript_logFC){
          if(plot_as_polygons){
            polygon(x = c(0*(1-bufferx)+bufferx/2+sex-1, 
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          0*(1-bufferx)+bufferx/2+sex-1),
                    y = c((0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[1,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[2,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1),
                    col = adjustcolor(ome_cols[1], col_opacity), border = NA)
          } else {
            segments(x0 = 0*(1-bufferx)+bufferx/2+sex-1, 
                     y0 = (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1, 
                     x1 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                     y1 = (data_subset_1$logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                     col = adjustcolor(ome_cols[1], col_opacity))
          }
          
        }
        if(plot_protein_logFC){
          if(plot_as_polygons){
            polygon(x = c(0*(1-bufferx)+bufferx/2+sex-1, 
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          0*(1-bufferx)+bufferx/2+sex-1),
                    y = c((0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[1,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[2,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1),
                    col = adjustcolor(ome_cols[2], col_opacity), border = NA)
          } else {
            segments(x0 = 0*(1-bufferx)+bufferx/2+sex-1, 
                     y0 = (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1, 
                     x1 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                     y1 = (data_subset_1$logFC_proteome - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                     col = adjustcolor(ome_cols[2], col_opacity))
          }
         
        }
        if(plot_logFC_diff){
          if(plot_as_polygons){
            polygon(x = c(0*(1-bufferx)+bufferx/2+sex-1, 
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2,
                          0*(1-bufferx)+bufferx/2+sex-1),
                    y = c((0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[1,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (data_subset_1_bounds[2,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                          (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1),
                    col = adjustcolor(ome_cols[5], col_opacity), border = NA)
          } else {
            segments(x0 = 0*(1-bufferx)+bufferx/2+sex-1, 
                     y0 = (0 - range_logFC[1]) / diff(range_logFC) * 2*(1-buffery)+buffery/2 + 1, 
                     x1 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                     y1 = (data_subset_1$logFC_diff_prot.minus.trans - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                     col = adjustcolor(ome_cols[5], col_opacity))
          }
          
        }
        
      }
      
      
      
      if(plot_trnscript_zscore){
        
        if(plot_as_polygons){
          polygon(x = c(sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2),
                  y = c((data_subset_1_bounds[1,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[1,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[2,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_1_bounds[2,"zscore"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1),
                  col = adjustcolor(ome_cols[3], col_opacity), border = NA)
        } else {
          segments(x0 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                   y0 = (data_subset_1$zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1, 
                   x1 = sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2, 
                   y1 = (data_subset_2$zscore - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                   col = adjustcolor(ome_cols[3], col_opacity))
        }
        
      }
      
      if(plot_protein_tscore){
        
        if(plot_as_polygons){
          polygon(x = c(sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2),
                  y = c((data_subset_1_bounds[1,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[1,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[2,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_1_bounds[2,"tscore_proteome"] - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1),
                  col = adjustcolor(ome_cols["tscore_proteome"], col_opacity), border = NA)
        } else {
          segments(x0 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                   y0 = (data_subset_1$tscore_proteome - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1, 
                   x1 = sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2, 
                   y1 = (data_subset_2$tscore_proteome - range_transcr_zscore[1]) / diff(range_transcr_zscore) * 2 *(1-buffery)+buffery/2 + 1,
                   col = adjustcolor(ome_cols["tscore_proteome"], col_opacity))
        }
        
      }
      
      if(plot_trnscript_logFC){
        
        if(plot_as_polygons){
          polygon(x = c(sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2),
                  y = c((data_subset_1_bounds[1,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[1,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[2,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_1_bounds[2,"logFC"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1),
                  col = adjustcolor(ome_cols[1], col_opacity), border = NA)
        } else {
          segments(x0 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                   y0 = (data_subset_1$logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1, 
                   x1 = sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2, 
                   y1 = (data_subset_2$logFC - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                   col = adjustcolor(ome_cols[1], col_opacity))
        }
        
      }
      if(plot_protein_logFC){
        
        if(plot_as_polygons){
          polygon(x = c(sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2),
                  y = c((data_subset_1_bounds[1,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[1,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[2,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_1_bounds[2,"logFC_proteome"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1),
                  col = adjustcolor(ome_cols[2], col_opacity), border = NA)
        } else {
          segments(x0 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                   y0 = (data_subset_1$logFC_proteome - range_logFC[1]) / diff(range_logFC) * 2 * (1-buffery)+buffery/2 + 1, 
                   x1 = sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2, 
                   y1 = (data_subset_2$logFC_proteome - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                   col = adjustcolor(ome_cols[2], col_opacity))  
        }
        
      }
      
      if(plot_logFC_diff){
        
        if(plot_as_polygons){
          polygon(x = c(sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2,
                        sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2),
                  y = c((data_subset_1_bounds[1,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[1,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_2_bounds[2,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                        (data_subset_1_bounds[2,"logFC_diff_prot.minus.trans"] - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1),
                  col = adjustcolor(ome_cols[5], col_opacity), border = NA)
        } else {
          segments(x0 = sex-1+c(1,2,4,8)[timepoint]/8*(1-bufferx)+bufferx/2, 
                   y0 = (data_subset_1$logFC_diff_prot.minus.trans - range_logFC[1]) / diff(range_logFC) * 2 * (1-buffery)+buffery/2 + 1, 
                   x1 = sex-1+c(1,2,4,8)[timepoint+1]/8*(1-bufferx)+bufferx/2, 
                   y1 = (data_subset_2$logFC_diff_prot.minus.trans - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1,
                   col = adjustcolor(ome_cols[5], col_opacity))  
        }
      }
    }
}

#disinclude timepoint 0 for corrs
transcr_subset <- transcr_subset[order(transcr_subset$sex),]
logFC_prots <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), logFC_proteome])
logFC_trans <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), logFC])
zscore_trans <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), zscore])
tscore_prots <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), tscore_proteome])
gene_ids <- sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), 
                                                     paste0(gene_id, "-", list(male = "\u2642", female = "\u2640")[sex])])[,1]
# sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), sex])
# sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), gene_id])
# sapply(1:4,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), gene.sex.tissue.time])



if(PiC){
  logFC_prots <- cbind(0, logFC_prots)
  logFC_trans <- cbind(0, logFC_trans)
  zscore_trans <- cbind(0, zscore_trans)
  tscore_prots <- cbind(0, tscore_prots)
  
  # logFC_prots <- t(sapply(1:nrow(logFC_prots), function(ri) sapply(2:ncol(logFC_prots), function(ci) logFC_prots[ri,ci] - logFC_prots[ri,ci-1])))
  # logFC_trans <- t(sapply(1:nrow(logFC_trans), function(ri) sapply(2:ncol(logFC_trans), function(ci) logFC_trans[ri,ci] - logFC_trans[ri,ci-1])))
  # zscore_trans <- t(sapply(1:nrow(zscore_trans), function(ri) sapply(2:ncol(zscore_trans), function(ci) zscore_trans[ri,ci] - zscore_trans[ri,ci-1])))
  
  logFC_prots <- t(sapply(1:nrow(logFC_prots), function(ri) sapply(2:ncol(logFC_prots), function(ci) (logFC_prots[ri,ci] - logFC_prots[ri,ci-1]) / sqrt(c(1,1,2,4)[ci-1]) )))
  logFC_trans <- t(sapply(1:nrow(logFC_trans), function(ri) sapply(2:ncol(logFC_trans), function(ci) (logFC_trans[ri,ci] - logFC_trans[ri,ci-1]) / sqrt(c(1,1,2,4)[ci-1]) )))
  zscore_trans <- t(sapply(1:nrow(zscore_trans), function(ri) sapply(2:ncol(zscore_trans), function(ci) (zscore_trans[ri,ci] - zscore_trans[ri,ci-1]) / sqrt(c(1,1,2,4)[ci-1]) )))
  tscore_prots <- t(sapply(1:nrow(tscore_prots), function(ri) sapply(2:ncol(tscore_prots), function(ci) (tscore_prots[ri,ci] - tscore_prots[ri,ci-1]) / sqrt(c(1,1,2,4)[ci-1]) )))
}

if(use_zscore_corrs){
  corrs_prots.trans <- sapply(1:nrow(zscore_trans), function(fg) cor(tscore_prots[fg,], zscore_trans[fg,]))
  corrs_trans_logFC.zscores <- sapply(1:nrow(zscore_trans), function(fg) cor(zscore_trans[fg,], logFC_trans[fg,]))
} else {
  corrs_prots.trans <- sapply(1:nrow(logFC_prots), function(fg) cor(logFC_prots[fg,], logFC_trans[fg,]))
  corrs_trans_logFC.zscores <- sapply(1:nrow(zscore_trans), function(fg) cor(zscore_trans[fg,], logFC_trans[fg,]))
}

hist_breaks <- -5:5/5
hist_corrs_trans_logFC.zscores <- hist(corrs_trans_logFC.zscores, breaks = hist_breaks, plot = F)
hist_corrs_prots.trans <- hist(corrs_prots.trans, breaks = hist_breaks, plot = F)

#get gene names sorted
gene_bin_inds <- (sapply(corrs_prots.trans, function(cr) sum(cr > hist_breaks)))
genes_in_bins <- sapply(1:(length(hist_breaks)-1), function(gbi) sort(unique(gene_ids[gene_bin_inds == gbi])))
genes_in_bins_nosex <- sapply(genes_in_bins, function(gb) (sapply(gb, function(g) substr(g, 1, nchar(g) - 2))))
for(gbi in 1:length(genes_in_bins)){
  dupes <- as.logical(duplicated(genes_in_bins_nosex[[gbi]]) + rev(duplicated(rev(genes_in_bins_nosex[[gbi]]))))
  # genes_in_bins[[gbi]][dupes] <- paste0(genes_in_bins_nosex[[gbi]][dupes], "-\u26A5")
  genes_in_bins[[gbi]][dupes] <- paste0(genes_in_bins_nosex[[gbi]][dupes], "-\u2642\u2640")
  genes_in_bins[[gbi]] <- sort(unique(genes_in_bins[[gbi]]))
}


#make a little sub-window for sex-proportions
buffer_dens <- 0.1
kernel_width_adjustment <- 1
beta_regularizer_params <- c(1.5, 1.5)
xl = 0; xr = 0.875; yb = 0.5; yt = 0.85
xl0 = xl; xr0 = xr; yb0 = yb; yt0 = yt
xl = xl + buffer_dens/2 * (xr0-xl0)
xr = xr - buffer_dens/2 * (xr0-xl0)
rect(xleft = xl0, xright = xr0, ybottom = yb0, ytop = yt0, lwd = 0.5, xpd = NA, col = "grey98")
sexes <- sapply(1:1,function(tmpt) transcr_subset[transcr_subset$comparison_group == paste0(c(1,2,4,8)[tmpt], "w"), sex])
m_pt_corrs <- corrs_prots.trans[sexes == "male"]
f_pt_corrs <- corrs_prots.trans[sexes == "female"]
m_pt_corrs_logit01 <- logit((m_pt_corrs + 1) / 2)
f_pt_corrs_logit01 <- logit((f_pt_corrs + 1) / 2)
m_pt_corrs_logit01_density <- (density(m_pt_corrs_logit01, adjust = kernel_width_adjustment))
f_pt_corrs_logit01_density <- (density(f_pt_corrs_logit01, adjust = kernel_width_adjustment))
m_pt_corrs_density <- list(x = invlogit(m_pt_corrs_logit01_density$x), y = m_pt_corrs_logit01_density$y / invlogit_deriv(m_pt_corrs_logit01_density$x))
f_pt_corrs_density <- list(x = invlogit(f_pt_corrs_logit01_density$x), y = f_pt_corrs_logit01_density$y / invlogit_deriv(f_pt_corrs_logit01_density$x))

#regularize with beta
f_pt_corrs_density$y <- f_pt_corrs_density$y * dbeta(x = f_pt_corrs_density$x, shape1 = beta_regularizer_params[1], shape2 = beta_regularizer_params[2])
m_pt_corrs_density$y <- m_pt_corrs_density$y * dbeta(x = m_pt_corrs_density$x, shape1 = beta_regularizer_params[1], shape2 = beta_regularizer_params[2])

#rescale to (-1,1)
m_pt_corrs_density$x <- m_pt_corrs_density$x * 2 - 1
f_pt_corrs_density$x <- f_pt_corrs_density$x * 2 - 1

#law of total prob
# m_pt_corrs_density$y  <- m_pt_corrs_density$y / sum(m_pt_corrs_density$y * 2 / length(m_pt_corrs_density$x)) #whoops these are not equally spaced on the raw scale
m_pt_corrs_density$y  <- m_pt_corrs_density$y / sum(m_pt_corrs_density$y * diff(c(0,m_pt_corrs_density$x)))
# f_pt_corrs_density$y  <- f_pt_corrs_density$y / sum(f_pt_corrs_density$y * 2 / length(f_pt_corrs_density$x)) #whoops these are not equally spaced on the raw scale
f_pt_corrs_density$y  <- f_pt_corrs_density$y / sum(f_pt_corrs_density$y * diff(c(0,f_pt_corrs_density$x)))

#add bookends for polygon
m_pt_corrs_density$x <- c(-1,m_pt_corrs_density$x,1)
f_pt_corrs_density$x <- c(-1,f_pt_corrs_density$x,1)
m_pt_corrs_density$y <- c(0,m_pt_corrs_density$y,0)
f_pt_corrs_density$y <- c(0,f_pt_corrs_density$y,0)
# plot(m_pt_corrs_density$x, m_pt_corrs_density$y, type = "l")
# lines(f_pt_corrs_density$x, f_pt_corrs_density$y, type = "l")
max_y <- max(c(f_pt_corrs_density$y, m_pt_corrs_density$y)) * 1.2
polygon(x = (m_pt_corrs_density$x + 1) / 2 * (xr - xl) + xl, 
        y = (m_pt_corrs_density$y) / max_y * (yt - yb) + yb, 
        col = adjustcolor(cols$Sex["male"], col_opacity), border = NA)
polygon(x = (f_pt_corrs_density$x + 1) / 2 * (xr - xl) + xl, 
        y = (f_pt_corrs_density$y) / max_y * (yt - yb) + yb, 
        col = adjustcolor(cols$Sex["female"], col_opacity), border = NA)
dens_tick_xlocs <- (c(-1,0,1) + 1) / 2 * (xr - xl) + xl
segments(x0 = dens_tick_xlocs, x1 = dens_tick_xlocs, y0 = yb, y1 = yb-0.015, lwd = 0.5)
text(labels = c(-1,0,1), x = dens_tick_xlocs, y = yb + 0.025, pos = 1, cex = 0.4)
dens_y_ndec <- ceiling(log(max_y))
dens_y_ndec <- 1
dens_tick_y <- seq(0, max_y*0.85, by = diff(round(seq(0, max_y, length.out = 4), dens_y_ndec))[1])
dens_tick_ylocs <- (dens_tick_y) / max_y * (yt - yb) + yb
segments(x0 = xr0, x1 = xr0 + 0.0125, y0 = dens_tick_ylocs, y1 = dens_tick_ylocs, lwd = 0.5)
text(labels = dens_tick_y, x = xr0 - 0.0225, y = dens_tick_ylocs, pos = 4, cex = 0.4)
text_cols(string = "\u2642\u2640 Corr. Density", srt = 270, x = xr0 + 0.065, y = (yt0 * 0.325 + yb0 * 0.675), 
          cols = c(cols$Sex[c("male", "female")], rep(1, nchar("\u2642\u2640 Corr. Density") - 2)),
          cex = 0.4, pos = 3)


#figure out height of histogram
histogram_range_rescaler <- 1.4 * max(1, max(hist_corrs_prots.trans$counts[1:5]) /  max(hist_corrs_prots.trans$counts) * 2)
range_hist_corrs_prots.trans <- c(0,max(hist_corrs_prots.trans$counts)) * histogram_range_rescaler
range_hist_corrs_trans_logFC.zscores <- c(0,max(hist_corrs_trans_logFC.zscores$counts)) * histogram_range_rescaler


#plot histogram blocks
for(bin in 1:(length(hist_corrs_prots.trans$breaks) - 1)){
  rect(xleft = (hist_corrs_prots.trans$breaks[bin] + 1)*(1-buffer_hist)+buffer_hist,
       xright = (hist_corrs_prots.trans$breaks[bin+1] + 1)*(1-buffer_hist)+buffer_hist,
       # ybottom = (-lower_hist_bin_by/range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2,
       ybottom = 0,
       ytop = (hist_corrs_prots.trans$counts[bin] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2,
       col = adjustcolor(ome_cols[5], 0.75),
       border = NA
        )
}
# for(bin in 1:(length(hist_corrs_trans_logFC.zscores$breaks) - 1)){
#   rect(xleft = (hist_corrs_trans_logFC.zscores$breaks[bin] + 1)*(1-buffer_hist)+buffer_hist,
#        xright = (hist_corrs_trans_logFC.zscores$breaks[bin+1] + 1)*(1-buffer_hist)+buffer_hist,
#        ybottom = 0*(1-buffer_hist)+buffer_hist/2,
#        ytop = (hist_corrs_trans_logFC.zscores$counts[bin] / range_hist_corrs_trans_logFC.zscores[2])*(1-buffer_hist)+buffer_hist/2,
#        col = adjustcolor(ome_cols[4], 0.75),
#        border = NA
#   )
# }



#put words in histogram blocks
for(gbi in 1:length(genes_in_bins)){
  
  #write number of genes of each sex above rectangle
  n_male <- length(grep(pattern = "\u2642", genes_in_bins[[gbi]]))
  n_female <- length(grep(pattern = "\u2640", genes_in_bins[[gbi]]))
  sex_summary_text <- paste0(n_male, "\u2642,", n_female, "\u2640")
  sex_summary_text_cols <- c(rep(1,nchar(n_male)), cols$Sex["male"], rep(1,nchar(n_female)+1), cols$Sex["female"])
  text_cols(string = sex_summary_text, pos = 3, cex = 0.75, cols = sex_summary_text_cols,
            x = (hist_corrs_prots.trans$mids[gbi] + 1)*(1-buffer_hist)+buffer_hist, 
            y = (hist_corrs_prots.trans$counts[gbi] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - 0.05)
  
  if(length(genes_in_bins[[gbi]]) == 0){
    next()
  }
  
  if(one_word_per_line){
    for(gi in 1:length(genes_in_bins[[gbi]])){
      bin_width <- diff((hist_corrs_prots.trans$mids[1:2] + 1)*(1-buffer_hist)+buffer_hist) * 0.85
      max_bin_height <- (hist_corrs_prots.trans$counts[gbi] / range_hist_corrs_prots.trans[2])*(1-buffer_hist) / length(genes_in_bins[[gbi]]) * 0.5
      text_cex <- min(bin_width / strwidth(genes_in_bins[[gbi]][gi], units = "user"),
                      max_bin_height / strheight(genes_in_bins[[gbi]][gi], units = "user"))
      if(length(grep(pattern = "-\u2642\u2640", genes_in_bins[[gbi]][gi])) > 0){
        text(labels = genes_in_bins[[gbi]][gi], col = cols$Sex["female"], pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist - bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
        text(labels = substr(genes_in_bins[[gbi]][gi], start = 1, nchar((genes_in_bins[[gbi]][gi]))-1), col = cols$Sex["male"], 
             pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist- bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
        text(labels = substr(genes_in_bins[[gbi]][gi], start = 1, nchar((genes_in_bins[[gbi]][gi]))-2), col = "white", 
             pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist- bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
      } else if(length(grep(pattern = "\u2640", genes_in_bins[[gbi]][gi])) > 0){
        text(labels = genes_in_bins[[gbi]][gi], col = cols$Sex["female"], pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist - bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
        text(labels = substr(genes_in_bins[[gbi]][gi], start = 1, nchar((genes_in_bins[[gbi]][gi]))-1), col = "white", 
             pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist- bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
      } else{
        text(labels = genes_in_bins[[gbi]][gi], col = cols$Sex["male"], pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist - bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
        text(labels = substr(genes_in_bins[[gbi]][gi], start = 1, nchar((genes_in_bins[[gbi]][gi]))-1), col = "white", 
             pos = 4, cex = text_cex, family = "Arial Unicode MS",
             x = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist- bin_width/6.5,
             y = ((hist_corrs_prots.trans$counts[gbi] - gi + 1.35) / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - max_bin_height*1.75)      
      }
      
    }
  } else {
    txt <- genes_in_bins[[gbi]]
    
    #truncate long gene names
    max_gene_name_length <- 10
    short_txt <- sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][-length(strsplit(stxt, "-")[[1]])])
    short_txt <- sapply(short_txt, function(stxt) paste0(substr(stxt, 1, max_gene_name_length), ifelse(nchar(stxt) > max_gene_name_length, ".", "")))
    txt <- sapply(1:length(txt), function(txti) paste0(short_txt[txti], "-", strsplit(txt[txti], "-")[[1]][length(strsplit(txt[txti], "-")[[1]])]))
    
    rect_coords <- list(x0 = (hist_corrs_prots.trans$breaks[gbi] + 1)*(1-buffer_hist)+buffer_hist,
                        x1 = (hist_corrs_prots.trans$breaks[gbi+1] + 1)*(1-buffer_hist)+buffer_hist,
                        y0 = (hist_corrs_prots.trans$counts[gbi] / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2 - 0.01,
                        y1 = 0*(1-buffer_hist)+buffer_hist/2)
    optimal_word_placement_inf <- find_optimal_cex_and_lines(txt = txt, rect_coords = rect_coords, rect_rescaling_ratio = 0.925)
    #1 is male, 2 female, 3 both
    txt_sexes <- as.numeric(nchar(txt) - nchar(sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][1])) == 3) + 2
    txt_sexes[txt_sexes == 2] <- as.numeric(sapply(txt, function(stxt) strsplit(stxt, "-")[[1]][2])[txt_sexes == 2] == "\u2640" ) + 1
    sex_cols_map <- list(male = cols$Sex[c("male")], female = cols$Sex[c("female")], both = cols$Sex[c("male", "female")])
    cols_list <- lapply(1:length(txt), function(word_i) as.vector(c(rep(0, nchar(txt[word_i]) - 1 - floor(txt_sexes[word_i]/3)), unlist(sex_cols_map[txt_sexes[word_i]]))))
    text_wrapped_words(txt, rect_coords, optimal_word_placement_inf, justified = T, col = "white", str_width_lefter_start_ratio = 0.15, rect_rescaling_ratio = 0.925,
                       cols_list = cols_list, multicolor_words = T)
    
  }
}

#zero line
segments(x0 = 0, x1 = 2, col = "grey50", lty = 2, lwd = 1.5,
         y0 = (0 - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1, y1 = (0 - range_logFC[1]) / diff(range_logFC) * 2 *(1-buffery)+buffery/2 + 1)

#histogram axes
hist_vert_ticks <- round(seq(0, max(hist_corrs_prots.trans$counts) * histogram_range_rescaler, length.out = 7))
hist_vert_ticks <- seq(0, max(hist_corrs_prots.trans$counts) * histogram_range_rescaler, by = diff(hist_vert_ticks)[1])
hist_vert_ticks <- hist_vert_ticks[-length(hist_vert_ticks)]
hist_vert_ticks_locs <- (hist_vert_ticks / range_hist_corrs_prots.trans[2])*(1-buffer_hist)+buffer_hist/2
segments(x0 = 0, x1 = -0.035, lty = 1, lwd = 1.5,
         y0 = hist_vert_ticks_locs, y1 = hist_vert_ticks_locs)
text(labels = hist_vert_ticks, x = -0.015, pos = 2, 
     y = hist_vert_ticks_locs, xpd = NA, cex = 0.75)
text(x = -0.125, y = 0.55, pos = 2, labels = "Count", srt = 90, xpd = NA)

segments(x0 = (hist_breaks + 1)*(1-buffer_hist)+buffer_hist, x1 = (hist_breaks + 1)*(1-buffer_hist)+buffer_hist, lty = 1, lwd = 1.5,
         y0 = 0, y1 = -0.035)
text(labels = hist_breaks, x = (hist_breaks + 1)*(1-buffer_hist)+buffer_hist, 
     y = -0.01, pos = 1, xpd = NA, cex = 0.75)
if(use_zscore_corrs){
  text( latex2exp::TeX("Pearson Correlation between Transcript & Protein Z-/T-Score Trajectories"), x = 1, y = -0.15, pos = 1, xpd = NA)
} else {
  text( latex2exp::TeX("Pearson Correlation between Transcript & Protein log_2FC Trajectories"), x = 1, y = -0.15, pos = 1, xpd = NA)  
}



#plot legend title
legend_params <- legend(x = 0, y = 1, col = adjustcolor(ome_cols[c(plot_trnscript_logFC, plot_protein_logFC, plot_trnscript_zscore, F, plot_logFC_diff, plot_protein_tscore, plot_TZ_diff)], col_opacity), 
                        legend = c(latex2exp::TeX("transcript log_2(FC)"), 
                                   latex2exp::TeX("protein log_2(FC)"), 
                                   "transcript z-score",
                                   "transcript z-score / logFC corr",
                                   latex2exp::TeX("$\\Delta_{log_2(FC)}$ (prot-trans)"),
                                   "protein t-score",
                                   latex2exp::TeX("protein_{t-score}-transcr_{z-score}"))
                        [c(plot_trnscript_logFC, plot_protein_logFC, plot_trnscript_zscore, F, plot_logFC_diff, plot_protein_tscore, plot_TZ_diff)], 
                        pch = 15, pt.cex = 2, yjust = 0, cex = 0.75, plot = F)
shadowtext(x =legend_params$rect$left  - 0.025, y = legend_params$rect$top, bg = "white", r = 0.1,
           labels = paste0("(", quantiles_to_use[1], ", ", quantiles_to_use[2], ")\nQuantile Range"), 
           col = "grey20", pos = 4, cex = 0.7, font = 2, xpd = NA)

#plot legend




legend(x = 0, y = 1, col = adjustcolor(ome_cols[c(plot_trnscript_logFC, plot_protein_logFC, plot_trnscript_zscore, F, plot_logFC_diff, plot_protein_tscore, plot_TZ_diff)], col_opacity), 
       legend = c(latex2exp::TeX("transcript log_2(FC)"), 
                  latex2exp::TeX("protein log_2(FC)"), 
                  "transcript z-score",
                  "transcript z-score / logFC corr",
                  latex2exp::TeX("$\\Delta_{log_2(FC)}$ (prot-trans)"),
                  "protein t-score",
                  latex2exp::TeX("protein_{t-score}-transcr_{z-score}"))
                    [c(plot_trnscript_logFC, plot_protein_logFC, plot_trnscript_zscore, F, plot_logFC_diff, plot_protein_tscore, plot_TZ_diff)], 
       pch = 15, pt.cex = 2, yjust = 0, cex = 0.75)

# if(plot_trnscript_zscore){
#   legend(x = 0, y = 1, col = adjustcolor(ome_cols[c(1,2,5,3)], col_opacity), 
#          legend = c(latex2exp::TeX("transcript log_2(FC)"), latex2exp::TeX("protein log_2(FC)"), latex2exp::TeX("$\\Delta_{log_2(FC)}$ (prot-trans)"), "transcript z-score"), 
#          pch = 15, pt.cex = 2, yjust = 0, cex = 0.75)
# } else {
#   legend(x = 0, y = 1, col = adjustcolor(ome_cols[c(1,2,5)], col_opacity), 
#          legend = c(latex2exp::TeX("transcript log_2(FC)"), latex2exp::TeX("protein log_2(FC)"), latex2exp::TeX("$\\Delta_{log_2(FC)}$ (prot-trans)")), 
#          pch = 15, pt.cex = 2, yjust = 0, cex = 0.75)
# }

#plot sex symnbols
text(labels = c("\u2640", "\u2642"), x = c(-0.025,0.975), y = c(2.85,2.85), pos = 4, cex = 2.5, col = sex_cols[c('female','male')], family = "Arial Unicode MS")

#label time axis
text(x = 0, y = 0.9125, pos = 2, "time", font = 4, cex = 1.25, xpd = NA)
segments(x0 = -0.2, x1 = 2, y0 = 0.85, y1 = 0.85, lwd = 0.5, col = "grey0", xpd = NA)

#plot cluster label
if(tissue == "t55-gastrocnemius"){
  text(labels = paste0("Cluster ", cli), x = 2.45, y = 3.325, pos = 3, cex = 4, xpd = NA, font = 2, family = "Courier", col = "grey20")
}

#plot dividing line between tissues
if(tissue == "t55-gastrocnemius"){
  segments(x0 = -0.25, x1 = 5, y0 = -0.45, y1 = -0.45, col = "grey30", lty = 2, lwd = 2.5, xpd = NA)
  segments(x0 = 2.4, x1 = 2.4, y0 = 3.25, y1 = -4, col = "grey30", lty = 2, lwd = 2.5, xpd = NA)
}

#plot subframes
rect(xleft = 0, ybottom = 0, xright = 2, ytop = 3, lwd = 2)
rect(xleft = xl0, xright = xr0, ybottom = yb0, ytop = yt0, lwd = 0.5, xpd = NA)
segments(x0 = 0, y0 = 1, x1 = 2, y1 = 1, lwd = 2)
segments(x0 = 1, y0 = 1, x1 = 1, y1 = 3, lwd = 2)


}
dev.off()
}



pdftools::pdf_combine(paste0("~/Documents/transcr_v_protein/", ifelse(plot_trnscript_zscore | plot_protein_tscore, "usingT-Z-Scores_", ""),
                             ifelse(proteins_determine_cluster_assignment, "proteins-determine-cluster-assignment_", ""), 
                             ifelse(swap_logFC_shrunkLogFC, "shrunkLogFC_", ""),"cluster_", 1:15,".pdf"), 
                      output = paste0("~/Documents/transcr_v_protein/", ifelse(plot_trnscript_zscore | plot_protein_tscore, "usingT-Z-Scores_", ""),
                                      ifelse(proteins_determine_cluster_assignment, "proteins-determine-cluster-assignment_", ""), 
                                      ifelse(swap_logFC_shrunkLogFC, "shrunkLogFC_", ""),
                                      "transcriptome-vs-proteome_clusters1-15_(", polygon_interval ,"-QuantileRange).pdf"))

#### transcr / prot concordance ####

transcriptome_contrasts <- lapply(paste0(c(1,2,4,8), "w"), function(ti) transcriptome[transcriptome$comparison_group == ti,])
sapply(colnames(transcriptome_contrasts[[1]]), function(ci) all(transcriptome_contrasts[[2]][,..ci] == transcriptome_contrasts[[1]][,..ci]))
sapply(colnames(transcriptome_contrasts[[1]]), function(ci) all(transcriptome_contrasts[[3]][,..ci] == transcriptome_contrasts[[2]][,..ci]))
sapply(colnames(transcriptome_contrasts[[1]]), function(ci) all(transcriptome_contrasts[[4]][,..ci] == transcriptome_contrasts[[3]][,..ci]))
transcriptome_contrasts <- list(transcriptome_contrasts[[1]], transcriptome_contrasts[[1]], transcriptome_contrasts[[2]], transcriptome_contrasts[[3]], transcriptome_contrasts[[4]]) #get time 0 in there
transcriptome_contrasts[[1]][,"logFC"] <- 0
transcriptome_contrasts[[1]][,"logFC_proteome"] <- 0

transcriptome_contrasts_vals <- lapply(2:5, function(ti) (transcriptome_contrasts[[ti]]$logFC - transcriptome_contrasts[[ti-1]]$logFC) / sqrt(c(1,1,2,4)[ti-1]))
proteome_contrasts_vals <- lapply(2:5, function(ti) (transcriptome_contrasts[[ti]]$logFC_proteome - transcriptome_contrasts[[ti-1]]$logFC_proteome) / sqrt(c(1,1,2,4)[ti-1]))

# transcriptome_proteome_contrasts <- transcriptome_contrasts[[1]]
# transcriptome_proteome_contrasts <- lapply(1:length(proteome_contrasts_vals), function(x) 
#   transcriptome_proteome_contrasts[,c("feature_ID", "tissue", "sex", "gene_id", "tissue_abbreviation")])
# names(transcriptome_proteome_contrasts) <- paste0(c(1,2,4,8), "w-", c(0,1,2,4), "w")
# for(i in 1:length(transcriptome_contrasts_vals)){
#   transcriptome_proteome_contrasts[[i]]$transcriptome_contrasts <- transcriptome_contrasts_vals[[i]]  
#   transcriptome_proteome_contrasts[[i]]$proteome_contrasts <- proteome_contrasts_vals[[i]]  
#   transcriptome_proteome_contrasts[[i]]$contrast_i <- i
# }

transcriptome_proteome_contrasts_1dt <- transcriptome_contrasts[[1]][,c("feature_ID", "tissue", "sex", "gene_id", "tissue_abbreviation")]
transcriptome_proteome_contrasts_1dt$C1_prot <- proteome_contrasts_vals[[1]]
transcriptome_proteome_contrasts_1dt$C1_tran <- transcriptome_contrasts_vals[[1]]
transcriptome_proteome_contrasts_1dt$C2_prot <- proteome_contrasts_vals[[2]]
transcriptome_proteome_contrasts_1dt$C2_tran <- transcriptome_contrasts_vals[[2]]
transcriptome_proteome_contrasts_1dt$C3_prot <- proteome_contrasts_vals[[3]]
transcriptome_proteome_contrasts_1dt$C3_tran <- transcriptome_contrasts_vals[[3]]
transcriptome_proteome_contrasts_1dt$C4_prot <- proteome_contrasts_vals[[4]]
transcriptome_proteome_contrasts_1dt$C4_tran <- transcriptome_contrasts_vals[[4]]

transcriptome_proteome_contrasts_1dt$corr <- sapply(1:nrow(transcriptome_proteome_contrasts_1dt), function(ri) 
  cor(unlist(transcriptome_proteome_contrasts_1dt[ri, paste0("C", 1:4, "_tran")]),
      unlist(transcriptome_proteome_contrasts_1dt[ri, paste0("C", 1:4, "_prot")])))

bad_inds <- which(is.na(transcriptome_proteome_contrasts_1dt$corr))
transcriptome_proteome_contrasts_1dt <- transcriptome_proteome_contrasts_1dt[-bad_inds,]

# plot(unlist(transcriptome_contrasts_vals), unlist(proteome_contrasts_vals))
# plot(unlist(transcriptome_contrasts_vals)[transcriptome_proteome_contrasts_1dt$sex == "female"], unlist(proteome_contrasts_vals)[transcriptome_proteome_contrasts_1dt$sex == "female"])

hist(transcriptome_proteome_contrasts_1dt$corr)
hist(transcriptome_proteome_contrasts_1dt$corr[transcriptome_proteome_contrasts_1dt$tissue == tissues[4]])


cluster_membership$gene_symbols <- f2gm$gene_symbol[match(cluster_membership$feature_ID, f2gm$feature_ID)]
transcriptome_proteome_contrasts_1dt$trans_cluster <- sapply(1:nrow(transcriptome_proteome_contrasts_1dt), function(ri)
       cluster_membership$cluster[intersect(intersect(which(transcriptome_proteome_contrasts_1dt$tissue_abbreviation[ri] == cluster_membership$tissue), 
       which(cluster_membership$ome == "TRNSCRPT")),
       which(transcriptome_proteome_contrasts_1dt$gene_id[ri] == cluster_membership$gene_symbols))[1]
       ])
# zero_clusters <- which(sapply(1:length(transcriptome_proteome_contrasts_1dt$trans_cluster), function(x) length(transcriptome_proteome_contrasts_1dt$trans_cluster[x][[1]]) < 1))
# transcriptome_proteome_contrasts_1dt$trans_cluster[zero_clusters] <- 0
# nonzero_clusters <- which(sapply(1:length(transcriptome_proteome_contrasts_1dt$trans_cluster), function(x) length(transcriptome_proteome_contrasts_1dt$trans_cluster[x][[1]]) > 0))
# length(nonzero_clusters)
# sum(cluster_membership$assay == "TRNSCRPT" & cluster_membership$tissue %in% tissue_abbr[tissues])
# diff is because same gene symbol maps to multiple feature IDs

#compare to x-correlation function
cross_correlations <- sapply(1:nrow(transcriptome_contrasts[[1]]), function(ri) {temp <- ccf(
  c(transcriptome_contrasts[[1]]$logFC[ri], transcriptome_contrasts[[2]]$logFC[ri], 
    transcriptome_contrasts[[3]]$logFC[ri], transcriptome_contrasts[[4]]$logFC[ri]),
  c(transcriptome_contrasts[[1]]$logFC_proteome[ri], transcriptome_contrasts[[2]]$logFC_proteome[ri], 
    transcriptome_contrasts[[3]]$logFC_proteome[ri], transcriptome_contrasts[[4]]$logFC_proteome[ri]),
    plot = F);
  temp$lag[which.max(c(temp$acf))]}
)

transcriptome_proteome_contrasts_1dt$xcorr_tlag_max <- (cross_correlations)[-bad_inds]
transcriptome_proteome_contrasts_1dt$xcorr_tlag_max[which(sapply(1:length(transcriptome_proteome_contrasts_1dt$xcorr_tlag_max), function(xc) length(transcriptome_proteome_contrasts_1dt$xcorr_tlag_max[xc][[1]]) < 1))] <- NA
transcriptome_proteome_contrasts_1dt$xcorr_tlag_max <- unlist(transcriptome_proteome_contrasts_1dt$xcorr_tlag_max)

cor(transcriptome_proteome_contrasts_1dt$xcorr_tlag_max, transcriptome_proteome_contrasts_1dt$corr, use = "com")
sapply(sort(unique(transcriptome_proteome_contrasts_1dt$xcorr_tlag_max)), function(lg){
  y <- (transcriptome_proteome_contrasts_1dt$corr[transcriptome_proteome_contrasts_1dt$xcorr_tlag_max == lg]);
  c(lag = lg, mean = mean(y, na.rm = T), var = var(y, na.rm = T))})


# transcriptome_proteome_contrasts_1dt[(transcriptome_proteome_contrasts_1dt$trans_cluster) == 1,] #hallelujah it is consistent

#between sex, within tissue
sex_matched_corrs <- do.call(rbind, lapply(tissues, function(tissuei){
  
  ms <- transcriptome_proteome_contrasts_1dt[(transcriptome_proteome_contrasts_1dt$tissue) == tissuei & (transcriptome_proteome_contrasts_1dt$sex) == "male",] 
  fs <- transcriptome_proteome_contrasts_1dt[(transcriptome_proteome_contrasts_1dt$tissue) == tissuei & (transcriptome_proteome_contrasts_1dt$sex) == "female",] 
  
  inters_genes <- intersect(ms$feature_ID, fs$feature_ID)
  ms <- ms[ms$feature_ID %in% inters_genes,]
  fs <- fs[fs$feature_ID %in% inters_genes,]
  
  ms <- ms[order(ms$feature_ID),]
  fs <- fs[order(fs$feature_ID),]
  
  return(data.frame(ms$tissue, ms$trans_cluster, ms$gene_id, fs$gene_id, ms$corr, fs$corr))
  
}))

all.equal(sex_matched_corrs[,1], sex_matched_corrs[,2])
cor(sex_matched_corrs[,c("ms.corr", "fs.corr")], use = "comp")
plot(sex_matched_corrs[,c("ms.corr", "fs.corr")], col = adjustcolor(1, 0.1), pch = 19)
hist(sex_matched_corrs$ms.corr)
hist(sex_matched_corrs$fs.corr)

#between tissue, between sex
common_genes <- do.call(c, lapply(tissues, function(tissuei) lapply(c("male", "female"), function(sexj)
  c(transcriptome_proteome_contrasts_1dt[(transcriptome_proteome_contrasts_1dt$tissue) == tissuei & (transcriptome_proteome_contrasts_1dt$sex) == sexj,"gene_id"])
  )))
common_genes <- lapply(common_genes, function(cg) cg[[1]])
common_genes <- sort(Reduce(intersect, common_genes))
tissue_sex_corrs <- do.call(cbind, lapply(tissues, function(tissuei){ sapply(c("male", "female"), function(sexj){
  sub <- transcriptome_proteome_contrasts_1dt[(transcriptome_proteome_contrasts_1dt$tissue) == tissuei & 
                                               (transcriptome_proteome_contrasts_1dt$sex) == sexj &
                                               transcriptome_proteome_contrasts_1dt$gene_id %in% common_genes,]
  sub <- sub[!duplicated(sub$gene_id),]
  return(sub$corr)
})}))

colnames(tissue_sex_corrs) <- paste0(rep(tissues, each = 2), "-", rep(c("male", "female"), 4))
str(tissue_sex_corrs)
cor(tissue_sex_corrs)
apply(tissue_sex_corrs, 2, mean)

#####

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

## put (absolute) correlations on the upper panels,
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

png(filename = "Documents/transcr_v_protein/pairs_plot.png", width = 1000, height = 1000)
pairs(tissue_sex_corrs, diag.panel = panel.hist, lower.panel = panel.cor, pch = 19, col = adjustcolor(1, 0.1), cex.cor = 80, prefix = "r =\n")
dev.off()
