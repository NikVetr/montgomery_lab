#load libraries
library(ks)
library(data.table)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
library(MotrpacBicQC)
library(plotrix)
library(ggplot2)
library(testit)
library(circlize)
library(jpeg)
library(foreach)
library(doParallel)
library(pracma)


#### define functions ####
source(file = "~/scripts/montgomery_lab/deg-trait_functions.R")

#general graphical parameters
motrpac_gtex_map = c('t30-blood-rna'='Whole Blood',
                     't52-hippocampus'='Brain - Hippocampus',
                     't53-cortex'='Brain - Cortex',
                     't54-hypothalamus'='Brain - Hypothalamus',
                     't55-gastrocnemius'='Muscle - Skeletal',
                     't56-vastus-lateralis'='Muscle - Skeletal',
                     't58-heart'='Heart - Left Ventricle',
                     't59-kidney'='Kidney - Cortex',
                     't60-adrenal'='Adrenal Gland',
                     't61-colon'='Colon - Transverse',
                     't62-spleen'='Spleen',
                     't63-testes'='Testis',
                     't64-ovaries'='Ovary',
                     't66-lung'='Lung',
                     't67-small-intestine'='Small Intestine - Terminal Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose - Subcutaneous')

cols = list(Tissue=tissue_cols[names(motrpac_gtex_map)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'

#### figure 1 ####
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
coloc_phenotypes <- stringr::str_replace_all(gwas_names, "imputed_", "")
ldsc_output_dir <- "~/repos/ldsc/output/"
ldsc_results_paths <- list.files(ldsc_output_dir)
ldsc_log_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "log")])
ldsc_results_paths <- paste0(ldsc_output_dir, ldsc_results_paths[grep(ldsc_results_paths, pattern = "results")])
cluster_names <- paste0("Cluster_", 1:15)
ldsc_results_paths <- lapply(gwas_names, function(gwas) ldsc_results_paths[grep(pattern = gwas, ldsc_results_paths)])
ldsc_log_paths <- lapply(gwas_names, function(gwas) ldsc_log_paths[grep(pattern = gwas, ldsc_log_paths)])
names(ldsc_results_paths) <- names(ldsc_log_paths) <- gwas_names


log_files <- as.data.frame(matrix(data = NA, ncol = 5, nrow = length(gwas_names) * length(cluster_names)))
colnames(log_files) <- c("h2", "h2se", "chi2", "gwas", "cluster")
log_files$cluster <- rep(1:length(cluster_names), length(gwas_names))
log_files$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
for(i in 1:nrow(log_files)){
  logfile <- readLines(ldsc_log_paths[[log_files$gwas[i]]][grep(pattern = paste0(cluster_names[log_files$cluster[i]], ".log"), ldsc_log_paths[[log_files$gwas[i]]])][1])
  h2 <- logfile[grep(logfile, pattern = "Total Observed scale h2")]
  log_files$h2[i] <- as.numeric(strsplit(h2, " ")[[1]][5])
  log_files$h2se[i] <- as.numeric(substr(x = strsplit(h2, " ")[[1]][6], start = 2, stop = nchar(strsplit(h2, " ")[[1]][6]) - 1))
  chi2 <- as.numeric(strsplit(logfile[grep(logfile, pattern = "Mean Chi")], split = " ")[[1]][3])
  log_files$chi2[i] <- chi2
}
log_files$gwas <- stringr::str_replace_all(log_files$gwas, "imputed_", "")

#subset to insteresting trait categories
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
log_files$trait_category <- traitwise_partitions$Category[match(log_files$gwas, traitwise_partitions$Tag)]
salient.categories <- c("Cardiometabolic", "Aging", "Anthropometric", "Immune", "Psychiatric-neurologic")
log_files <- log_files[log_files$trait_category %in% salient.categories,]

total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])

#do the gcor
ldsc_directory <- "~/repos/ldsc/"

#let's get the gwas summary stats into a format that munge_sumstats.py can process
traitnames <- traits_with_satisfactory_heritaility
gcors <- lapply(traitnames, function(x) integer())
names(gcors) <- traitnames  
for(gi in traitnames){
  logfile <- readLines(paste0(ldsc_directory, "output/signed_gcor/imputed_", gi, "_pairwise-Gcorrs.log"))
  logfile <- logfile[(grep(x = logfile, pattern = "Summary of Genetic Correlation Results")+1):(length(logfile)-3)]
  writeLines(logfile, con = paste0(ldsc_directory, "temp.txt"))
  gcors[[gi]] <- fread(paste0(ldsc_directory, "temp.txt"))
  file.remove(paste0(ldsc_directory, "temp.txt"))
  
  gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
  gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
  gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
  gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
}

#construct genetic correlation matrix
gcor_mat <- diag(length(traitnames))
colnames(gcor_mat) <- rownames(gcor_mat) <- traitnames
for(ri in rownames(gcor_mat)){
  for(ci in colnames(gcor_mat)){
    gcor_mat[ri, ci] <- gcors[[ri]]$rg[match(ci, gcors[[ri]]$p2)]  
  }
}
diag(gcor_mat) <- rep(1, length(traitnames))

#remove impossible correlations
gcor_mat <- gcor_mat[traits_with_satisfactory_heritaility,traits_with_satisfactory_heritaility]
troublesome_inds_NA_corrs <- as.numeric(names(which(table(which(is.na(gcor_mat), arr.ind = T)) > 1)))
troublesome_inds_NA_corrs_traits <- rownames(gcor_mat)[troublesome_inds_NA_corrs]
if(length(troublesome_inds_NA_corrs) > 0){
  gcor_mat <- gcor_mat[-troublesome_inds_NA_corrs, -troublesome_inds_NA_corrs]
}
troublesome_inds_imposs_corrs <- as.numeric(names(which(table(which(gcor_mat > 1 | gcor_mat < -1, arr.ind = T)) > 1)))
troublesome_inds_imposs_corrs_traits <- rownames(gcor_mat)[troublesome_inds_imposs_corrs]
if(length(troublesome_inds_imposs_corrs) > 0){
  gcor_mat <- gcor_mat[-troublesome_inds_imposs_corrs, -troublesome_inds_imposs_corrs]
}

rownames(gcor_mat) <- colnames(gcor_mat) <- traitname_map$new_Phenotype[match(rownames(gcor_mat), traitname_map$Tag)]

#perform hypothesis testing
gcor_mat_se <- diag(length(traitnames))
colnames(gcor_mat_se) <- rownames(gcor_mat_se) <- traitnames
for(ri in rownames(gcor_mat_se)){
  for(ci in colnames(gcor_mat_se)){
    gcor_mat_se[ri, ci] <- gcors[[ri]]$se[match(ci, gcors[[ri]]$p2)]
  }
}
diag(gcor_mat_se) <- rep(NA, length(traitnames))
rownames(gcor_mat_se) <- colnames(gcor_mat_se) <- traitname_map$new_Phenotype[match(rownames(gcor_mat_se), traitname_map$Tag)]
gcor_mat_se <- gcor_mat_se[rownames(gcor_mat), rownames(gcor_mat)]
gcor_mat_z <- abs(gcor_mat / gcor_mat_se)
zscore_sig_thresh <- abs(qnorm(0.025  / choose(dim(gcor_mat)[1], 2), 0, 1))
gcor_mat_sig <- gcor_mat_z > zscore_sig_thresh
diag(gcor_mat_sig) <- T
sum(gcor_mat_sig, na.rm = T) - nrow(gcor_mat_sig)

#get a few last graphical parameters
unique_trait_categories <- unique(traitwise_partitions$Category)
category_colors <- RColorBrewer::brewer.pal(length(unique_trait_categories), "Dark2")
names(category_colors) <- sort(salient.categories)
cols$category <- category_colors
# cols$category <- cols$Tissue[1:length(unique_trait_categories)+1]
# names(cols$category) <- unique_trait_categories

#### actually plot genetic correlations ####

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/figure1_high-level-overview.pdf"), 
                     width = 1300 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mar = c(6,5,3,8), mfrow = c(1,2), xpd = NA)
# par(mar = c(6,5,3,4), xpd = NA)
# layout(mat = t(as.matrix(c(rep(1,5), rep(2,6)))))


#plotting params
traits <- rownames(gcor_mat)
rate = 0.001
exp_dist_cols <- round(cumsum(c(1, dexp(1:100, rate = rate) / min(dexp(1:100, rate = rate)))))
heatcols <- viridis::magma(max(exp_dist_cols))[exp_dist_cols]
heatcols <- RColorBrewer::brewer.pal(11, "RdBu")[-c(1,11)]
heatcols <- rev(colorRampPalette(heatcols)(101))
# plot(1:length(heatcols), 1:length(heatcols), cex = 2, pch = 19, col = heatcols)

plot(1,1,xlim = c(0,length(traits)), ylim = c(0,length(traits)), xpd = NA,
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sorted_row_inds <- order(cmdscale(1-gcor_mat, k = 1))
for(rowi in 1:length(traits)){
  text(labels = traits[sorted_row_inds[rowi]], x = length(traits), y = length(traits) - rowi + 1, 
       col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
         match(traits[sorted_row_inds[rowi]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
       , pos = 4, xpd = NA, cex = 0.75, font = 1)
  for(colj in 1:length(traits)){
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = colj - 1/2, xr = colj + 1/2, pch = 15, cex = 1, border = NA,
         col = heatcols[round((gcor_mat[sorted_row_inds[rowi], sorted_row_inds[colj]] + 1) / 2 * 100) + 1])
    
    if(gcor_mat_sig[sorted_row_inds[rowi], sorted_row_inds[colj]]){
      points(y = length(traits) - rowi + 1, x = colj, pch = 19, col = "white", cex = 0.4)
    }
    
    if(rowi == 1){
      text(labels = traits[sorted_row_inds[colj]], x = colj + 1, y = 0.2, 
           col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
             match(traits[sorted_row_inds[colj]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
           , pos = 2, srt = 45, xpd = NA, cex = 0.6, font = 1)
    }
    
    
    
    
    
  }
}

xl = -2.5; xr = -0.5; yb = length(traitnames) / 3; yt = length(traitnames) / 1.5
ydisp <- -3
rect(xleft = xl, xright = xr, col = heatcols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)] + ydisp, yt= seq(yb, yt, length.out = length(heatcols)+1)[-1] + ydisp)
rect(xleft = xl, xright = xr, ybottom = yb + ydisp, ytop = yt + ydisp)
text(labels = -5:5/5, x = xl, pos = 2, y = seq(yb, yt, length.out = 11) + ydisp, cex = 1)
text(x = mean(c(xl, xr)), y = yt + ydisp - 0.25, labels = latex2exp::TeX("$r_g$"), pos = 3, cex = 2, font = 2)
text("Genetic Correlation Matrix", x = length(traits) / 2 + 0.5, y = length(traits) + 0.5, pos = 3, cex = 2.5, font = 2)

points(rep((xl+xr)/2 + 0.5, length(salient.categories)), 1:length(salient.categories)*1.05, pch = 15, col = cols$category[salient.categories], cex = 1.75)
text(rep((xl+xr)/2 + 0.5, length(salient.categories)), 1:length(salient.categories)*1.05, pos = 2, col = 1, 
     labels = gsub("-.*", "", salient.categories), cex = 0.75)

text(labels = latex2exp::TeX(paste0("• mark p-val < $10^{", round(log10(0.025  / choose(dim(gcor_mat)[1], 2)), 2), "}$")), 
     x = -0.5, y = yt + ydisp + nrow(gcor_mat)/5, cex = 0.75, srt = 90)
text(labels = "bonferroni", x = -1.5, y = yt + ydisp + nrow(gcor_mat)/5, cex = 0.75, srt = 90, font = 3, family = "Arial")

fig_label(text = "a)", region = "plot", cex = 3, shrinkX = -20, shrinkY = 1.065)


#now do the tissue transcriptome correlation matrix
zcor = as.matrix(read.table("~/data/smontgom/zcor_transcriptome_pass1b.tsv"))
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

#figure parameters
nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])
nice_names <- sapply(names(colours$Tissue), function(tissue) nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)])

corr_thresh <- 0.3

axis.length <- 1.5
center_rescaler <- 1.25
inner_shifter <- 0.96
no_concentric_arcs <- F
opacity_concentric_arcs <- 1
opacity_nonconcentric_arcs <- 0.5
outer_shifter <- 1.04
line_weight_power <- 2
line_weight_multiplier <- 15
numbers_in_squares <- F
tissue_names_not_colors <- T
tissue_name_cex <- 0.29
across_relationships <- T
adjacent_relationships <- T
opacity_multiplier_within_tissue <- 0.35

tissues <- as.list(names(colours$Tissue))
tissues <- unlist(tissues)


tissues_to_include <- tissues
print(tissues_to_include)

par(mar = c(7,10,5.5,6), xpd = NA)

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
et$opacity[et$tiss1 == et$tiss2] <- et$opacity[et$tiss1 == et$tiss2] * opacity_multiplier_within_tissue
et$opacity[et$concentric] <- opacity_concentric_arcs

et_wt <- lapply(1:4, function(time) et[et$time1 == time & et$time2 == time & (et$tiss1 %in% tissues_to_include | et$tiss2 %in% tissues_to_include),])
# et_wt <- lapply(et_wt, function(et_sub) et_sub[1:3,])
if(no_concentric_arcs){
  et_wt <- lapply(1:4, function(time) et_wt[[time]][et_wt[[time]]$tiss1 != et_wt[[time]]$tiss2,])
}

#get within-tissue opacities multiplied
for(i in 1:length(et_wt)){
  et_wt[[i]]$opacity[et_wt[[i]]$tiss1 == et_wt[[i]]$tiss2] <- et_wt[[i]]$opacity[et_wt[[i]]$tiss1 == et_wt[[i]]$tiss2] * opacity_multiplier_within_tissue
}


plot(1,1,xlim = c(-2,2), ylim = c(-2,2), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
centers <- list(c(-1,1), c(1,1), c(1,-1), c(-1,-1))
centers <- lapply(centers, function(x) x * center_rescaler)
for(time in 1:4){
  if(nrow(et_wt[[time]]) == 0){next()}
  for(ri in 1:nrow(et_wt[[time]])){
    arc(t1 = aas[[time]][et_wt[[time]][ri,"theta1"]], t2 = aas[[time]][et_wt[[time]][ri,"theta2"]], r1 = et_wt[[time]][ri,"r1"], center = centers[[time]] * outer_shifter,
        r2 = et_wt[[time]][ri,"r2"], lwd = et_wt[[time]]$weight[ri]^line_weight_power * line_weight_multiplier, col = adjustcolor(et_wt[[time]]$color[ri], et_wt[[time]]$opacity[ri]), 
        random_selfing = F, clockwise_selfing = (et_wt[[time]][ri,"sex1"] == "1") * c(1,-1,1,-1)[time]*-1 + rev(c(0,1,0,1))[time], self_adjust = 0.3)
  }
  for(sex in 1:2){
    # line(t = aas[[time]][sex], 0, axis.length, col = colours$Sex[sex], lwd = 4, center = centers[[time]])
    # for(tissue in 1:n_tiss){
    #   polarp(t = aas[[time]][sex], r = rs[tissue], col = colours$Tissue[tissue], pch = 19, cex = 0.85, center = centers[[time]])
    # }
  }
}

max_outer_radii <- sapply(et_wt, function(etwt) max(etwt$r1 + etwt$r2) / 2)
week_label_locs <- sapply(1:4, function(w) polar2cart(r = max_outer_radii[w], t = (-(0:3)*pi/2 + 3*pi/4)[w]) + centers[[w]])
week_label_locs[week_label_locs == Inf] <- -center_rescaler * 1.2
week_label_locs[week_label_locs == -Inf] <- center_rescaler * 1.2
timelab_nudge <- 0.25

for(i in 1:4){
  text(c("\u2642", "\u2640", "\u2642", "\u2640")[i], col = c(colours$Sex, colours$Sex)[i], 
       x = c(0,1,0,-1)[i]*(2 - sqrt(2) + axis.length)*1.15 - c(2,0,-2.15,0)[i] * 0.06, 
       y = c(1,0,-1,0)[i]*(2 - sqrt(2) + axis.length)*1.275 - c(-1,-1,-1,-1)[i] * 0.06, 
       cex = 4.5, srt = c(45,90,225,270)[i], pos = 1)
  # shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = c(-1,1,1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, 
  #            y = c(1,1,-1,-1)[i]*(1 + axis.length/sqrt(2)) * 0.99 * center_rescaler, cex = 4,
  #            srt = c(45,-45,235,135)[i], pos = 1)
  
  #plot weeks
  shadowtext(paste0(c(1,2,4,8), "w")[i], col = colours$Time[i], x = week_label_locs[1,i] + c(-1,1,1,-1)[i]*timelab_nudge, 
             y = week_label_locs[2,i] + c(1,1,-1,-1)[i]*timelab_nudge, cex = 4,
             srt = c(45,-45,235,135)[i], pos = 1)
  
}

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
# rs_bt <- c(rs, rs + max(rs) + diff(rs)[1]) - diff(rs)[1] + sqrt(2) * center_rescaler - axis_length
# polarp(r = rs_bt, t = rep(pi/4, length(rs_bt)), center = c(-2.5,-0), cex = , pch = 1, col = "black")



# rs_bt <- seq(sqrt(2) * center_rescaler - max(rs), 
#              sqrt(2) * center_rescaler, 
#              length.out = n_tiss)
# rs_bt <- c(rs_bt, rs_bt + max(rs_bt) - min(rs_bt) + diff(rs_bt)[1])

# rs_bt <- rs_bt - rs_bt[1] + diff(rs_bt)[1]

et$tpairs <- sapply(1:nrow(et), function(ri) paste0(sort(c(et$time1[ri], et$time2[ri])), collapse = ""))
et$close_sex <- NA
et$close_sex[et$tpairs == "12"] <- 1
et$close_sex[et$tpairs == "23"] <- 2
et$close_sex[et$tpairs == "34"] <- 1
et$close_sex[et$tpairs == "14"] <- 2

# et$r1 <- rs_bt[match(et$tiss1, names(colours$Tissue)) + (as.numeric(et$sex1) - 1) * n_tiss]
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

#get between-tissue opacities multiplied
for(i in 1:length(et_bt)){
  et_bt[[i]]$opacity[et_bt[[i]]$tiss1 == et_bt[[i]]$tiss2] <- et_bt[[i]]$opacity[et_bt[[i]]$tiss1 == et_bt[[i]]$tiss2] * opacity_multiplier_within_tissue
}

#hack to fix sex mixup -- TODO find originl bug
# for(time in 1:4){
#   temp <- et_bt[[time]]$sex1
#   et_bt[[time]]$sex1 <- et_bt[[time]]$sex2
#   et_bt[[time]]$sex2 <- temp
# }

# et_bt <- lapply(et_bt, function(et_sub) et_sub[1:10,])
# et_bt <- lapply(1:4, function(et_sub) et_bt[[et_sub]])

if(adjacent_relationships){
  
  for(time in 1:4){
    if(nrow(et_bt[[time]]) == 0){next()}
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




# time = 1
# ri = 1
# arc(t1 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta1"]], 
#     t2 = aas_bt[[time]][[et_bt[[time]][ri,"tiw"]]][et_bt[[time]][ri,"theta2"]], 
#     r1 = et_bt[[time]][ri,"r1"], 
#     r2 = et_bt[[time]][ri,"r2"],
#     center = centers_bt[[time]][et_bt[[time]][ri,"tiw"]][[1]],
#     lwd = et_bt[[time]]$weight[ri]^1.5 * 5, 
#     col = adjustcolor(et_bt[[time]]$color[ri], 0.85)
# )

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
      for(tissue in 1:n_tiss){
        cartesian_coordinates <- polar2cart(aas[[time]][sex], rs[tissue]) + centers[[time]]
        if((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)){
          name_color = colours$Sex[sex]
        } else {
          name_color = colours$Sex[sex]
        }
        text(nice_names[tissue], x = cartesian_coordinates[1], y = cartesian_coordinates[2], 
             col = adjustcolor(name_color, ifelse(tissues[tissue] %in% tissues_to_include & 
                                                    !((nice_names[tissue] == "TESTES" & sex == 2) | (nice_names[tissue] == "OVARY" & sex == 1)), 1, 0.5)), 
             srt = c(-45,45,-45,45)[time], adj = 0.5, cex = 1 / strwidth(nice_names[tissue]) * tissue_name_cex, font = 2)   
      }
    }
  }
}

tissues <- names(nice_names)
if(across_relationships){
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
      # cartesian_coordinates_1 <- (cartesian_coordinates_1 - centers[as.numeric(et_bt_across$time1[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time1[ri])][[1]]
      # cartesian_coordinates_2 <- (cartesian_coordinates_2 - centers[as.numeric(et_bt_across$time2[ri])][[1]]) * inner_shifter + centers[as.numeric(et_bt_across$time2[ri])][[1]]
      # cartesian_coordinates_1 <- cartesian_coordinates_1 * inner_shifter
      # cartesian_coordinates_2 <- cartesian_coordinates_2 * inner_shifter
      segments(x0 = cartesian_coordinates_1[1], y0 = cartesian_coordinates_1[2], x1 = cartesian_coordinates_2[1], y1 = cartesian_coordinates_2[2],
               lwd = et_bt_across$weight[ri]^line_weight_power * line_weight_multiplier,
               col = adjustcolor(et_bt_across$color[ri], et_bt_across$opacity[ri]))
    }
  }
}

xl <- -3.2; yb <- -0.425; xr <- -2.9; yt <- 0.575;
# et_wt[[time]]$weight[ri]^1 * 4
line_weights <- c(0:10/10)^line_weight_power * line_weight_multiplier
lws <- line_weights / 96 / (par("pin")[1]  / 4)
# segments(xl,yb,xl,yt,lwd = 10)
corresponding_heights <- seq(0,abs(yt - yb)/2,length.out = 11)


# line_weight_power <- 2
# line_weight_multiplier <- 10

hoffset <- 0 + c(5:0/5, 1:5/5)^line_weight_power * line_weight_multiplier / 300
voffset_rhos <- -0.1
text(labels = round(seq(-1, 1, length.out = 11), 2), y = seq(yb, yt, length.out = 11) + voffset_rhos, 
     x = xl - 0.0075 + hoffset, pos = 4, las=2, cex=0.9)
text(labels = latex2exp::TeX(paste0("|$\\rho$|")), y = yt - 0.03 + voffset_rhos, x = (xl) - (xr-xl)*0.25, pos = 3, cex = 2)
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

fig_label(text = "b)", region = "plot", cex = 3, shrinkX = 1.5, shrinkY = 1.32)

dev.off()



#### figure 2: relative expression cdf ####
load("~/data/smontgom/relative_expression_motrpac_gtex")
nnmap <- as.data.frame(bic_animal_tissue_code[,4:5])

qs2use <- 1:9999/10000
qs2use <- c(1/10000, 1:99/100, 9999/10000)
qs2use <- round(invlogit(seq(-9,9, length.out = 75)), 4)

grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/motrpac_companion/figures/figure2_relative-expression_quantile-function.pdf"), 
                     width = 1200 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mfrow = c(2,4), xpd = NA, mar = c(6,4.5,4,1))
for(type in c(3,4)){
  
  EZ_PZ <- relative_expression_data[[type]]
  EZ_PZ <- lapply(EZ_PZ, function(x) x[match(paste0(sprintf("%.2f", qs2use*100), "%"), paste0(sprintf("%.2f", as.numeric(gsub("%", "", rownames(x)))), "%")),])
  
  for(ti in paste0(2^(0:3), "w")){
    f_p <- 0.4
    f_x <- 2
    ylims = c(min(sort(EZ_PZ[[ti]])[sort(EZ_PZ[[ti]]) != -Inf]),max(sort(EZ_PZ[[ti]],T)[sort(EZ_PZ[[ti]],T) != Inf]))
    ylims <- squish_middle_x(ylims, f_x)
    
    plot(100,100,xlim = c(0,1.25), ylim = ylims, xpd=NA, ylab = "Standardized Effect Size (SD)", xlab = "", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)
    text("Quantile", x = 0.5, y = ylims[1] - diff(ylims)/6, pos = 1, cex = 1.25)
    if(ti == "2w"){
      text(latex2exp::TeX(paste0("Ratio of Exercise DE to \\sqrt{", ifelse(type == 1, "Phenotypic", "Genetic"), " Variance in $log_2$(Gene Expression)}")), x = 1.35, y = ylims[2] + diff(ylims) / 15, pos = 3, cex = 3)}
    if(ti == "1w"){
      fig_label(ifelse(type == 1, "a)", "b)"), region = "figure", pos = "topleft", cex = 3.5)
      }
    
    xylocs_tissue_names <- cbind(rep(1.05, ncol(EZ_PZ[[ti]])), redistribute(as.numeric(squish_middle_x(tail(EZ_PZ[[ti]], 1), f_x)), diff(ylims) / 24))
    rownames(xylocs_tissue_names) <- colnames(EZ_PZ[[ti]]); colnames(xylocs_tissue_names) <- c("x", "y")
    
    #horiz axis
    segments(x0 = 0, y0 = ylims[1] - diff(ylims) / 100, x1 = 1, y1 = ylims[1] - diff(ylims) / 100, lwd = 2)
    segments(x0 = 0:10/10, y0 = ylims[1] - diff(ylims) / 100, x1 = 0:10/10, y1 = ylims[1] - diff(ylims) / 50, lwd = 2, xpd = NA)
    horiz_axis_labels <- round(unsquish_middle_p(0:10/10, f_p), 3)
    horiz_axis_labels[1] <- 0; horiz_axis_labels[length(horiz_axis_labels)] <- 1;
    text(labels = horiz_axis_labels, x = 0:10/10 + 0.035, y = rep(ylims[1] - diff(ylims) / 25, 10), pos = 2, xpd = NA, srt = 90)
    segments(x0 = 0.5, y0 = ylims[1], x1 = 0.5, y1 = ylims[2], lwd = 2, lty = 2, col = "grey50")
    segments(x0 = 0:10/10, y0 = ylims[1], x1 = 0:10/10, y1 = ylims[2], lwd = 1, lty = 3, col = "grey75")
    
    #vert axis
    segments(x0 = -1/100, y0 = ylims[1] + diff(ylims)/100, x1 = -1/100, y1 = ylims[2], lwd = 2)
    segments(x0 = -1/100, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = -1/50,
             y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 2)
    text(labels = round(unsquish_middle_x(seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), f_x), 2),
         x = -1/50, y = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), pos = 2, xpd = NA)
    segments(x0 = 0, y0 = 0, x1 = 1, y1 = 0, lwd = 2, lty = 2, col = "grey50")
    segments(x0 = 0, y0 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), x1 = 1, 
             y1 = seq(ylims[1] + diff(ylims)/100, ylims[2], length.out = 10), lwd = 1, lty = 3, col = "grey75")
    
    
    for(tissue in colnames(EZ_PZ[[ti]])){
      lines(squish_middle_p(qs2use, f_p), squish_middle_x(EZ_PZ[[ti]][,tissue], f_x), col = cols$Tissue[tissue])
      text(nnmap$abbreviation[match(tissue, nnmap$tissue_name_release)], x = xylocs_tissue_names[tissue,"x"], y = xylocs_tissue_names[tissue,"y"], pos = 4, xpd = NA,
           col = cols$Tissue[tissue])
      segments(x0 = squish_middle_p(tail(qs2use, 1), f_p), x1 = xylocs_tissue_names[tissue,"x"]+0.01,
               y0 = squish_middle_x(tail(EZ_PZ[[ti]][,tissue], 1), f_x), y1 = xylocs_tissue_names[tissue,"y"],
               lty = 3, col = cols$Tissue[tissue])
    }
    shadowtext(x = -0.01, y = ylims[2] - diff(ylims) / 8, labels = ti, 
               cex = 5, col = cols$Time[ti], pos = 4, r = 0.2)
  }
}

dev.off()