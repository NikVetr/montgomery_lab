# knitr::opts_knit$set(root.dir = '/oak/stanford/groups/smontgom/nicolerg/MOTRPAC/PASS_ANALYSIS/PASS1B/RNA/eqtl-coloc')

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


#### define functions ####

plotMatrix <- function(mobject, size, location, lwd = 2, lwd_inner = 1.5, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])), lwd = lwd_inner)
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])), lwd = lwd_inner)
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}

draw.contour<-function(a,alpha=0.95,plot.dens=FALSE, line.width=2, line.type=1, limits=NULL, density.res=300,spline.smooth=-1,...){
  ##a is a list or matrix of x and y coordinates (e.g., a=list("x"=rnorm(100),"y"=rnorm(100)))
  ## if a is a list or dataframe, the components must be labeled "x" and "y"
  ## if a is a matrix, the first column is assumed to be x, the second y
  ##alpha is the contour level desired
  ##if plot.dens==TRUE, then the joint density of x and y are plotted,
  ##   otherwise the contour is added to the current plot.
  ##density.res controls the resolution of the density plot
  
  ##A key assumption of this function is that very little probability mass lies outside the limits of
  ## the x and y values in "a". This is likely reasonable if the number of observations in a is large.
  
  require(MASS)
  require(ks)
  if(length(line.width)!=length(alpha)){
    line.width <- rep(line.width[1],length(alpha))
  }
  
  if(length(line.type)!=length(alpha)){
    line.type <- rep(line.type[1],length(alpha))
  }
  
  if(is.matrix(a)){
    a=list("x"=a[,1],"y"=a[,2])
  }
  ##generate approximate density values
  if(is.null(limits)){
    limits=c(range(a$x),range(a$y))
  }
  f1<-kde2d(a$x,a$y,n=density.res,lims=limits)
  
  ##plot empirical density
  if(plot.dens) image(f1,...)
  
  if(is.null(dev.list())){
    ##ensure that there is a window in which to draw the contour
    plot(a,type="n",xlim=limits[1:2],ylim=limits[3:4],...)
  }
  
  ##estimate critical contour value
  ## assume that density outside of plot is very small
  
  zdens <- rev(sort(f1$z))
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)])
  for(cont.level in 1:length(alpha)){
    ##This loop allows for multiple contour levels
    crit.val <- zdens[max(which(Czdens<=alpha[cont.level]))]
    
    ##determine coordinates of critical contour
    b.full=contourLines(f1,levels=crit.val)
    for(c in 1:length(b.full)){
      ##This loop is used in case the density is multimodal or if the desired contour
      ##  extends outside the plotting region
      b=list("x"=as.vector(unlist(b.full[[c]][2])),"y"=as.vector(unlist(b.full[[c]][3])))
      
      ##plot desired contour
      line.dat<-xspline(b,shape=spline.smooth,open=TRUE,draw=FALSE)
      lines(line.dat,lty=line.type[cont.level],lwd=line.width[cont.level])
    }
  }
}

minwhich <- function(x){min(which(x))}
maxwhich <- function(x){max(which(x))}

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

#gsutil = '~/google-cloud-sdk/bin/gsutil'
#source('/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/integrative/outliers_covariates/pi1_cook_fx.R')
#source('/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/tools/get_fx.R')

#### load data ####

# get tissue labels 
# load MotrpacBicQC conversion tables 
bic_animal_tissue_code = as.data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code = bic_animal_tissue_code[tissue_name_release!='']
bic_animal_tissue_code[,my_tissue := tolower(gsub(" Powder", "", bic_tissue_name))]
bic_animal_tissue_code[,my_tissue := gsub(' ','_',my_tissue)]
bic_animal_tissue_code[my_tissue == 'blood_rna', my_tissue := 'paxgene_rna']
tissue_codes = bic_animal_tissue_code[, tissue_name_release]
tissue_codes = tissue_codes[tissue_codes != '']

## Define data paths

#download data to local drive
# system("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/motrpac/mawg_data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt /Users/nikolai/data/smontgom/")

GTEx_logo <- readJPEG("~/Documents/GTEx_logo.jpg")
eqtl = '~/data/smontgom/GTEx_Analysis_v8_eQTL'
#deg = 'gs://mawg-data/pass1b-06/transcript-rna-seq/dea/'
deg = '~/data/smontgom/dea/'
#map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t')
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t', header=T)
map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
map <- unique(map[,c("feature_ID", "human_ensembl_gene", "human_gene_symbol")])
gwas = '~/data/smontgom/imputed_gwas_hg38_1.1'
coloc = '~/data/smontgom/results_enloc_priors'

motrpac_gtex_map = c('t30-blood-rna'='Whole_Blood',
                     't52-hippocampus'='Brain_Hippocampus',
                     't53-cortex'='Brain_Cortex',
                     't54-hypothalamus'='Brain_Hypothalamus',
                     't55-gastrocnemius'='Muscle_Skeletal',
                     't56-vastus-lateralis'='Muscle_Skeletal',
                     't58-heart'='Heart_Left_Ventricle',
                     't59-kidney'='Kidney_Cortex',
                     't60-adrenal'='Adrenal_Gland',
                     't61-colon'='Colon_Transverse',
                     't62-spleen'='Spleen',
                     't63-testes'='Testis',
                     't64-ovaries'='Ovary',
                     't66-lung'='Lung',
                     't67-small-intestine'='Small_Intestine_Terminal_Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose_Subcutaneous')

## Download and format RNA-seq differential expression analysis results for easy loading
# timewise_dea_list = list()
# training_dea_list = list()
# # for (file in system(sprintf('gsutil ls %s | grep "transcript-rna-seq_dea_20201028.txt"', deg), intern=T)){
# #   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
# #   dea_list[[tissue_code]] = dl_read_gcp(file, sep='\t')
# # }
# for (file in list.files(deg, full.names=T, pattern='timewise-dea')){
#   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
#   timewise_dea_list[[tissue_code]] = fread(file, sep='\t', header=T)
# }
# for (file in list.files(deg, full.names=T, pattern='training-dea')){
#   tissue_code = gsub('_.*','',gsub('.*pass1b-06_','',file))
#   training_dea_list[[tissue_code]] = fread(file, sep='\t', header=T)
# }
# training_dea_list = lapply(training_dea_list, function(x) x[,removed_samples := as.character(removed_samples)])
# timewise_dea_list = lapply(timewise_dea_list, function(x) x[,removed_samples := as.character(removed_samples)])
# training_dea = rbindlist(training_dea_list)
# timewise_dea = rbindlist(timewise_dea_list)
# # remove aorta (BAT contamination in F 1w,2w samples)
# training_dea = training_dea[tissue!='t65-aorta']
# timewise_dea = timewise_dea[tissue!='t65-aorta']
# # adjust tests across all tissues
# training_dea[,adj_p_value := p.adjust(p_value, method='BH')]
# table(training_dea[,adj_p_value] < 0.05, training_dea[,tissue])
# rna_dea = list(training_dea = training_dea,
#                timewise_dea = timewise_dea)
# save(rna_dea, file='/oak/stanford/groups/smontgom/shared/motrpac/shared_rdata/rna_dea_20210114.RData')

if(!file.exists("~/data/smontgom/rna_dea_20210114.RData")){
  system("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/motrpac/shared_rdata/rna_dea_20210114.RData /Users/nikolai/data/smontgom/",
       show.output.on.console = TRUE)
}

if(!exists("rna_dea")){
  load('~/data/smontgom/rna_dea_20210114.RData')
  load('~/data/smontgom/dea/transcript_rna_seq_20210804.RData')
  rna_dea$training_dea <- as.data.table(transcript_rna_seq$training_dea)
  rna_dea$timewise_dea <- as.data.table(transcript_rna_seq$timewise_dea)
  rm(transcript_rna_seq)
}
# list with two data.tables:
# training_dea: use `adj_p_value` to select genes that are differential due to training (use 0.1 (10% FDR) threshold for now)
# timewise_dea: cross-reference list of DEGs with this table to ID sex- and time- specific log fold-changes 

## Intersect DA results with GTEx eQTLs 
if(!exists("deg_eqtl_list")){
  deg_eqtl_list = list()
  for(motrpac_tissue in unique(rna_dea$timewise_dea$tissue)){
    
    if(!motrpac_tissue %in% names(motrpac_gtex_map)){next}
    
    cat(paste0(motrpac_tissue, "\n"))
    
    # read in eQTLs
    gtex_tissue = motrpac_gtex_map[[motrpac_tissue]]
    gtex_egene = fread(sprintf('%s/%s.v8.egenes.txt.gz',eqtl,gtex_tissue), sep='\t', header=T)
    gtex_egene[, human_ensembl_gene := gsub('\\..*','',gene_id)]
    gtex_egene$human_gene_symbol <- gtex_egene$gene_name
    
    
    
    # match human genes with rat ensembl genes 
    gtex_egene = merge(gtex_egene, map, by='human_gene_symbol')
    gtex_motrpac = merge(gtex_egene, rna_dea$timewise_dea[tissue == motrpac_tissue], by='feature_ID')
    gtex_motrpac$abs_slope <- abs(gtex_motrpac$slope)
    
    cat(paste0("prop of feature_IDs matched: ", round(length(unique(gtex_motrpac$feature_ID)) / 
                   length(unique(rna_dea$timewise_dea[tissue == motrpac_tissue]$feature_ID)), 3),
                 "\nprop gene symbols in map: ", round(mean(gtex_egene$human_gene_symbol %in% map$human_gene_symbol), 2), "\n"))
    
    deg_eqtl_list[[motrpac_tissue]] = gtex_motrpac
    # gtex_motrpac$gene_symbol
    
    # de-duplicate gene entries, slowly but straightforwardly, while also preserving ensemble IDs
    # dup_gene_IDs <- table(gtex_motrpac$gene_id)
    # dup_gene_IDs <- names(dup_gene_IDs)[dup_gene_IDs > ifelse(any(motrpac_tissue == c("t63-testes", "t64-ovaries")), 4, 8)]
    # for(dgid in dup_gene_IDs){
    #   dgid_inds <- which(gtex_motrpac$gene_id == dgid)
    #   # print(length(unique(gtex_motrpac[dgid_inds,]$slope)))
    #   # print(length(unique(gtex_motrpac[dgid_inds,]$logFC)))
    #   a <- deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == dgid,]
    #   print(paste0(dgid, " has a max unique entry count of ", max(sapply(colnames(a), function(i) length(unique(a[[i]]))))))
    # }
  }
}

# gene_ID_map <- unique.data.frame(as.data.frame(do.call(rbind, deg_eqtl_list)[,c("feature_ID","human_ensembl_gene")]))
# write.table(gene_ID_map, file = "~/data/smontgom/motrpac_geneID_map.txt", col.names = T)
# deg_eqtl_list[["t68-liver"]]
# for(i in 1:length(names(motrpac_gtex_map))){
#   print(nrow(deg_eqtl_list[[names(motrpac_gtex_map)[i]]]))
# }
# 
# length(unique(deg_eqtl_list[[motrpac_tissue]]$gene_id))
# length(deg_eqtl_list[[motrpac_tissue]]$gene_id)
# 
# deg_eqtl_list[[motrpac_tissue]]$logFC
# deg_eqtl_list[[motrpac_tissue]]$slope
# 
# ind1 <- sample(1:nrow(deg_eqtl_list[[motrpac_tissue]]), 1)
# ind2 <- sample(which(deg_eqtl_list[[motrpac_tissue]]$gene_id == deg_eqtl_list[[motrpac_tissue]]$gene_id[ind1]), 1)
# if(length(ind2) > 0){
#   colnames(deg_eqtl_list[[motrpac_tissue]])[which(abs(as.numeric(deg_eqtl_list[[motrpac_tissue]][ind1,]) -
#                                                       as.numeric(deg_eqtl_list[[motrpac_tissue]][ind2,])) > 1E-10)]
# }

#prop of duplicate gene names
# sapply(names(motrpac_gtex_map), function(name) length(deg_eqtl_list[[name]]$gene_id) / length(unique(deg_eqtl_list[[name]]$gene_id)))
# 
# sapply(names(motrpac_gtex_map), function(name) length(deg_eqtl_list[[name]]$gene_id) / length(unique(deg_eqtl_list[[name]]$gene_id)))
# 
# sapply(1:ncol(deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == "ENSG00000070087.13",]), 
#        function(x) length(unique(deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == "ENSG00000070087.13",x])))
# 
#some quick EDA
# counts <- table(deg_eqtl_list[[tss]]$gene_id)
# ind <- sample(which(counts > 8), 1)
# a <- deg_eqtl_list[[tss]][deg_eqtl_list[[tss]]$gene_id == names(ind),]
# max(sapply(colnames(a), function(i) length(unique(a[[i]]))))
# 
# for(i in colnames(a)){
# 
# }
# dev.off()
# par(mfrow = c(4,4))
# a = 0
# for(tss in names(motrpac_gtex_map)){
#   counts <- table(deg_eqtl_list[[tss]]$gene_id)
#   counts <- counts[counts > 8]
#   print(paste0(tss, " = ", length(counts)))
# }
# 
# 
# tss <- sample(names(motrpac_gtex_map), size = 1)
# par(mfrow = c(4,4))
# for(tss in names(motrpac_gtex_map)){
# inds <- sample(x = length(deg_eqtl_list[[tss]]$logFC), 1E4, replace = F)
# plot(deg_eqtl_list[[tss]]$logFC[inds], deg_eqtl_list[[tss]]$slope[inds], main = tss)
# # plot(deg_eqtl_list[[tss]]$shrunk_logFC[inds], deg_eqtl_list[[tss]]$logFC[inds], main = tss)
# }



# this now includes the intersection of GTEx and MoTrPAC expressed genes annotated in both species
# there are some duplicate gene IDs (both human and rat). figure out how to deal with these
# direction of DEG: logFC
# direction of eQTL: "slope" ("variant_id", "alt")
#save(deg_eqtl_list, file='rdata/deg_eqtl_list_20201109.RData')

calc_gene_intersect = F
if(calc_gene_intersect){
  sapply(names(deg_eqtl_list), function(tissue_1) sapply(names(deg_eqtl_list), function(tissue_2) 
    round(length(intersect(deg_eqtl_list[[tissue_1]]$gene_name, deg_eqtl_list[[tissue_2]]$gene_name)) / 
    length(union(deg_eqtl_list[[tissue_1]]$gene_name, deg_eqtl_list[[tissue_2]]$gene_name)) * 100, 2)))
}

# add in merged gonads
if(!any(names(deg_eqtl_list) == "t1000-gonads")){
  deg_eqtl_list$`t1000-gonads` = rbind(deg_eqtl_list$`t63-testes`, deg_eqtl_list$`t64-ovaries`)
}

#compare old and new rna_dea files
plot(quantile(old_rna_dea$timewise_dea$logFC, 1:99 / 100), quantile(rna_dea$timewise_dea$logFC, 1:99 / 100))



cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"


#### plot basic figure ####

#getting some plotting code up -- just the raw DE x eQTL first
sign_filter <- T #for the motrpac DE
use_abs_slope <- T
sign_filter_alpha <- 0.1
n_genes_to_label_by_eQTL_effectSizes <- 5
n_genes_to_label_by_logFC <- 5
use_overall_fdr_threshold_across_sex_time <- T

# p_val_cols = colorRampPalette(c("blue", "red"))(100)
p_val_cols = viridis::magma(440)[c(1:20*10, 200+1:20*5, 300+1:20*4, 380+1:20*2, 420+1:20)]



plot_basic_figure = T
if(plot_basic_figure){
  
grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_eQTLxDE", ifelse(use_overall_fdr_threshold_across_sex_time, "", "_timesexthresh"),".pdf"), width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
# png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
par(mfrow = c(18, 1), mar = c(4,4,4,8), xpd = T)
for(tissue in names(deg_eqtl_list)){
# for(tissue in names(deg_eqtl_list)[1]){
  
  # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))

  min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
  
  # if(tissue != names(deg_eqtl_list)[1]){break}
  
  plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  
  #legend for dots and dashes
  stponr <- c(1,7) #set_to_put_on_first_row
  legend(x = 3.7165, y = ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 1.26, 0.25775), col = c(1,1,"lightgray"), lwd = c(NA, 1, 2), lty = c(NA, 1, 2), pch = c(19, NA, NA), 
     legend = c("Point Estimates", "±1 SE", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
  
  rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
  rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
  segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
  segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
  shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
             cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
  mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 2, line = 4) #horiz axis label
  mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression (log_{2}-fold change)"), cex = 2, line = -4.25, side = 2) #vert axis label
  
  #note filtration scheme
  text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
      x = 3.37, y = 2.18, pos = 4, cex = 1.75)
  text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
  
  #plot gtex p-value legend
  xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
  rect(
    xl,
    head(seq(yb,yt,(yt-yb)/100),-1),
    xr,
    tail(seq(yb,yt,(yt-yb)/100),-1),
    col=p_val_cols
  )
  text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
  text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
  addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
  # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
  
  
  
  
  for(sex in 1:2){
  # for(sex in 1){
    
    sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
    
    for(timepoint in 1:4){
    # for(timepoint in 1){
      segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
      segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
      shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                 cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
      
      #plot sex symbol
      text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
      
      timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
      
      if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0 & use_overall_fdr_threshold_across_sex_time){
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
      } else {
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
      }
      
      #incompatible gonads message
      if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
        text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
        next()
      }
      
      if(sign_filter){
        delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
      } else{
        delsub <- deg_eqtl_list[[tissue]][intersect(sex_inds, timepoint_inds),]
      }
      
      
      if(!sign_filter){
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
      } else{
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope - deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se), 
                             max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope + deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope_se))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope - deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se), 
                             max(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope + deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC - deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se), 
                           max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC + deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se))
      }
      tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
      
      #add horizontal line to mark 0
      if(timepoint == 1){
        segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                 y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                 y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
      }
        
      #point estimates
      points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
             y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
      
      print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
      
      #vertical standard-error bars
      segments(x0 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               x1 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y0 = ((delsub$logFC + delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
               y1 = ((delsub$logFC - delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      
      #horizontal standard-error bars
      segments(x0 = ((delsub$abs_slope + delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               x1 = ((delsub$abs_slope - delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y0 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
               y1 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      
      #timepoint axes
      if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
        tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
        tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
        axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
        mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
      }
      
      #sex axes
      if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
        # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
        tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
        tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
        axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
        mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
      }
      
      #label genes
      ptl <- order(delsub$abs_slope, decreasing = T)[1:n_genes_to_label_by_eQTL_effectSizes] #points to label
      ptl <- c(ptl, order(delsub$logFC, decreasing = T)[1:n_genes_to_label_by_logFC]) #points to label
      text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
           y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
           labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
           cex = 1.5, pos = 4)
      
      
    }
      
  }
  # dev.off()
}
dev.off()

}

#### connect the dots by sex ####

#now let's connect the dots! first by sex
sign_filter_alpha <- 0.1
use_abs_slope <- T
n_genes_to_label_by_contrastMagnitude <- 5
cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

plot_sexdots_figure = F
if(plot_sexdots_figure){
  
  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE_connectTheDots_sex.pdf", width = 2000 / 72, height = 6000 / 72, family="Arial Unicode MS")
  par(mfrow = c(16, 1), mar = c(4,4,4,8), xpd = T)
  
  for(tissue in setdiff(names(deg_eqtl_list), c("t63-testes", "t64-ovaries"))){
  # for(tissue in names(deg_eqtl_list)[1]){
    print(tissue)
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,4), ylim = c(0,1.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    #legend for arrows
    text(labels = "\u2642", x = 3.75, y = 0.075, pos = 4, cex = 4, col = cols$Sex[1])
    text(labels = "\u2640", x = 3.925, y = 0.075, pos = 4, cex = 4, col = cols$Sex[2])
    arrows(x0 = 3.815, y0 = 0.075, x1 = 3.935, y1 = 0.075, lwd = 2)
    rect(xleft = 3.75, ybottom = 0, xright = 4, ytop = 0.175, lty = 3, lwd = 2)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 1.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 1, xright = 4, ytop = 1.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 3) #week title
    shadowtext(x = 2, y = 1.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #tissue labels
    mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 2, line = 3.5, side = 1) #horiz axis label
    mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression"), cex = 2, line = -1.5, side = 2) #vert axis label
    mtext(text = latex2exp::TeX("(log_{2}-fold change)"), cex = 2, line = -4.5, side = 2) #vert axis label
    
    #plot gtex p-value legend
    xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 1;
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/100),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/100),-1),
      col=p_val_cols
    )
    text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
    text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
    addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
    # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.37, y = 1.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 1.18, pos = 4, cex = 1.75, font = 2)
    
    for(timepoint in 1:4){
      # for(timepoint in 1){
      segments(x0 = timepoint, x1 = timepoint, y0 = 1.15, y1 = 1, lwd = 3) #timepoint divisor title
      segments(x0 = timepoint, x1 = timepoint, y0 = 1, y1 = 0, lwd = 2) #timepoint divisor title
      shadowtext(x = timepoint - 0.5, y = 1, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                 cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
      
      #find indices
      timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
      motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < 0.05)
      motrpac_signif_genes <- unique(deg_eqtl_list[[tissue]]$gene_name[intersect(timepoint_inds, motrpac_signif_inds)])
      motrpac_signif_inds <- unlist(sapply(motrpac_signif_genes, function(signif_gene) which(deg_eqtl_list[[tissue]]$gene_name == signif_gene)))
      male_inds <- which(deg_eqtl_list[[tissue]]$sex == "male")
      female_inds <- which(deg_eqtl_list[[tissue]]$sex == "female")
      sex_inds <- c(male_inds, female_inds)
      if(tissue == "t1000-gonads"){
        motrpac_signif_genes <- intersect(motrpac_signif_genes, 
                                intersect(deg_eqtl_list[[tissue]]$gene_name[male_inds], deg_eqtl_list[[tissue]]$gene_name[female_inds]))
        motrpac_signif_inds <- unlist(sapply(motrpac_signif_genes, function(signif_gene) which(deg_eqtl_list[[tissue]]$gene_name == signif_gene)))
      }
      
      
      
      if(sign_filter){
        delsub_m <- deg_eqtl_list[[tissue]][intersect(intersect(male_inds, timepoint_inds), motrpac_signif_inds),]
        delsub_m <- delsub_m[order(delsub_m$gene_name),]
        delsub_f <- deg_eqtl_list[[tissue]][intersect(intersect(female_inds, timepoint_inds), motrpac_signif_inds),]
        delsub_f <- delsub_f[order(delsub_f$gene_name),]
        if(!all(delsub_f$gene_name == delsub_m$gene_name)){
          print("incompatible gene sets?") #ah due to the duplication
          #let's just deduplicate naively, keeping only the first entry
          tablef <- table(delsub_f$gene_name)
          delsub_f <- delsub_f[-unlist(sapply(names(tablef)[tablef > 1], function(x) which(delsub_f$gene_name == x)[-1])),]
          tablem <- table(delsub_m$gene_name)
          delsub_m <- delsub_m[-unlist(sapply(names(tablem)[tablem > 1], function(x) which(delsub_m$gene_name == x)[-1])),]
          if(!all(delsub_f$gene_name == delsub_m$gene_name)){
            break()
          }
        }
        delsub <- rbind(delsub_m, delsub_f)
      } else{
        # delsub <- deg_eqtl_list[[tissue]][intersect(sex_inds, timepoint_inds),]
      }
      
      
      if(!sign_filter){
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
      } else{
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope), 
                             max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope), 
                             max(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC), 
                           max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC))
      }
      
      tarb <- 0.075 #total_axis_rescaling_buffer, so points don't fall on lines
      
      #add horizontal line to mark 0
      if(timepoint == 1){
        segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                 y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
                 y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb)
      }
      
      #arrows connecting M & F points
      arrows(x0 = (delsub_m$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
             x1 = (delsub_f$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
             y0 = (delsub_m$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             y1 = (delsub_f$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_m$pval_beta) / min_logpval * 100)], alpha.f = 0.75), length = 0.2, lwd = 1.75)
      
      # #point estimates
      # points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #        y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #        pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
      # 
      # #vertical standard-error bars
      # segments(x0 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          x1 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          y0 = ((delsub$logFC + delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
      #          y1 = ((delsub$logFC - delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #          grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      # 
      # #horizontal standard-error bars
      # segments(x0 = ((delsub$abs_slope + delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          x1 = ((delsub$abs_slope - delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          y0 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
      #          y1 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #          grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      
      #timepoint axes
      if(T){
        tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
        tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
        axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.15)
        mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = 0)
      }
      
      #vertical axis
      if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
        # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
        tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
        tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb
        axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
        mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
      }
      
      #label genes
      ptl <- order(abs(delsub_m$logFC - delsub_f$logFC), decreasing = T)[1:n_genes_to_label_by_contrastMagnitude] #points to label
      text(x = (delsub_m$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
           y = (delsub_m$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
           labels = delsub_m$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_m$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
           cex = 1.5, pos = c(1,3)[as.integer(delsub_m$logFC[ptl] - delsub_f$logFC[ptl] > 0) + 1])
      text(x = (delsub_f$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
           y = (delsub_f$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
           labels = delsub_f$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_f$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
           cex = 1.5, pos = c(1,3)[as.integer(delsub_f$logFC[ptl] - delsub_m$logFC[ptl] > 0) + 1])
      
      
    }
    
  
    # dev.off()
  }
  dev.off()
  
}

#### connect the dots by time ####

#now let's connect the dots! now by time
sign_filter_alpha <- 0.1
use_abs_slope <- T
jitter_eQTLs <- F; jitterVariance <- 0.0001
n_genes_to_label_by_contrastMagnitude <- 5
cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

plot_timedots_figure = F
if(plot_timedots_figure){
  
  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE_connectTheDots_time.pdf", width = 2000 / 72 / 1 / 1.75, height = 6000 / 72 * 18 / 16, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,4,4,8), xpd = T)
  
  for(tissue in names(deg_eqtl_list)){
  # for(tissue in names(deg_eqtl_list)[1]){
    print(tissue)
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,2), ylim = c(0,1.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    # legend for arrows
    arrowShift <- 0.025
    text(labels = "1W", x = 1.525 + arrowShift, y = 0.05, pos = 4, cex = 2, col = 1)
    text(labels = "2W", x = 1.65 + arrowShift, y = 0.05, pos = 4, cex = 2, col = 1)
    text(labels = "4W", x = 1.775 + arrowShift, y = 0.05, pos = 4, cex = 2, col = 1)
    text(labels = "8W", x = 1.9 + arrowShift, y = 0.05, pos = 4, cex = 2, col = 1)
    
    arrows(x0 = 1.592 + arrowShift, y0 = 0.055, x1 = 1.655 + arrowShift, y1 = 0.055, lwd = 2, length = 0.1, col = cols$Time[2])
    arrows(x0 = 1.655 + arrowShift, y0 = 0.055, x1 = 1.592 + arrowShift, y1 = 0.055, lwd = 2, length = 0.05, col = cols$Time[2], angle = 90)
    arrows(x0 = 1.592 + 0.125 + arrowShift, y0 = 0.055, x1 = 1.655 + 0.125 + arrowShift, y1 = 0.055, lwd = 2, length = 0.1, col = cols$Time[3])
    arrows(x0 = 1.592 + 0.125*2 + arrowShift, y0 = 0.055, x1 = 1.655 + 0.125*2 + arrowShift, y1 = 0.055, lwd = 2, length = 0.1, col = cols$Time[4])
    rect(xleft = 1.525 + arrowShift, ybottom = 0, xright = 2, ytop = 0.125, lty = 3, lwd = 2)
    
    rect(xleft = 0, ybottom = 0, xright = 2, ytop = 1.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 1, xright = 2, ytop = 1.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 2, y0 = 1, y1 = 1, lwd = 3) #week title
    shadowtext(x = 1, y = 1.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #tissue labels
    mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 1.5, line = 3, side = 1) #horiz axis label
    mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression"), cex = 1.5, line = 0.5, side = 2) #vert axis label
    mtext(text = latex2exp::TeX("(log_{2}-fold change)"), cex = 1.5, line = -1.5, side = 2) #vert axis label
    
    #plot gtex p-value legend
    xl <- 2.025; yb <- 0; xr <- 2.065; yt <- 1;
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/100),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/100),-1),
      col=p_val_cols
    )
    text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 2.06, pos = 4, las=2, cex=1.5)
    text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
    addImg(GTEx_logo, x = xr + 0.009, y = yt + 0.125, width = 0.1) #much better
    # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 1.4175, y = 1.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 1.91, y = 1.18, pos = 4, cex = 1.75, font = 2)
    
    for(sex in 1:2){
      
      # for(timepoint in 1){
      segments(x0 = sex, x1 = sex, y0 = 1.15, y1 = 1, lwd = 3) #timepoint divisor title
      segments(x0 = sex, x1 = sex, y0 = 1, y1 = 0, lwd = 2) #timepoint divisor title
      shadowtext(x = sex - 0.5, y = 1 + ifelse(sex == 1, 0, 0.01), labels = c("\u2642", "\u2640")[sex], 
                 cex = ifelse(sex == 1, 5.5, 4.5), col = cols$Sex[sex], pos = 3) #timepoint labels
      
      if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
        text(x = 0.5 + sex - 1, y = 0.5 , labels = "(not available)", cex = 2)
        next()
      }
      
      #find indices
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < 0.05)
      motrpac_signif_genes <- unique(deg_eqtl_list[[tissue]]$gene_name[intersect(sex_inds, motrpac_signif_inds)])
      motrpac_signif_inds <- unlist(sapply(motrpac_signif_genes, function(signif_gene) which(deg_eqtl_list[[tissue]]$gene_name == signif_gene)))
      timepoint_inds_1 <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[1])
      timepoint_inds_2 <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[2])
      timepoint_inds_3 <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[3])
      timepoint_inds_4 <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[4])
      timepoint_inds <- c(timepoint_inds_1, timepoint_inds_2, timepoint_inds_3, timepoint_inds_4)
      
      if(sign_filter){
        delsub_t1 <- deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds_1, sex_inds), motrpac_signif_inds),]
        delsub_t1 <- delsub_t1[order(delsub_t1$gene_name),]
        delsub_t2 <- deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds_2, sex_inds), motrpac_signif_inds),]
        delsub_t2 <- delsub_t2[order(delsub_t2$gene_name),]
        delsub_t3 <- deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds_3, sex_inds), motrpac_signif_inds),]
        delsub_t3 <- delsub_t3[order(delsub_t3$gene_name),]
        delsub_t4 <- deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds_4, sex_inds), motrpac_signif_inds),]
        delsub_t4 <- delsub_t4[order(delsub_t4$gene_name),]
        
        if(!all(c(delsub_t1$gene_name == delsub_t2$gene_name,
                  delsub_t1$gene_name == delsub_t3$gene_name,
                  delsub_t1$gene_name == delsub_t4$gene_name
                  ))){
          print("incompatible gene sets?") #ah due to the duplication
          #let's just deduplicate naively, keeping only the first entry
          table_t1 <- table(delsub_t1$gene_name)
          delsub_t1 <- delsub_t1[-unlist(sapply(names(table_t1)[table_t1 > 1], function(x) which(delsub_t1$gene_name == x)[-1])),]
          
          table_t2 <- table(delsub_t2$gene_name)
          delsub_t2 <- delsub_t2[-unlist(sapply(names(table_t2)[table_t2 > 1], function(x) which(delsub_t2$gene_name == x)[-1])),]
          
          table_t3 <- table(delsub_t3$gene_name)
          delsub_t3 <- delsub_t3[-unlist(sapply(names(table_t3)[table_t3 > 1], function(x) which(delsub_t3$gene_name == x)[-1])),]
          
          table_t4 <- table(delsub_t4$gene_name)
          delsub_t4 <- delsub_t4[-unlist(sapply(names(table_t4)[table_t4 > 1], function(x) which(delsub_t4$gene_name == x)[-1])),]
          
          if(!all(c(delsub_t1$gene_name == delsub_t2$gene_name,
                    delsub_t1$gene_name == delsub_t3$gene_name,
                    delsub_t1$gene_name == delsub_t4$gene_name
          ))){
            break()
          }
        }
        delsub <- rbind(delsub_m, delsub_f)
      } else{
        # delsub <- deg_eqtl_list[[tissue]][intersect(sex_inds, timepoint_inds),]
      }
      
      if(jitter_eQTLs){
        delsub_t1$abs_slope <- delsub_t1$abs_slope + rnorm(length(delsub_t1$abs_slope), 0, sqrt(jitterVariance))
        delsub_t2$abs_slope <- delsub_t2$abs_slope + rnorm(length(delsub_t2$abs_slope), 0, sqrt(jitterVariance))
        delsub_t3$abs_slope <- delsub_t3$abs_slope + rnorm(length(delsub_t3$abs_slope), 0, sqrt(jitterVariance))
        delsub_t4$abs_slope <- delsub_t4$abs_slope + rnorm(length(delsub_t4$abs_slope), 0, sqrt(jitterVariance))
      }  
      
      
      
      if(!sign_filter){
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
      } else{
        if(!use_abs_slope){
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope), 
                             max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope))
        } else {
          slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope), 
                             max(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope))
        }
        logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC), 
                           max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC))
      }
      tarb <- 0.075 #total_axis_rescaling_buffer, so points don't fall on lines
      
      #add horizontal line to mark 0
      if(timepoint == 1){
        segments(x0 = 0, x1 = 2, lwd = 2, col = "lightgray", lty = 2,
                 y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
                 y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb)
      }
      
      #arrows connecting timepoints
      arrows(x0 = (delsub_t2$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             x1 = (delsub_t1$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             y0 = (delsub_t2$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             y1 = (delsub_t1$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             col = cols$Time[2], length = 0.05, lwd = 2, angle = 90)
      
      arrows(x0 = (delsub_t1$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             x1 = (delsub_t2$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             y0 = (delsub_t1$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             y1 = (delsub_t2$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             col = cols$Time[2], length = 0.1, lwd = 2)
      
      arrows(x0 = (delsub_t2$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             x1 = (delsub_t3$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             y0 = (delsub_t2$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             y1 = (delsub_t3$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             col = cols$Time[3], length = 0.1, lwd = 2)
      
      arrows(x0 = (delsub_t3$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             x1 = (delsub_t4$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
             y0 = (delsub_t3$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             y1 = (delsub_t4$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
             col = cols$Time[4], length = 0.1, lwd = 2)
      
      # #point estimates
      # points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #        y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #        pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
      # 
      # #vertical standard-error bars
      # segments(x0 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          x1 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          y0 = ((delsub$logFC + delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
      #          y1 = ((delsub$logFC - delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #          grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      # 
      # #horizontal standard-error bars
      # segments(x0 = ((delsub$abs_slope + delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          x1 = ((delsub$abs_slope - delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
      #          y0 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
      #          y1 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
      #          grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
      
      #timepoint axes
      if(T){
        tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
        tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb
        axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.15)
        mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = 0)
      }
      
      #vertical axis
      if(sex == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
        # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
        tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
        tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb
        axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -4)
        mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -3)
      }
      
      #label genes
      ptl <- order(abs(delsub_t1$logFC - delsub_t2$logFC) + abs(delsub_t2$logFC - delsub_t3$logFC) + abs(delsub_t3$logFC - delsub_t4$logFC), 
                   decreasing = T)[1:n_genes_to_label_by_contrastMagnitude] #points to label
      max_FC <- apply(cbind(delsub_t1$logFC, delsub_t2$logFC, delsub_t3$logFC, delsub_t4$logFC), 1, max)
      min_FC <- apply(cbind(delsub_t1$logFC, delsub_t2$logFC, delsub_t3$logFC, delsub_t4$logFC), 1, min)
      text(x = (delsub_t1$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
           y = (max_FC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
           labels = delsub_t1$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_t1$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
           cex = 1.5, pos = 3)
      text(x = (delsub_t4$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + sex - 1 + tarb,
           y = (min_FC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + tarb,
           labels = delsub_t4$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_t4$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
           cex = 1.5, pos = 1)
      
      
    }
    
    
    # dev.off()
  }
  dev.off()
  
}


##### draw contour envelopes figure ####

#now let's draw some envelopes

sign_filter_alpha <- 0.1
use_abs_slope <- T
prop_in_envelope <- 0.85
cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

if(!exits("deg_eqtl_list_all")){
  deg_eqtl_list_all <- do.call(rbind, deg_eqtl_list)
  deg_eqtl_list_all$tissue[deg_eqtl_list_all$tissue == "t64-ovaries"] <- "t1000-gonads"
  deg_eqtl_list_all$tissue[deg_eqtl_list_all$tissue == "t63-testes"] <- "t1000-gonads"
  deg_eqtl_list_all$tissue_abbreviation[deg_eqtl_list_all$tissue_abbreviation == "TESTES"] <- "GONADS"
  deg_eqtl_list_all$tissue_abbreviation[deg_eqtl_list_all$tissue_abbreviation == "OVARIES"] <- "GONADS"
}

plot_tissueEnvelope_figure = F
if(plot_tissueEnvelope_figure){
  
  for(g in 1:2){
    
    if(g == 1){
      grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/basic_figure_tissueEnvelope.pdf", width = 2000 / 72, height = 600 / 72, family="Arial Unicode MS")
      # png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
      par(mfrow = c(1, 1), mar = c(4,2,4,8), xpd = T)
    
      min_logpval <- floor(min(log(deg_eqtl_list_all$pval_beta))) #gtex p-val
      
      # if(tissue != names(deg_eqtl_list)[1]){break}
      
      plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    }  
      #legend for dots and dashes
      # stponr <- c(1,7) #set_to_put_on_first_row
      # legend(x = 3.7165, y = ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 1.26, 0.25775), col = c(1,1,"lightgray"), lwd = c(NA, 1, 2), lty = c(NA, 1, 2), pch = c(19, NA, NA), 
      #        legend = c("Point Estimates", "±1 SE", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    if(g == 2){
      rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
      rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
      segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
      segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
      # shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
      #            cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
      mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 2, line = 2, side = 1) #horiz axis label
      mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression (log_{2}-fold change)"), cex = 2, line = -2, side = 2) #vert axis label
      
      #note filtration scheme
      text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
           x = 3.0825, y = 2.19, pos = 4, cex = 1.75)
      text(labels = sign_filter_alpha, x = 3.91, y = 2.19, pos = 4, cex = 1.75, font = 2)
      
      #note polygon envelope boundaries
      text(labels = paste0("tetragons drawn to contain approximately ", prop_in_envelope * 100, "% of the mass"), x = 0, y = 2.19, pos = 4, cex = 1.75)
      
      #draw tissue legend
      plotMatrix(t(t(sapply(names(cols$Tissue), function(name) paste0(strsplit(x = name, split = "-")[[1]][-1], collapse = "-")))), size = c(0.225,1.25), 
                 location = c(4.025,0), title = F, rownames = F, colnames = F, cex = 1.1)
      rect(xleft = 4.26, 
          xright = 4.26 + 1.25 / 36,
          ybottom = 0 + 0:17 * 1.25 / 18,
          ytop = 0 + 1:18 * 1.25 / 18,
          col = rev(grDevices::adjustcolor(cols$Tissue, alpha.f = 0.5)))
      text(x = 4.075, y = 1.25, labels = "Tissue", pos = 3, cex = 1.5, font = 2)
      
      
      #plot gtex p-value legend
      # xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
      # rect(
      #   xl,
      #   head(seq(yb,yt,(yt-yb)/100),-1),
      #   xr,
      #   tail(seq(yb,yt,(yt-yb)/100),-1),
      #   col=p_val_cols
      # )
      # text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
      # text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
      # addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
      # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    }
      for(sex in 1:2){
        
        sex_inds <- which(deg_eqtl_list_all$sex == c("male", "female")[sex])
        
        for(timepoint in 1:4){
          
          cat(paste0(c("M", "F")[sex],timepoint, " "))
          
          if(g == 2){
            segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
            segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
            shadowtext(x = timepoint - 0.5, y = 1.99, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                       cex = 2.5, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
            
            #plot sex symbol
            text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825) - 0.025, pos = 4, cex = 3, col = cols$Sex[sex])
          }
          
          timepoint_inds <- which(deg_eqtl_list_all$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
          
          if(length(deg_eqtl_list_all$selection_fdr) > 0){
            motrpac_signif_inds <- which(deg_eqtl_list_all$selection_fdr < sign_filter_alpha)
          } else {
            motrpac_signif_inds <- which(deg_eqtl_list_all$adj_p_value < sign_filter_alpha)
          }
          
          if(!sign_filter){
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list_all[timepoint_inds,]$slope), max(deg_eqtl_list_all[sex_inds,]$slope))
            } else {
              slope_rescale <- c(min(deg_eqtl_list_all[timepoint_inds,]$abs_slope), max(deg_eqtl_list_all[sex_inds,]$abs_slope))
            }
            logFC_rescale <- c(min(deg_eqtl_list_all[sex_inds,]$logFC), max(deg_eqtl_list_all[timepoint_inds,]$logFC))
          } else{
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list_all[intersect(timepoint_inds, motrpac_signif_inds),]$slope), 
                                 max(deg_eqtl_list_all[intersect(sex_inds, motrpac_signif_inds),]$slope))
            } else {
              slope_rescale <- c(min(deg_eqtl_list_all[intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope), 
                                 max(deg_eqtl_list_all[intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope))
            }
            logFC_rescale <- c(min(deg_eqtl_list_all[intersect(sex_inds, motrpac_signif_inds),]$logFC), 
                               max(deg_eqtl_list_all[intersect(sex_inds, motrpac_signif_inds),]$logFC))
          }
          tarb <- 0.025 #total_axis_rescaling_buffer, so points don't fall on lines
          
          for(tissue in unique(deg_eqtl_list_all$tissue)){
            
            tissue_inds <- which(deg_eqtl_list_all$tissue == tissue)
            
            delsub <- deg_eqtl_list_all[intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), tissue_inds),]
            
            if(g == 2){
              #add horizontal line to mark 0
              if(timepoint == 1){
                segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                         y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                         y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
              }
            }
            
            if(g == 1){
            #compute bivariate kde
            bvkde <- MASS::kde2d(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                                 y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, n = 100)
            # bvkde <- KernSmooth::bkde2D(x = cbind((delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
            #                      (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb), bandwidth = c(0.00001, 0.0002))
            # names(bvkde) <- c("x", "y", "z")
            # landscape <- cbind((delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
            #                    (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
            # bandw <- ks::Hpi(landscape,binned=T)
            # bvkde <- kde(landscape,H=bandw)
            # bvkde <- kde.truncate(bvkde, matrix(c(-Inf, Inf, 0, Inf),2,2, byrow = F))
            
            #find appropriate contour height to contain prop_in_envelope of mass
            mass_grid <- diff(bvkde$y)[1] * diff(bvkde$x)[1] * bvkde$z
            dens_seq <- seq(max(bvkde$z), min(bvkde$z), length.out = 1E3)
            mass_contained = 0
            index <- 1
            while(mass_contained < prop_in_envelope & index < length(dens_seq)){
              highwater <- dens_seq[index]
              mass_contained <- sum(mass_grid[bvkde$z > highwater])
              index = index + 1
            }
            
            # contour(x = bvkde$x, y = bvkde$y, z = bvkde$z, levels = highwater, col = cols$Tissue[names(cols$Tissue) == tissue],
            #         labels = "", add = T, lwd = 2)
            # .filled.contour(bvkde$x, bvkde$y, bvkde$z, levels = highwater, col = cols$Tissue[names(cols$Tissue) == tissue])
            # plot(bvkde,cont=c(prop_in_envelope * 100),lty=2, lwd=0, display="filled.contour2",
            #      col=c(NA,grDevices::adjustcolor(cols$Tissue[names(cols$Tissue) == tissue], alpha.f = 0.5)),
            #      add=T)
            # holy crap how hard is it to get a single contour-filling function around here :[
            # let's try drawing the vertices of a polygon
            ptd <- bvkde$z > highwater
            
            tops <- apply(ptd, 2, minwhich)
            tops <- cbind(which(tops < Inf), tops[tops < Inf])
            
            rights <- apply(ptd, 1, maxwhich)
            rights <- cbind(rights[rights > -Inf], which(rights > -Inf))
            
            bottoms <- apply(ptd, 2, maxwhich)
            bottoms <- cbind(rev(which(bottoms > -Inf)), rev(bottoms[bottoms > -Inf]))
            
            lefts <- apply(ptd, 1, minwhich)
            lefts <- cbind(rev(lefts[lefts < Inf]), rev(which(lefts < Inf)))
            
            vs <- rbind(tops[which.min(tops[,2]),], rights[which.max(rights[,1]),], bottoms[which.max(bottoms[,2]),], lefts[which.min(lefts[,1]),])
            # vs <- rbind(tops, rights, apply(bottoms, 2, rev), lefts, tops[1,])
            polygon(x = bvkde$x[vs[,1]], y = bvkde$y[vs[,2]], col = grDevices::adjustcolor(cols$Tissue[names(cols$Tissue) == tissue], alpha.f = 0.5), 
                    border = grDevices::adjustcolor(cols$Tissue[names(cols$Tissue) == tissue], alpha.f = 0.5))
            
            }
            
            
            if(g == 2){
            
          #timepoint axes
          if(sex == 1 & tissue == unique(deg_eqtl_list_all$tissue)[1]){
            tick_vals <- trunc(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10) * 10) / 10
            tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
            axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -1.25)
            mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1.15, line = -0.5)
          }
          
          #sex axes
          if(timepoint == 1 & tissue == unique(deg_eqtl_list_all$tissue)[1]){
            # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
            tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
            tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
            axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -4.74)
            mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1.15, line = -4.1)
          }
            }
          
        }
        
      }
    
      }
    
    if(g == 2){
    dev.off()
    } else {
      box(which = "plot", lwd = 4, col = "white")
    }
  }
}

#### colocalization figure -- just the colocs ####


if(!file.exists('~/data/smontgom/coloc_list.RData')){
  coloc_list = list()
  pp4_threshold <- 0.01
  i = 1
  for(tissue in unname(motrpac_gtex_map)){
    # download all coloc results for this tissue 
    # only save coloc with PP4 > 0.5
    motrpac_tissues = names(motrpac_gtex_map)[motrpac_gtex_map==tissue]
    for(cfile in list.files(sprintf('%s/%s', coloc, tissue), full.names = T)){
      cat(gsub('__PM__.*','',basename(cfile)), '\n')
      c = fread(cfile, sep='\t', header=T)
      c = c[p4 > pp4_threshold]
      c[,gtex_tissue := tissue]
      
      if(length(motrpac_tissues)>1){
        clist = list()
        for(j in motrpac_tissues){
          cj <- copy(c)
          cj[,motrpac_tissue := j]
          clist[[j]] = cj
        }
        c = rbindlist(clist)
      }else{
        c[,motrpac_tissue := motrpac_tissues]
      }
      c[,gwas_trait := gsub('__PM__.*','',basename(cfile))]
      coloc_list[[i]] = c
      i = i+1
    }
  }
  coloc_list <- as.data.frame(do.call(rbind, coloc_list))
  coloc_list <- coloc_list[,c("gene_id", "p4", "motrpac_tissue", "gwas_trait")]
  save(coloc_list, file=paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
}

load('~/data/smontgom/coloc_list.RData')
colocs <- coloc_list
gonocs <- coloc_list[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
gonocs$motrpac_tissue <- "t1000-gonads"
colocs <- rbind(colocs, gonocs)
# table(coloc_list$motrpac_tissue)


coloc_phenotypes <- sort(unique(colocs$gwas_trait))
n_genes_to_label_by_coloc_PP4 <- 5
plot_coloc_figures = F
if(plot_coloc_figures){
  
  if(!exists("cl")){
    cl <- makeCluster(16, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()
  
  foreach(coloc_phenotype=coloc_phenotypes, .packages = c("latex2exp")) %dopar% {
  # for(coloc_phenotype in coloc_phenotypes){
    print(coloc_phenotype)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/coloc_filter/", coloc_phenotype, "_eQTLxDE.pdf"), width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
    # png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
    par(mfrow = c(18, 1), mar = c(4,4,4,8), xpd = T)
    for(tissue in names(deg_eqtl_list)){
      
      # for(tissue in names(deg_eqtl_list)[1]){
      # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
      
      coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
      coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
      
      
      pp4s <- unname(unlist(sapply(coloc_genes, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue])))
      if(length(pp4s) > 0){
        pp4s[pp4s > (1 - 1E-10)] <- (1 - 1E-10) #*sigh* totes getting those probs of 1 in 10B lol
        min_logpp4 <- floor(min(log10(1-pp4s))) #gtex p-val
      }
      # if(tissue != names(deg_eqtl_list)[1]){break}
      
      plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
      
      #legend for dots and dashes
      stponr <- c(1,7) #set_to_put_on_first_row
      legend(x = 3.7165, y = ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 1.26, 0.25775), col = c(1,1,"lightgray"), lwd = c(NA, 1, 2), lty = c(NA, 1, 2), pch = c(19, NA, NA), 
             legend = c("Point Estimates", "±1 SE", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
      
      rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
      rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
      segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
      segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
      shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
                 cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
      mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 2, line = 4) #horiz axis label
      mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression (log_{2}-fold change)"), cex = 2, line = -4.25, side = 2) #vert axis label
      
      #note filtration scheme
      text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
           x = 3, y = 2.18, pos = 4, cex = 1.75)
      text(labels = paste0(sign_filter_alpha, " & Colocalization PP4 > 0.5"), x = 3 + 0.54, y = 2.18, pos = 4, cex = 1.75, font = 2)
      text(labels = paste0("Top ", n_genes_to_label_by_coloc_PP4, " Colocalization 'PP4's labeled"), x = 0, y = 2.18, pos = 4, cex = 1.75, font = 2)
      
      #plot gtex p-value legend
      xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
      rect(
        xl,
        head(seq(yb,yt,(yt-yb)/100),-1),
        xr,
        tail(seq(yb,yt,(yt-yb)/100),-1),
        col=p_val_cols
      )
      text(labels = c(0:-min_logpp4), y = seq(yb, yt - 0.025, length.out = length(0:-min_logpp4)), x = 4.07, pos = 4, las=2, cex=1.5)
      text(labels = latex2exp::TeX("-log_{10}(1-PP4)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
      # addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
      # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
      
      
      
      
      for(sex in 1:2){
        
        
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        for(timepoint in 1:4){
          
          segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
          segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
          shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                     cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
          
          #plot sex symbol
          text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
          
          timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
          
          if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
          } else {
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
          }
          
          #incompatible gonads message
          if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
            text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
            next()
          }
          
          
          
          if(sign_filter){
            delsub <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
            if(length(delsub$feature_ID) == 0){
              text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "nothing passes our filters :[", cex = 2)
              next()
            }
            delsub$pp4 <- sapply(delsub$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1])
            
          } else{
            delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), coloc_inds),]
          }
          
          
          if(!sign_filter){
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
            } else {
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
            }
            logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
          } else{
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$slope - deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$slope_se), 
                                 max(deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$slope + deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$slope_se))
            } else {
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$abs_slope - deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$slope_se), 
                                 max(deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$abs_slope + deg_eqtl_list[[tissue]][intersect(intersect(timepoint_inds, motrpac_signif_inds), coloc_inds),]$slope_se))
            }
            logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$logFC - deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$logFC_se), 
                               max(deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$logFC + deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, motrpac_signif_inds), coloc_inds),]$logFC_se))
          }
          
          #correct logFC_rescale that do not straddle zero
          if(all(logFC_rescale > 0)){
            logFC_rescale[1] <- -0.1
          } else if(all(logFC_rescale < 0)){
            logFC_rescale[2] <- 0.1
          }
          
          tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
          
          #add horizontal line to mark 0
          if(timepoint == 1){
            segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                     y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                     y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
          }
          
          #point estimates
          points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                 pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log10(1-delsub$pp4) / min_logpp4 * 100)], alpha.f = 0.75), cex = 1.5)
          
          
          
          
          
          # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
          
          #vertical standard-error bars
          segments(x0 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = ((delsub$logFC + delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
                   y1 = ((delsub$logFC - delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   col = grDevices::adjustcolor(p_val_cols[ceiling(log10(1-delsub$pp4) / min_logpp4 * 100)], alpha.f = 0.75), lwd = 1.5)
          
          #horizontal standard-error bars
          segments(x0 = ((delsub$abs_slope + delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = ((delsub$abs_slope - delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
                   y1 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   col = grDevices::adjustcolor(p_val_cols[ceiling(log10(1-delsub$pp4) / min_logpp4 * 100)], alpha.f = 0.75), lwd = 1.5)
          
          #timepoint axes
          if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
            tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
            tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
            axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
            mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
          }
          
          #sex axes
          if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
            # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
            tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
            tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
            axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
            mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
          }
          
          #label genes
          ptl <- order(delsub$pp4, decreasing = T)[1:n_genes_to_label_by_coloc_PP4] #points to label
          text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
               y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
               labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log10(1-delsub$pp4[ptl]) / min_logpp4 * 100)], alpha.f = 0.75), 
               cex = 1.5, pos = 4)
          
          
        }
        
      }
      # dev.off()
    }
    dev.off()
  } 
  system("cd Documents/DExEQTL/; zip -r eQTLxDE_ColocsFilter.zip coloc_filter")
}


p_val_cols = viridis::magma(440)[c(1:20*10, 200+1:20*5, 300+1:20*4, 380+1:20*2, 420+1:20)]

#### basic figure, z-scores ####

plot_basic_figure_zcores = F
if(plot_basic_figure_zcores){
  
  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE_zscoreEdition.pdf", width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
  # png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,4,4,8), xpd = T)
  # for(tissue in names(deg_eqtl_list)){
    for(tissue in names(deg_eqtl_list)){
    
    print(tissue)
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    #legend for dots and dashes
    stponr <- c(1,7) #set_to_put_on_first_row
    legend(x = 3.7165, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
           legend = c("Point Estimates", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change, z-score)"), cex = 1.75, line = 3, side = 1) #horiz axis label
    mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression (log_{2}-fold change, z-score)"), cex = 1.75, line = -4.25, side = 2) #vert axis label
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.37, y = 2.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    #plot gtex p-value legend
    xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/100),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/100),-1),
      col=p_val_cols
    )
    text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
    text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
    addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
    # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    
    
    
    
    for(sex in 1:2){
      # for(sex in 1){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        # for(timepoint in 1){
        segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
        segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
        shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                   cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
        
        #plot sex symbol
        text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
        
        if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
        } else {
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        }
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
          next()
        }
        
        if(sign_filter){
          delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
          delsub$abs_slope <- delsub$abs_slope / delsub$slope_se
          delsub$logFC <- delsub$logFC / delsub$logFC_se
          delsub$logFC_se <- 1
          delsub$slope_se <- 1
        } else{
          delsub <- deg_eqtl_list[[tissue]][intersect(sex_inds, timepoint_inds),]
        }
        
        
        # if(!sign_filter){
        #   if(!use_abs_slope){
        #     slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
        #   } else {
        #     slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
        #   }
        #   logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
        # } else{
        #   if(!use_abs_slope){
        #     slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope - deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se), 
        #                        max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope + deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope_se))
        #   } else {
            slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope / deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se - 1), 
                               max(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope / deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se + 1))
          # }
          logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC / deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se   - 1), 
                             max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC / deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se + 1))
        # }
        tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
        
        #add horizontal line to mark 0
        if(timepoint == 1){
          segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                   y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        }
        
        #point estimates
        points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
        
        # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
        
        #timepoint axes
        if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
          tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
          tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
          axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
          mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        }
        
        #sex axes
        if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
          # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
          tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
          tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
          axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
          mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
        }
        
        #label genes
        ptl <- order(delsub$abs_slope, decreasing = T)[1:n_genes_to_label_by_eQTL_effectSizes] #points to label
        text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
             y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
             labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
             cex = 1.5, pos = 4)
        
        
      }
      
    }
    # dev.off()
  }
  dev.off()
  
}

#### ratio figure ####

divide_one_by_the_other = F
n_genes_to_label_by_ratioMagnitude <- 10
if(divide_one_by_the_other){
  
  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/divide_one_by_the_other.pdf", width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,0,4,8), xpd = T)
  for(tissue in names(deg_eqtl_list)){
    # for(tissue in names(deg_eqtl_list)[1]){
    
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    #legend for dots and dashes
    stponr <- c(99) #set_to_put_on_first_row
    legend(x = 3.5155, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
           legend = c("Point Estimate Ratio (log2-scale)", "Same Magnitude DE as LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    mtext(text = latex2exp::TeX("Differential Expression ÷ Lead eVariant eQTL Effect Size (log_{2}ratio)"), cex = 2, line = 3.5, side = 1) #horiz axis label
    text(labels = "density                genes", cex = 1.5, x = -0.0225, y = 0.5, font = 3, srt = 90) #horiz axis label
    text(labels = "density                genes", cex = 1.5, x = -0.0225, y = 1.5, font = 3, srt = 90) #horiz axis label
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.3775, y = 2.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    #plot gtex p-value legend
    xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/100),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/100),-1),
      col=p_val_cols
    )
    text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
    text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
    addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
    # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    
    for(sex in 1:2){
      # for(sex in 1){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        # for(timepoint in 1){
        segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
        segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
        shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                   cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
        
        #plot sex symbol
        text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
        
        if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
        } else {
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        }
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
          next()
        }
        
        
          delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
          delsub$abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) * sign(delsub$logFC)
          delsub$logFC <- 0
        
          slope_rescale <- c(min((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                    deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) * 
                                   sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)), 
                               max((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                      deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) *
                                     sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)))
          
          logFC_rescale <- c(-1,1)
        
        tarb <- 0.125 #total_axis_rescaling_buffer, so points don't fall on lines
        
        #add horizontal line to names
        if(timepoint == 1){
          segments(x0 = 0, x1 = 4, lwd = 1, col = "lightgray", lty = 1,
                   y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        }
        
        # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
        
        #plot KDE for ratios
        kde1d <- density(x = delsub$abs_slope, bw = 0.2)
        lines(x = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb, lwd = 3, xpd = F,
              y = ((-kde1d$y / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb)
        segments(x0 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 x1 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 y0 = ((0 / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                 y1 = ((-kde1d$y / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                 col = "grey90", lwd = 2, xpd = F)
        
        #add vertical line to mark zero
        segments(y0 = sex - 1, y1 = sex, lwd = 2, col = "lightgray", lty = 2,
                 x0 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 x1 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb)
        
        #point estimates
        points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
        
        #timepoint axes
        if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
          tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
          tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
          axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
          mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        }
        
        # #sex axes
        # if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
        #   # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
        #   tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
        #   tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
        #   axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
        #   mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
        # }
        
        #label genes
        ptl <- order(abs(delsub$abs_slope), decreasing = T)[1:n_genes_to_label_by_ratioMagnitude] #points to label
        text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
             y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
             labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
             cex = 1.5, pos = 4, srt = 65)
        
        
      }
      
    }
    # dev.off()
  }
  dev.off()
  
}


#### coloc figure -- all points with circles ####

p_val_cols = viridis::magma(100, begin = 0, end = 0.8)
# plot(1:length(p_val_cols), 1:length(p_val_cols), col = p_val_cols, pch = 19)
coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)

coloc_phenotypes <- sort(unique(colocs$gwas_trait))
n_genes_to_label_by_coloc_PP4 <- 5
plot_coloc_allPoints_figures = F
if(plot_coloc_allPoints_figures){
  
  if(!exists("cl")){
    cl <- makeCluster(16, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()

  foreach(coloc_phenotype=coloc_phenotypes, .packages = c("latex2exp")) %dopar% {
  # for(coloc_phenotype in coloc_phenotypes[51]){
    print(coloc_phenotype)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/coloc_filter_allPoints/", coloc_phenotype, "_eQTLxDE.pdf"), width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
    # png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
    par(mfrow = c(18, 1), mar = c(4,4,4,12), xpd = T)
    for(tissue in names(deg_eqtl_list)){
    
      min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
        
      # for(tissue in names(deg_eqtl_list)[1]){
      # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
      
      coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
      coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
      
      pp4s <- unname(unlist(sapply(coloc_genes, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue])))
      if(length(pp4s) > 0){
        pp4s[pp4s > (1 - 1E-10)] <- (1 - 1E-10) #*sigh* totes getting those probs of 1 in 10B lol
        min_logpp4 <- floor(min(log10(1-pp4s))) #gtex p-val
      }
      # if(tissue != names(deg_eqtl_list)[1]){break}
      
      plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
      
      #legend for dots and dashes
      stponr <- c(1,7) #set_to_put_on_first_row
      legend(x = 3.7125, y = ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 1.26, 0.25775), col = c(1,1,"lightgray"), lwd = c(NA, 1, 2), lty = c(NA, 1, 2), pch = c(19, NA, NA), 
             legend = c("Point Estimates", "±1 SE", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
      
      rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
      rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
      segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
      segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
      shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
                 cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
      mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change)"), cex = 2, line = 4) #horiz axis label
      mtext(text = latex2exp::TeX("MoTrPAC Diff. Expression (log_{2}-fold change)"), cex = 2, line = -4.25, side = 2) #vert axis label
      
      #note filtration scheme
      text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
           x = 3, y = 2.18, pos = 4, cex = 1.75)
      text(labels = paste0(sign_filter_alpha, " & Colocalization PP4 > 0.5"), x = 3 + 0.54, y = 2.18, pos = 4, cex = 1.75, font = 2)
      text(labels = paste0("Top ", n_genes_to_label_by_coloc_PP4, " Colocalization 'PP4's labeled"), x = 0, y = 2.18, pos = 4, cex = 1.75, font = 2)
      
      #plot colocs pp4 legend
      xl <- 4.15; yb <- 0; xr <- 4.2; yt <- 1;
      rect(
        xl,
        head(seq(yb,yt,(yt-yb)/50),-1),
        xr,
        tail(seq(yb,yt,(yt-yb)/50),-1),
        col=coloc_cols
      )
      text(labels = c(0:-min_logpp4), y = seq(yb, yt - 0.025, length.out = length(0:-min_logpp4)), x = xr - 0.005, pos = 4, las=2, cex=1.5)
      text(labels = latex2exp::TeX("-log_{10}(1-PP_{4})"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
      # addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
      # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
      
      #plot gtex p-value legend
      xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
      rect(
        xl,
        head(seq(yb,yt,(yt-yb)/100),-1),
        xr,
        tail(seq(yb,yt,(yt-yb)/100),-1),
        col=p_val_cols
      )
      text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
      text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
      addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
      
      
      for(sex in 1:2){
        
        
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        for(timepoint in 1:4){
          
          segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
          segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
          shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                     cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
          
          #plot sex symbol
          text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
          
          timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
          
          if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
          } else {
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
          }
          
          #incompatible gonads message
          if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
            text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
            next()
          }
          
          
          
          delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
          delsub_colocs <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
          if(length(delsub_colocs$feature_ID) != 0){
            delsub_colocs$pp4 <- sapply(delsub_colocs$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1]) 
          }
          
          
          if(!sign_filter){
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$slope), max(deg_eqtl_list[[tissue]][sex_inds,]$slope))
            } else {
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][timepoint_inds,]$abs_slope), max(deg_eqtl_list[[tissue]][sex_inds,]$abs_slope))
            }
            logFC_rescale <- c(min(deg_eqtl_list[[tissue]][sex_inds,]$logFC), max(deg_eqtl_list[[tissue]][timepoint_inds,]$logFC))
          } else{
            if(!use_abs_slope){
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope - deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se), 
                                 max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope + deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$slope_se))
            } else {
              slope_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope - deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se), 
                                 max(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope + deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$slope_se))
            }
            logFC_rescale <- c(min(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC - deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se), 
                               max(deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC + deg_eqtl_list[[tissue]][intersect(sex_inds, motrpac_signif_inds),]$logFC_se))
          }
          
          #correct logFC_rescale that do not straddle zero
          if(all(logFC_rescale > 0)){
            logFC_rescale[1] <- -0.1
          } else if(all(logFC_rescale < 0)){
            logFC_rescale[2] <- 0.1
          }
          
          tarb <- 0.075 #total_axis_rescaling_buffer, so points don't fall on lines
          
          #add horizontal line to mark 0
          if(timepoint == 1){
            segments(x0 = 0, x1 = 4, lwd = 2, col = "lightgray", lty = 2,
                     y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                     y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
          }
          
          #point estimates
          points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 y = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                 pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
          
          # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
          
          #vertical standard-error bars
          segments(x0 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = ((delsub$logFC + delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
                   y1 = ((delsub$logFC - delsub$logFC_se) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
          
          #horizontal standard-error bars
          segments(x0 = ((delsub$abs_slope + delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = ((delsub$abs_slope - delsub$slope_se) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
                   y1 = (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), lwd = 1.5)
          
          #timepoint axes
          if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
            tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)
            tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
            axis(1, at = tick_locs, labels = rep("", 10), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
            mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
          }
          
          #sex axes
          if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
            # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
            tick_vals <- trunc(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5))
            tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
            axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
            mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
          }
          
          #label genes by eQTL
          ptl <- order(delsub$abs_slope, decreasing = T)[1:n_genes_to_label_by_eQTL_effectSizes] #points to label
          text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
               y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
               labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
               cex = 1.5, pos = 4)
          
          #label genes by colocalization
          if(length(delsub_colocs$feature_ID) != 0){
            ptl <- order(delsub_colocs$pp4, decreasing = T)[1:n_genes_to_label_by_coloc_PP4] #points to label
            shadowtext(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
                 y = (delsub_colocs$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
                 labels = delsub_colocs$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub_colocs$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
                 cex = 1.5, pos = 4, bg = 1)
          }
          
          #circle colocs
          if(length(delsub_colocs$feature_ID) != 0){
            for(np in 1:5){
              points(x = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                     y = (delsub_colocs$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                     pch = 19, col = grDevices::adjustcolor(coloc_cols[ceiling(log10(1-delsub_colocs$pp4) / min_logpp4 * 50)], alpha.f = 1 / np), cex = 1.6 * np)       
            }
            
          }
          
          #label genes by colocalization
          if(length(delsub_colocs$feature_ID) != 0){
            ptl <- order(delsub_colocs$pp4, decreasing = T)[1:min(n_genes_to_label_by_coloc_PP4,nrow(delsub_colocs))] #points to label
            shadowtext(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb + 0.005,
                 y = (delsub_colocs$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
                 labels = delsub_colocs$gene_name[ptl], col = coloc_cols[ceiling(log10(1-delsub_colocs$pp4[ptl]) / min_logpp4 * 50)], srt = 45, 
                 cex = 1.6, pos = 4, font = 4)
            text(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb + 0.005,
                 y = (delsub_colocs$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
                 labels = delsub_colocs$gene_name[ptl], col = grDevices::adjustcolor("white", 0.5), srt = 45, 
                 cex = 1.6, pos = 4, font = 4)
          }
          
          
        }
        
      }
      # dev.off()
    }
    dev.off()
  } 
  # system("cd Documents/DExEQTL/; zip -r eQTLxDE_ColocsFilter_allPoints.zip coloc_filter_allPoints")
}

#### coloc figure -- percentile plot ####

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)

coloc_phenotypes <- sort(unique(colocs$gwas_trait))
plot_coloc_percentile_figure = F

# percentiles <- array(data = NA, dim = c(length(names(deg_eqtl_list)), 2, 4, length(coloc_phenotypes)), 
#                      dimnames = list(names(deg_eqtl_list), c("male", "female"), paste0(c(1,2,4,8), "w"), coloc_phenotypes)) #dims are tissue, sex, time, phenotype
if(!file.exists("~/data/smontgom/coloc_DEeQTL_percentiles")){
  percentiles <- list()
  for(tissue in names(deg_eqtl_list)){
    motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
    print(tissue)
    for(coloc_phenotype in coloc_phenotypes){
      coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
      coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
      for(sex in 1:2){
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        for(timepoint in 1:4){
          timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
          delsub_noncolocs <- deg_eqtl_list[[tissue]][setdiff(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
          delsub_colocs <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
          if(length(delsub_colocs$feature_ID) != 0){
            delsub_colocs$pp4 <- sapply(delsub_colocs$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1])
            absLogFC_percentile <- sapply(delsub_colocs$logFC, function(x) sum(abs(x) > abs(delsub_noncolocs$logFC)) / length(delsub_noncolocs$logFC))
            absSlope_percentile <- sapply(delsub_colocs$abs_slope, function(x) sum(x > delsub_noncolocs$abs_slope) / length(delsub_noncolocs$abs_slope))
            # percentiles[tissue, c("male", "female")[sex], paste0(c(1,2,4,8), "w")[timepoint], coloc_phenotype] <- list(rbind(absLogFC_percentile, absSlope_percentile))
            percentiles[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]][[coloc_phenotype]] <- cbind(absLogFC = absLogFC_percentile, 
                                                                                                                                absSlope = absSlope_percentile, 
                                                                                                                                pp4 = delsub_colocs$pp4,
                                                                                                                                phenotype = which(coloc_phenotype == coloc_phenotypes))
            # print(c(coloc_phenotype, absLogFC_percentile, absSlope_percentile))
          }
          
        }
      }
    }
  } 
  save(percentiles, file = "~/data/smontgom/coloc_DEeQTL_percentiles")
}

load("~/data/smontgom/coloc_DEeQTL_percentiles")
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
if(plot_coloc_percentile_figure){
  
  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/coloc_percentiles_eQTLxDE.pdf", width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
  # png(filename = "~/Documents/DExEQTL/basic_figure_eQTLxDE.png", width = 2000, height = 10000 * 18 / 17, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,4,4,8), xpd = T)
  for(tissue in names(deg_eqtl_list)){
    
    plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    # #legend for dots and dashes
    # stponr <- c(1,7) #set_to_put_on_first_row
    # legend(x = 3.7165, y = ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 1.26, 0.25775), col = c(1,1,"lightgray"), lwd = c(NA, 1, 2), lty = c(NA, 1, 2), pch = c(19, NA, NA), 
    #        legend = c("Point Estimates", "±1 SE", "DE = 0 LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    mtext(text = latex2exp::TeX("lead eVariant eQTL effect size (absolute value, log_{2}-fold change percentile)"), cex = 2, line = 4) #horiz axis label
    mtext(text = latex2exp::TeX("MoTrPAC Differential Expression"), cex = 1.75, line = -1.5, side = 2) #vert axis label
    mtext(text = latex2exp::TeX("(absolute value, log_{2}-fold change percentile)"), cex = 1.75, line = -4.25, side = 2) #vert axis label
    
    #note size legend
    points(x = rep(4.075, 15), cex = 1:15, col = grDevices::adjustcolor("black", 0.75), y = cumsum(1:15) / 120 * 2, pch = 19)
    text(x = 4.075 + 1:15 / 250, y = cumsum(1:15) / 120 * 2, labels = round(1-10^(-10*(1:15)/50), 3), pos = 4, cex = seq(1, 3, length.out = 15))
    text(x = 4.15, font = 3, y = 2.125, pos = 3, labels = "Coloc. PP4", cex = 3)
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.37, y = 2.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    
    for(sex in 1:2){

      for(timepoint in 1:4){
        segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
        segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
        shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                   cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
        
        #plot sex symbol
        text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
          next()
        }
        
        #hack-y rescaling
        slope_rescale <- c(0,1)
        logFC_rescale <- c(0,1)
        tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
        
        
        #timepoint axes
        if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
          tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 11), 2)
          tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
          axis(1, at = tick_locs, labels = rep("", 11), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
          mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        }
        
        #sex axes
        if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
          # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
          tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 2)
          tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
          axis(2, at = tick_locs, labels = rep("", 5), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.35)
          mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
        }
        
        #retrieve appropriate percentiles
        if(length(percentiles[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]]) == 0){
          next()
        }
        if(length(percentiles[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]])[1] == 1){
          percsub <- as.data.frame(percentiles[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]])
        } else {
          percsub <- as.data.frame(do.call(rbind, percentiles[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]]))
        }
        # point estimates
        points(x = (percsub$absSlope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y = (percsub$absLogFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               pch = 19, col = grDevices::adjustcolor(pheno_cols[percsub$phenotype], 0.9), cex = (log10(1-percsub$pp4) / -10 * 50))
        text(labels = percsub$phenotype, x = (percsub$absSlope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb + (log10(1-percsub$pp4) / -200),
               y = (percsub$absLogFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + (log10(1-percsub$pp4) / -80),
               col = "white", srt = -45, cex = (log10(1-percsub$pp4) / -1))
        
        # grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75)

        
        #label genes
        # ptl <- order(delsub$abs_slope, decreasing = T)[1:n_genes_to_label_by_eQTL_effectSizes] #points to label
        # text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.005,
        #      y = (delsub$logFC[ptl] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.025,
        #      labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75), 
        #      cex = 1.5, pos = 4)
        
        
      }
      
    }
    # dev.off()
  }
  dev.off()

  grDevices::cairo_pdf(filename = "~/Documents/DExEQTL/coloc_percentiles_eQTLxDE_colorLegend.pdf", width = 1000 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
  plot(1,1,xlim = c(0,1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  text(paste0(1:61, ": ", coloc_phenotypes), col = pheno_cols, x = 0.5, y = seq(10, 0, length.out = length(pheno_cols)), pos = 3, cex = 2)
  text(labels = "Phenotype Color Legend", x = 0.5, y = 10.25, pos = 3, cex = 4, xpd = NA)
  dev.off()
  
}


#### DE / eQTL vs. PP4 figure ####

pp4_threshold <- 0.01
if(!exists("coloc_list") & !exists("colocs")){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
}
gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
gonocs$motrpac_tissue <- "t1000-gonads"
colocs <- rbind(colocs, gonocs)
coloc_phenotypes <- sort(unique(colocs$gwas_trait))

# coloc_phenotypes <- sort(unique(colocs$gwas_trait))
# deg_eqtl_list_colocs <- deg_eqtl_list
# deg_eqtl_list_colocs <- lapply(names(deg_eqtl_list_colocs), function(tissue) cbind(as.data.frame(deg_eqtl_list_colocs[[tissue]]), 
#                       as.data.frame(matrix(data = NA, nrow = nrow(deg_eqtl_list_colocs[[tissue]]), ncol = length(coloc_phenotypes), 
#                         dimnames = list(rep("", nrow(deg_eqtl_list_colocs[[tissue]])), coloc_phenotypes)))))
# names(deg_eqtl_list_colocs) <- names(deg_eqtl_list)
# for(tissue in names(deg_eqtl_list)){
#   matching_colocs <- lapply(sort(unique(deg_eqtl_list_colocs[[tissue]]$gene_id)), function(gene) which(colocs$gene_id == gene & colocs$motrpac_tissue == tissue))
#   rbind(colocs$gwas_trait[matching_colocs], colocs$p4[matching_colocs])
# }

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)

plot_DE_EQTL_PP4_figure = F

pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
if(plot_DE_EQTL_PP4_figure){
  
  if(!exists("cl")){
    cl <- makeCluster(8, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()
  foreach(coloc_phenotype=coloc_phenotypes, .packages = c("latex2exp")) %dopar% {
  # for(coloc_phenotype in coloc_phenotypes[2]){
    
  grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/dobtoaPP4/dobtoaPP4_", coloc_phenotype,".pdf"), 
                       width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,0,4,8), xpd = T)
  for(tissue in names(deg_eqtl_list)){

    coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
    coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
    
    
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    #legend for dots and dashes
    stponr <- c(99) #set_to_put_on_first_row
    legend(x = 3.5155, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
           legend = c("Point Estimate Ratio (log2-scale)", "Same Magnitude DE as LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    mtext(text = latex2exp::TeX("Abs. Differential Expression ÷ Abs. Lead eVariant eQTL Effect Size (log_{2}ratio)"), cex = 2, line = 3.5, side = 1) #horiz axis label
    text(labels = "density", cex = 1.5, x = -0.0225, y = 0.11, font = 3, srt = 90) #horiz axis label
    text(labels = "density", cex = 1.5, x = -0.0225, y = 1.11, font = 3, srt = 90) #horiz axis label
    text(labels = latex2exp::TeX("-log_{10}(1-PP_{4})"), cex = 3, x = -0.1, y = 0.6, font = 3, srt = 90) #horiz axis label
    text(labels = latex2exp::TeX("-log_{10}(1-PP_{4})"), cex = 3, x = -0.1, y = 1.6, font = 3, srt = 90) #horiz axis label
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.3775, y = 2.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    #plot gtex p-value legend
    xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/100),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/100),-1),
      col=p_val_cols
    )
    text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
    text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
    addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
    # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
    
    for(sex in 1:2){
      # for(sex in 1){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        # for(timepoint in 1){
        segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
        segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
        shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                   cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
        
        #plot sex symbol
        text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
        
        if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
        } else {
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        }
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
          next()
        }
        
        
        delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
        delsub$abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) #* sign(delsub$logFC)
        delsub$logFC <- 0
        
        delsub_colocs <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
        if(length(delsub_colocs$feature_ID) != 0){
          delsub_colocs$abs_slope <- (abs(delsub_colocs$logFC) - delsub_colocs$abs_slope) #* sign(delsub$logFC)
          delsub_colocs$logFC <- 0
          delsub_colocs$pp4 <- sapply(delsub_colocs$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1]) 
          delsub_colocs$pp4[delsub_colocs$pp4 > (1 - 1E-6)] <- (1 - 1E-6)
        }
        
        # slope_rescale <- c(min((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
        #                           deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) * 
        #                          sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)), 
        #                    max((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
        #                           deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) *
        #                          sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)))
        slope_rescale <- c(min((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                  deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope)), 
                           max((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                  deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope)))
        
        logFC_rescale <- c(-1,6) #that the maximum PP4 be 1-1E-6
        
        tarb <- 0.125 #total_axis_rescaling_buffer, so points don't fall on lines
        
        #add horizontal line to names
        if(timepoint == 1){
          segments(x0 = 0, x1 = 4, lwd = 1, col = "black", lty = 1,
                   y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        }
        
        # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
        
        #plot KDE for ratios
        kde1d <- density(x = delsub$abs_slope, bw = 0.2)
        lines(x = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb, lwd = 3, xpd = F,
              y = ((-kde1d$y / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb)
        segments(x0 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 x1 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 y0 = ((0 / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                 y1 = ((-kde1d$y / max(kde1d$y)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                 col = "grey90", lwd = 2, xpd = F)
        
        #add vertical line to mark zero
        segments(y0 = sex - 1, y1 = sex, lwd = 2, col = "lightgray", lty = 2,
                 x0 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 x1 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb)
        
        #point estimates
        lessY <- 0.025
        points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
               y = -lessY + (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
               pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
        
        #timepoint axes
        if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
          if(any(slope_rescale > 0) & any(slope_rescale < 0)){
            tick_vals <- round(seq(from = 0, to = max(abs(slope_rescale)) * sign(slope_rescale[which.max(abs(slope_rescale))]), length.out = 7), 2)  
            tick_vals <- c(rev(seq(from = 0, to = min(abs(slope_rescale)) * sign(slope_rescale[which.min(abs(slope_rescale))]), 
                               by = -diff(tick_vals)[1])), tick_vals)
            tick_vals[-which(diff(tick_vals) < 1E-6)]
          } else {
            tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)  
          }
          tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
          axis(1, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
          mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        }
        
        # vertical lines for colocalizing genes
        if(length(delsub_colocs$feature_ID) != 0){
          segments(x0 = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   y1 = (-(log10(1-delsub_colocs$pp4)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   lwd = 3)
        } else {
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = paste0("(no colocalizations w/ PP4 > ", pp4_threshold,")"), cex = 2)
        }
        
        # points(x = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
        #          y = (-(log10(1-delsub_colocs$pp4)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        
        #sex axes
        if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
          # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
          tick_vals <- (seq(from = 0, to = logFC_rescale[2], length.out = 7))
          tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
          axis(2, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.5)
          mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6)
        }
        
        ### label genes
        if(length(delsub_colocs$feature_ID) != 0){
          n_genes_to_label_by_PP4 <- 3
          ptl <- order((delsub_colocs$pp4), decreasing = T)[1:n_genes_to_label_by_PP4] #points to label
          text(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.01,
               y = (-(log10(1-delsub_colocs$pp4[ptl])) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.0125,
               labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75),
               cex = 1.5, pos = 4, srt = 65)
        }
        
        
      }
      
    }
    # dev.off()
  }
  dev.off()
  }
}

#### DE / eQTL vs. PP4 sum-figure ####

pp4_threshold <- 0.1
if(!exists("coloc_list") & !exists("colocs")){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
}
gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
gonocs$motrpac_tissue <- "t1000-gonads"
colocs <- rbind(colocs, gonocs)
coloc_phenotypes <- sort(unique(colocs$gwas_trait))

# coloc_phenotypes <- sort(unique(colocs$gwas_trait))
# deg_eqtl_list_colocs <- deg_eqtl_list
# deg_eqtl_list_colocs <- lapply(names(deg_eqtl_list_colocs), function(tissue) cbind(as.data.frame(deg_eqtl_list_colocs[[tissue]]), 
#                       as.data.frame(matrix(data = NA, nrow = nrow(deg_eqtl_list_colocs[[tissue]]), ncol = length(coloc_phenotypes), 
#                         dimnames = list(rep("", nrow(deg_eqtl_list_colocs[[tissue]])), coloc_phenotypes)))))
# names(deg_eqtl_list_colocs) <- names(deg_eqtl_list)
# for(tissue in names(deg_eqtl_list)){
#   matching_colocs <- lapply(sort(unique(deg_eqtl_list_colocs[[tissue]]$gene_id)), function(gene) which(colocs$gene_id == gene & colocs$motrpac_tissue == tissue))
#   rbind(colocs$gwas_trait[matching_colocs], colocs$p4[matching_colocs])
# }

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)

plot_DE_EQTL_PP4_figure_sum = T
color_by_category <- T
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))
catnames <- sapply(categories, function(x) strsplit(x, split = c("-"))[[1]][1])
catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
catnames["Anthropometric"] <- "Anthro"
catnames["Cardiometabolic"] <- "Cardio"

coloc_discr_cols <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(9, "Set1")[c(1,2,8,9)])
swap_elements <- function(x,e1,e2){nx <- names(x); xn <- x; xn[e1] <- x[e2]; xn[e2] <- x[e1]; names(xn) <- nx; return(xn)}
names(coloc_discr_cols) <- categories
coloc_discr_cols <- swap_elements(coloc_discr_cols, "Hair morphology", "Blood")
coloc_discr_cols["Hair morphology"] <- "black"
coloc_discr_cols["Allergy"] <- "#CEFA05"
n_genes_to_label_by_pp4 <- 3

pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
if(plot_DE_EQTL_PP4_figure_sum){
  
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/dobtoaPP4", ifelse(color_by_category, "_byCategory"),"_sum.pdf"), 
                         width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
    par(mfrow = c(18, 1), mar = c(4,0,4,8), xpd = T)
    for(tissue in names(deg_eqtl_list)[14]){
      
      # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
      
      min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
      
      # if(tissue != names(deg_eqtl_list)[1]){break}
      
      plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
      
      #color legend for categories
      if(color_by_category){
        legend(x = 3.5375, y = 2, legend = catnames, lwd = 5, col = coloc_discr_cols[names(catnames)], cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 2)
      }
      
      #legend for dots and dashes
      stponr <- c(99) #set_to_put_on_first_row
      legend(x = 3.5155, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
             legend = c("Point Estimate Ratio (log2-scale)", "Same Magnitude DE as LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
      
      rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
      rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
      segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
      segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
      shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
                 cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
      mtext(text = latex2exp::TeX("Abs. Differential Expression ÷ Abs. Lead eVariant eQTL Effect Size (log_{2}ratio)"), cex = 2, line = 3.5, side = 1) #horiz axis label
      text(labels = "density", cex = 1.5, x = -0.0225, y = 0.085, font = 3, srt = 90) #horiz axis label
      text(labels = "density", cex = 1.5, x = -0.0225, y = 1.085, font = 3, srt = 90) #horiz axis label
      text(labels = latex2exp::TeX("Pr(n_{colocs} ≥ 1)"), cex = 2.5, x = -0.1, y = 0.575, font = 3, srt = 90) #horiz axis label
      text(labels = latex2exp::TeX("Pr(n_{colocs} ≥ 1)"), cex = 2.5, x = -0.1, y = 1.575, font = 3, srt = 90) #horiz axis label
      # text(labels = latex2exp::TeX(paste0("-log_{10}\\prod_{}^{",length(pheno_cols),"}(1-PP4)")), cex = 2.5, x = -0.1, y = 1.6, font = 3, srt = 90) #horiz axis label
      
      #note filtration scheme
      text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
           x = 3.3775, y = 2.18, pos = 4, cex = 1.75)
      text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
      
      #plot gtex p-value legend
      xl <- 4.025; yb <- 0; xr <- 4.075; yt <- 2;
      rect(
        xl,
        head(seq(yb,yt,(yt-yb)/100),-1),
        xr,
        tail(seq(yb,yt,(yt-yb)/100),-1),
        col=p_val_cols
      )
      text(labels = round(seq(from = 0, to = -min_logpval, length.out = 10)), y = seq(yb, yt - 0.025, length.out = 10), x = 4.07, pos = 4, las=2, cex=1.5)
      text(labels = latex2exp::TeX("-log_e(p-value)"), pos = 4, x = xl - 0.0175, y = yt + 0.035, cex = 2, font = 2)
      addImg(GTEx_logo, x = xr - 0.01, y = yt + 0.125, width = 0.1) #much better
      # text(labels = "GTEx", pos = 4, x = xl - 0.0175, y = yt + 0.125, cex = 2, font = 4)
      
      for(sex in 1:2){
        # for(sex in 1){
        
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        for(timepoint in 1:4){
          
          
          # for(timepoint in 1){
          segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
          segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
          shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                     cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
          
          #plot sex symbol
          text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
          
          timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
          
          if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha)
          } else {
            motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
          }
          
          #incompatible gonads message
          if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
            text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
            next()
          }
          
          
          delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
          delsub$abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) #* sign(delsub$logFC)
          delsub$logFC <- 0
          
          # slope_rescale <- c(min((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
          #                           deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) * 
          #                          sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)), 
          #                    max((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
          #                           deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope) *
          #                          sign(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC)))
          slope_rescale <- c(min((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                    deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope)), 
                             max((abs(deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$logFC) - 
                                    deg_eqtl_list[[tissue]][intersect(timepoint_inds, motrpac_signif_inds),]$abs_slope)))
          
          kdey_rescale <- c(3,12,2,4,12,4,10,4,12,12,2,2,1,10,9,2,12,10) 
          logFC_rescale <- c(-kdey_rescale[which(tissue == names(deg_eqtl_list))] / 6, 
                             kdey_rescale[which(tissue == names(deg_eqtl_list))]) 
          
          
          tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
          
          #add horizontal line to names
          if(timepoint == 1){
            segments(x0 = 0, x1 = 4, lwd = 1, col = "black", lty = 1,
                     y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                     y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
          }
          
          # print(paste0(tissue, ", ", c("male", "female")[sex], ", timepoint ", timepoint, ", min = ", round(min(delsub$logFC), 3), ", max = ", round(max(delsub$logFC), 3)))
          
          #plot KDE for ratios
          kde1d <- density(x = delsub$abs_slope, bw = 0.2)
          kde1d$y <- kde1d$y[kde1d$x > min(delsub$abs_slope) & kde1d$x < max(delsub$abs_slope)]
          kde1d$x <- kde1d$x[kde1d$x > min(delsub$abs_slope) & kde1d$x < max(delsub$abs_slope)]
          
          lines(x = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb, lwd = 3, xpd = F,
                y = ((-kde1d$y / max(kde1d$y) * 0.8*kdey_rescale[which(tissue == names(deg_eqtl_list))] / 6) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb)
          segments(x0 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = (kde1d$x / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   y0 = ((0 / max(kde1d$y) * -logFC_rescale[1]) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                   y1 = ((-kde1d$y / max(kde1d$y)  * -0.8*logFC_rescale[1]) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb,
                   col = "grey90", lwd = 2, xpd = F)
          
          #add vertical line to mark zero
          segments(y0 = sex - 1, y1 = sex, lwd = 2, col = "lightgray", lty = 2,
                   x0 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                   x1 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb)
          
          #prop genes on either side of line
          prop_right <- round(sum(delsub$abs_slope > 0) / length(delsub$abs_slope) * 100)
          prop_left <- 100 - prop_right
          text(x = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb, 
               y = ((-max(kde1d$y) / 15 * kdey_rescale[which(tissue == names(deg_eqtl_list))] / 6) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb, 
               labels = prop_right, col = "grey35", srt = 270, pos = 4)
          text(x = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb, 
               y = ((-max(kde1d$y) / 15 * kdey_rescale[which(tissue == names(deg_eqtl_list))] / 6) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1.05 + tarb, 
               labels = prop_left, col = "grey35", srt = 90, pos = 2)
          
          #timepoint axes
          if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
            if(any(slope_rescale > 0) & any(slope_rescale < 0)){
              tick_vals <- round(seq(from = 0, to = max(abs(slope_rescale)) * sign(slope_rescale[which.max(abs(slope_rescale))]), length.out = 7), 2)  
              tick_vals <- c(rev(seq(from = 0, to = min(abs(slope_rescale)) * sign(slope_rescale[which.min(abs(slope_rescale))]), 
                                     by = -diff(tick_vals)[1])), tick_vals)
              tick_vals[-which(diff(tick_vals) < 1E-6)]
            } else {
              tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)  
            }
            tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
            axis(1, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
            mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
          }
          
          # vertical lines for colocalizing genes
          curr_heights <-  rep(0, nrow(delsub))
          print(tissue)
          
          if(color_by_category){
            coloc_phenotypes_to_color <- coloc_phenotypes[order(traitwise_partitions[match(coloc_phenotypes,traitwise_partitions[,1]),2])]
          } else {
            coloc_phenotypes_to_color <- coloc_phenotypes
          }
          
          for(coloc_phenotype in coloc_phenotypes_to_color){
            cat(paste0(timepoint, " "))
            coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
            coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
            delsub_colocs <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
            if(length(delsub_colocs$feature_ID) != 0){
              delsub_colocs$abs_slope <- (abs(delsub_colocs$logFC) - delsub_colocs$abs_slope) #* sign(delsub$logFC)
              delsub_colocs$logFC <- 0
              delsub_colocs$pp4 <- sapply(delsub_colocs$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1]) 
              delsub_colocs$pp4[delsub_colocs$pp4 > (1 - 1E-6)] <- (1 - 1E-6)
              colocalizing_genes <- sapply(delsub_colocs$gene_id, function(gid) which(delsub$gene_id == gid)[1])
              
              if(color_by_category){
                color_to_use <- coloc_discr_cols[traitwise_partitions[match(coloc_phenotype,traitwise_partitions[,1]),2]]
              } else {
                color_to_use <- pheno_cols[which(coloc_phenotype == coloc_phenotypes)]
              }
              
              segments(x0 = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                       x1 = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                       y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + curr_heights[colocalizing_genes],
                       y1 = (-(log10(1-delsub_colocs$pp4)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + curr_heights[colocalizing_genes],
                       lwd = 1, col = color_to_use)
              curr_heights[colocalizing_genes] <- curr_heights[colocalizing_genes] + (-(log10(1-delsub_colocs$pp4)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) - 
                (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) 
            }
          }
          
          #point estimates
          lessY <- 0.025
          points(x = (delsub$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
                 y = -lessY + (delsub$logFC / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                 pch = 19, col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta) / min_logpval * 100)], alpha.f = 0.75), cex = 1.5)
          
          
          # points(x = (delsub_colocs$abs_slope / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
          #          y = (-(log10(1-delsub_colocs$pp4)) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
          
          #sex axes
          if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
            # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
            tick_vals <- round(seq(from = 0, to = logFC_rescale[2], by = 10^-floor(log10(logFC_rescale[2])-1)))
            if(length(tick_vals) > 4){tick_vals <- round(seq(from = 0, to = logFC_rescale[2], by = 10^-floor(log10(logFC_rescale[2])-1)*3))}
            tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
            axis(2, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.02, line = -7.5)
            mtext(text = c(0, sapply(tick_vals[-1], function(num) as.expression(bquote(1-10^-.(num))))), side = 2, at = tick_locs, cex = 1, line = -6)
            
            
          }
          
          ### label genes
          if(any(curr_heights > 0)){
            ptl <- order(curr_heights, decreasing = T)[1:n_genes_to_label_by_pp4] #points to label
            
            text(x = (delsub$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb + 0.02,
                 y = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + curr_heights[ptl] + 0.025, 
                 labels = delsub$gene_name[ptl], col = 1, cex = 1, pos = 3, srt = 45)
          }
          
          
        }
        
      }
      # dev.off()
    }
    dev.off()
  }



#### summarizing by trait figure ####

pp4_threshold <- 0.1
# if(!exists("coloc_list") & !exists("colocs")){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
# }
gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
gonocs$motrpac_tissue <- "t1000-gonads"
colocs <- rbind(colocs, gonocs)
coloc_phenotypes <- sort(unique(colocs$gwas_trait))


if(!file.exists(paste0("~/Documents/DExEQTL/probs_noColocs_byTrait_", pp4_threshold, "_threshold"))){
  probs_noColocs <- list()
  for(tissue in names(deg_eqtl_list)){
    
    for(sex in 1:2){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
        
        if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        } else {
          motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        }
        
        delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
        delsub$abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) #* sign(delsub$logFC)
        delsub$logFC <- 0
        
        # vertical lines for colocalizing genes
        prob_noColocs <-  rep(0, length(coloc_phenotypes))
        names(prob_noColocs) <- coloc_phenotypes
        print(tissue)
        for(coloc_phenotype in coloc_phenotypes){
          cat(paste0(timepoint, " "))
          coloc_genes <- intersect(deg_eqtl_list[[tissue]]$gene_id, colocs$gene_id[colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype])
          coloc_inds <- as.vector(unname(unlist(sapply(coloc_genes, function(cg) which(deg_eqtl_list[[tissue]]$gene_id == cg)))))
          delsub_colocs <- deg_eqtl_list[[tissue]][intersect(intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds), coloc_inds),]
          if(length(delsub_colocs$feature_ID) != 0){
            delsub_colocs$abs_slope <- (abs(delsub_colocs$logFC) - delsub_colocs$abs_slope) #* sign(delsub$logFC)
            delsub_colocs$pp4 <- sapply(delsub_colocs$gene_id, function(cg) colocs$p4[colocs$gene_id == cg & colocs$motrpac_tissue == tissue & colocs$gwas_trait == coloc_phenotype][1]) 
            delsub_colocs$pp4[delsub_colocs$pp4 > (1 - 1E-6)] <- (1 - 1E-6)
            prob_noColocs[coloc_phenotype] <- -sum(log10(1-delsub_colocs$pp4))
            
          }
        }
        probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]] <- prob_noColocs 
        
      }
      
    }
    # dev.off()
  }
  save(probs_noColocs, file = paste0("~/Documents/DExEQTL/probs_noColocs_byTrait_", pp4_threshold, "_threshold"))
}

load(file = paste0("~/Documents/DExEQTL/probs_noColocs_byTrait_", pp4_threshold, "_threshold"))
# coloc_phenotypes <- sort(unique(colocs$gwas_trait))
# deg_eqtl_list_colocs <- deg_eqtl_list
# deg_eqtl_list_colocs <- lapply(names(deg_eqtl_list_colocs), function(tissue) cbind(as.data.frame(deg_eqtl_list_colocs[[tissue]]), 
#                       as.data.frame(matrix(data = NA, nrow = nrow(deg_eqtl_list_colocs[[tissue]]), ncol = length(coloc_phenotypes), 
#                         dimnames = list(rep("", nrow(deg_eqtl_list_colocs[[tissue]])), coloc_phenotypes)))))
# names(deg_eqtl_list_colocs) <- names(deg_eqtl_list)
# for(tissue in names(deg_eqtl_list)){
#   matching_colocs <- lapply(sort(unique(deg_eqtl_list_colocs[[tissue]]$gene_id)), function(gene) which(colocs$gene_id == gene & colocs$motrpac_tissue == tissue))
#   rbind(colocs$gwas_trait[matching_colocs], colocs$p4[matching_colocs])
# }

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)
padded_nums <- as.character(1:length(coloc_phenotypes))
padded_nums <- sapply(padded_nums, function(pn) paste0(rep("0", max(nchar(padded_nums)) - nchar(pn)), pn))

plot_phenotype_summary_figure = F
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
if(plot_phenotype_summary_figure){
  
  
  grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/byPhenotype_colocalization_summary_", length(coloc_phenotypes), "traits.pdf"), 
                       width = 2000 / 72, height = 10000 / 72 *18 / 17, family="Arial Unicode MS")
  par(mfrow = c(18, 1), mar = c(4,0,4,0), xpd = T)
  for(tissue in names(deg_eqtl_list)){
    
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    min_logpval <- floor(min(log(deg_eqtl_list[[tissue]]$pval_beta))) #gtex p-val
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    plot(1,1,xlim = c(0,4), ylim = c(0,2.15), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    #legend for dots and dashes
    stponr <- c(99) #set_to_put_on_first_row
    legend(x = 3.5155, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
           legend = c("Point Estimate Ratio (log2-scale)", "Same Magnitude DE as LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    
    rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
               cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    mtext(text = paste0("Phenotype Indices for the ", length(coloc_phenotypes), " GWAS Phenotypes w/ Colocalization Probabilities > ", pp4_threshold), cex = 2, line = 0, side = 1) #horiz axis label
    text(labels = "trait", cex = 1.5, x = -0.015, y = 0.05, font = 3, srt = 90) #horiz axis label
    text(labels = "trait", cex = 1.5, x = -0.015, y = 1.05, font = 3, srt = 90) #horiz axis label
    text(labels = latex2exp::TeX("-log_{10}Pr(No Colocs)"), cex = 2.5, x = -0.075, y = 0.565, font = 3, srt = 90) #horiz axis label
    text(labels = latex2exp::TeX("-log_{10}Pr(No Colocs)"), cex = 2.5, x = -0.075, y = 1.565, font = 3, srt = 90) #horiz axis label
    # text(labels = latex2exp::TeX(paste0("-log_{10}\\prod_{}^{",length(pheno_cols),"}(1-PP4)")), cex = 2.5, x = -0.1, y = 1.6, font = 3, srt = 90) #horiz axis label
    
    #note filtration scheme
    text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
         x = 3.3775, y = 2.18, pos = 4, cex = 1.75)
    text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    for(sex in 1:2){
      # for(sex in 1){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        # for(timepoint in 1){
        segments(x0 = timepoint, x1 = timepoint, y0 = 2.15, y1 = 2, lwd = 3) #timepoint divisor title
        segments(x0 = timepoint, x1 = timepoint, y0 = 2, y1 = 0, lwd = 2) #timepoint divisor title
        shadowtext(x = timepoint - 0.5, y = 2, labels = paste0(c(1,2,4,8), "w")[timepoint], 
                   cex = 4, col = cols$Time[timepoint], pos = 3, r = ifelse(timepoint == 4, 0.05, 0.2)) #timepoint labels
        
        #plot sex symbol
        text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          text(x = 0.5 + timepoint - 1, y = 0.5 + sex - 1, labels = "(not available)", cex = 2)
          next()
        }
        
       
        
        
        #plot KDE for ratios
        
        #add vertical line to mark zero
        # segments(y0 = sex - 1, y1 = sex, lwd = 2, col = "lightgray", lty = 2,
        #          x0 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
        #          x1 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb)
         
        # #timepoint axes
        # if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
        #   if(any(slope_rescale > 0) & any(slope_rescale < 0)){
        #     tick_vals <- round(seq(from = 0, to = max(abs(slope_rescale)) * sign(slope_rescale[which.max(abs(slope_rescale))]), length.out = 7), 2)  
        #     tick_vals <- c(rev(seq(from = 0, to = min(abs(slope_rescale)) * sign(slope_rescale[which.min(abs(slope_rescale))]), 
        #                            by = -diff(tick_vals)[1])), tick_vals)
        #     tick_vals[-which(diff(tick_vals) < 1E-6)]
        #   } else {
        #     tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)  
        #   }
        #   tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
        #   axis(1, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
        #   mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        # }
        
        # vertical lines for colocalizing traits
        trait_locs <- (1:length(coloc_phenotypes) / length(coloc_phenotypes) / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
        trait_probs <- probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]]
        trait_inds <- padded_nums[order(trait_probs, decreasing = F)]
        trait_probs <- sort(trait_probs, decreasing = F)
        # if(tissue == names(deg_eqtl_list)[5] & sex == 1 & timepoint == 4){print(max(trait_probs)); stop("enough of that");}
        
        max_all_trait_probs_for_row <- max(unlist(lapply(1:4, function(tps) probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[tps]]])))
        logFC_rescale <- c(-0.075 * max(max_all_trait_probs_for_row), 1 * max(max_all_trait_probs_for_row))
        if(all(abs(logFC_rescale) < 1E-6)){
          logFC_rescale <- c(-0.075,1)
        }
        slope_rescale <- c(0,1)
        tarb <- 0.05 #total_axis_rescaling_buffer, so points don't fall on lines
      
        segments(x0 = trait_locs, x1 = trait_locs, 
                 y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb, 
                 y1 = (trait_probs / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        
        #label indices of traits
        text(x = trait_locs, labels = trait_inds, pos = 3, col = rep(c("black", "black", "darkred", "darkred"), length(coloc_phenotypes) / 4),
             y = (logFC_rescale[1] / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 
               rep(c(-0.005,-0.0475), length(coloc_phenotypes) / 2))
        
        #add horizontal line to names
        if(timepoint == 1){
          segments(x0 = 0, x1 = 4, lwd = 1, col = "black", lty = 1,
                   y0 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb,
                   y1 = (0 / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb)
        }
        
        #sex axes
        if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
          # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
          tick_vals <- round(seq(from = 0, to = logFC_rescale[2], length.out = 4), max(-floor(log10(logFC_rescale[2])), 1))
          tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
          axis(2, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.75)
          mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6.75)
        }
        
        ### label genes
        # if(length(delsub_colocs$feature_ID) != 0){
        #   n_genes_to_label_by_PP4 <- 3
        #   ptl <- order((delsub_colocs$pp4), decreasing = T)[1:n_genes_to_label_by_PP4] #points to label
        #   text(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.01,
        #        y = (-(log10(1-delsub_colocs$pp4[ptl])) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.0125,
        #        labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75),
        #        cex = 1.5, pos = 4, srt = 65)
        # }
        
        
      }
      
    }
    # dev.off()
  }
  dev.off()
  
  grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/coloc_percentiles_eQTLxDE_colorLegend_", length(coloc_phenotypes), "traits.pdf"), width = 1000 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
  plot(1,1,xlim = c(0,1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  text(paste0(1:length(coloc_phenotypes), ": ", coloc_phenotypes), col = pheno_cols, x = 0.5, y = seq(10, 0, length.out = length(pheno_cols)), pos = 3, cex = 2)
  text(labels = "Phenotype Color Legend", x = 0.5, y = 10.25, pos = 3, cex = 4, xpd = NA)
  dev.off()
}

#### summarizing by trait figure -- nicole edition ####

change_names <- F
change_names_in_plot = T
pp4_threshold <- 0.01
if(!exists("coloc_list") & !exists("colocs") | last_import_pp4t != pp4_threshold | change_names != names_changed){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
  gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
  gonocs$motrpac_tissue <- "t1000-gonads"
  colocs <- rbind(colocs, gonocs)
  coloc_phenotypes <- sort(unique(colocs$gwas_trait))
  last_import_pp4t <- pp4_threshold
  trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
  if(change_names){
    traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
    colocs$gwas_trait <- traitname_map[match(colocs$gwas_trait, traitname_map[,1]),2]
    coloc_phenotypes <- sort(unique(colocs$gwas_trait))
    names_changed = T
  } else {
    names_changed = F
  }
}

load(file = paste0("~/Documents/DExEQTL/probs_noColocs_byTrait_", pp4_threshold, "_threshold"))
# coloc_phenotypes <- sort(unique(colocs$gwas_trait))
# deg_eqtl_list_colocs <- deg_eqtl_list
# deg_eqtl_list_colocs <- lapply(names(deg_eqtl_list_colocs), function(tissue) cbind(as.data.frame(deg_eqtl_list_colocs[[tissue]]), 
#                       as.data.frame(matrix(data = NA, nrow = nrow(deg_eqtl_list_colocs[[tissue]]), ncol = length(coloc_phenotypes), 
#                         dimnames = list(rep("", nrow(deg_eqtl_list_colocs[[tissue]])), coloc_phenotypes)))))
# names(deg_eqtl_list_colocs) <- names(deg_eqtl_list)
# for(tissue in names(deg_eqtl_list)){
#   matching_colocs <- lapply(sort(unique(deg_eqtl_list_colocs[[tissue]]$gene_id)), function(gene) which(colocs$gene_id == gene & colocs$motrpac_tissue == tissue))
#   rbind(colocs$gwas_trait[matching_colocs], colocs$p4[matching_colocs])
# }

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)
padded_nums <- as.character(1:length(coloc_phenotypes))
padded_nums <- (sapply(padded_nums, function(pn) paste0(paste0(rep("0", max(nchar(padded_nums)) - nchar(pn)), collapse = ""), pn)))

plot_phenotype_summary_figure_nicole = T
nicole_mods <- T
arnold_mods <- F
partition_by_category <- T
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))
arnold_bayes <- png::readPNG(source = "~/Pictures/ArnoldBayes.png")
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
reorder_vertical <- T

if(change_names){
  for(i1 in seq_along(probs_noColocs)){
    for(i2 in seq_along(probs_noColocs[[1]])){
      for(i3 in seq_along(probs_noColocs[[1]][[1]])){
        names(probs_noColocs[[i1]][[i2]][[i3]]) <- traitname_map[match(names(probs_noColocs[[i1]][[i2]][[i3]]), traitname_map[,1]),2]
      }
    }
  }
}


if(plot_phenotype_summary_figure_nicole){
  
  for(timepoint in 1:4){
  
  if(reorder_vertical){
    order_traits <- order(apply(sapply(coloc_phenotypes, (function(coloc_phenotype) sapply(1:2, function(sex) sapply(1:4, function(timepoint) sapply(names(deg_eqtl_list), function(tissue) 
        probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]][[coloc_phenotype]]
      ))))), 2, sum), decreasing = T)
  } else {
    order_traits <- 1:length(coloc_phenotypes)
  }
    
  if(partition_by_category){
    order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])], sort(order_traits))]
    trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])] ))
    trait_category_counts$cN <- cumsum(trait_category_counts$N)
  }
    
  max_prob <- max(unlist(lapply(names(deg_eqtl_list), function(tissue) lapply(1:2, function(sex) probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]])))) + 1
    
  nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
  grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/byPhenotype_colocalization_summary_", 
                                         length(coloc_phenotypes), "traits_", paste0(c(1,2,4,8), "w")[timepoint],"_nicole", ifelse(arnold_mods, "_AB", ""),".pdf"), 
                       width = 1500 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
  par(mfrow = c(1, 1), mar = c(4,3,2,3), xpd = NA)
  
  plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
  ylocs <- seq(10, 0, length.out = length(coloc_phenotypes))

  #put boxes around categories, if we're doing that  
  if(partition_by_category){
    segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
    for(traitcat in 1:(nrow(trait_category_counts)-1)){
      segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
    }
    catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
    catnames <- as.vector(trait_category_counts$V1)
    catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
    catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
    catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
    catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
    catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/3
    text(x = -0.1, y = catylocs, labels = catnames, pos = 2, srt = 90, col = "grey35")
  }
  
  #plot names of traits
  if(change_names_in_plot){
    coloc_phenotypes_newnames <- traitname_map[match(coloc_phenotypes, traitname_map[,1]),2]
  } else {
    coloc_phenotypes_newnames <- coloc_phenotypes
  }
  
  if(nicole_mods){
    text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = 1, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  } else {
    text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
  }
  
  text(labels = "Phenotypes", x = nameloc, y = 10.15, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
  text(labels = "Posterior Probability of Colocalization in at Least One Gene", x = 1.325, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
  
  #guiding lines for traits
  if(nicole_mods){
   segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5) 
  }
  
  #axes
  max_horiz_axis <- 2.05
  segments(x0 = nameloc, x1 = nameloc, y0 = 10.2, y1 = -0.1, lwd = 2)
  segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
  
  #horiz axis ticks and nums
  segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
           x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
           y0 = -0.1, y1 = -0.15, lwd = 2)
  text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
       y = -0.135, pos = 1, cex = 1.25,
       labels =  c(0, sapply(seq(1, max_prob, by = 1), function(num) as.expression(bquote(1-10^-.(num))))))
  
  #let's get minor tick marks in there too
  minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 1))))
  minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
  minticks <- minticks[minticks < max_horiz_axis]
  segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
  
  #dashed lines for ticks
  if(nicole_mods){
    segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
           x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
  } else {
    segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
             x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
  }
  #plot legend for colors etc.
  if(nicole_mods){rect(xleft = 1.97, ybottom = 8.35, ytop = 10.1, xright = 2.2, border = NA, col = "white")}
  text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
       y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
  points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
  
  text(labels = c("Male", "Female"),
       y = seq(8.5, 8.425, length.out = 2), x = 2, pos = 4)
  points(x = rep(1.9925, 2), y = seq(8.505, 8.43, length.out = 2), col = cols$Sex, pch = c(18, 16), cex = 2)
  
  for(tissue in names(deg_eqtl_list)){
    
    # grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/basic_figure_", tissue,".pdf"), width = 2000 / 72, height = 10000 / 72 / length(names(deg_eqtl_list)))
    
    
    
    # if(tissue != names(deg_eqtl_list)[1]){break}
    
    # #legend for dots and dashes
    # stponr <- c(99) #set_to_put_on_first_row
    # legend(x = 3.5155, y = 1.194 - ifelse(any(tissue == names(deg_eqtl_list)[stponr]), 0, 1.00225), col = c(1,"lightgray"), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
    #        legend = c("Point Estimate Ratio (log2-scale)", "Same Magnitude DE as LFC"), cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
    # 
    # rect(xleft = 0, ybottom = 0, xright = 4, ytop = 2.15, lwd = 3) #whole plot
    # rect(xleft = 0, ybottom = 2, xright = 4, ytop = 2.15, col = rgb(0,0,0.1, alpha = 0.05), border = NA) #whole plot
    # segments(x0 = 0, x1 = 4, y0 = 2, y1 = 2, lwd = 3) #week title
    # segments(x0 = 0, x1 = 4, y0 = 1, y1 = 1, lwd = 2) #sex divisor
    # shadowtext(x = 2, y = 2.15, labels = stringr::str_to_title(paste0(strsplit(tissue, "-")[[1]][-1], collapse = " ")), 
    #            cex = 5, col = cols$Tissue[tissue], pos = 3, r = 0.2) #timepoint labels
    # mtext(text = paste0("Phenotype Indices for the ", length(coloc_phenotypes), " GWAS Phenotypes w/ Colocalization Probabilities > ", pp4_threshold), cex = 2, line = 0, side = 1) #horiz axis label
    # text(labels = "trait", cex = 1.5, x = -0.015, y = 0.05, font = 3, srt = 90) #horiz axis label
    # text(labels = "trait", cex = 1.5, x = -0.015, y = 1.05, font = 3, srt = 90) #horiz axis label
    # text(labels = latex2exp::TeX("-log_{10}Pr(No Colocs)"), cex = 2.5, x = -0.075, y = 0.565, font = 3, srt = 90) #horiz axis label
    # text(labels = latex2exp::TeX("-log_{10}Pr(No Colocs)"), cex = 2.5, x = -0.075, y = 1.565, font = 3, srt = 90) #horiz axis label
    # # text(labels = latex2exp::TeX(paste0("-log_{10}\\prod_{}^{",length(pheno_cols),"}(1-PP4)")), cex = 2.5, x = -0.1, y = 1.6, font = 3, srt = 90) #horiz axis label
    
    #note filtration scheme
    # text(labels = latex2exp::TeX("Results filtered at MoTrPAC FDR-$\\alpha$ = "), 
    #      x = 3.3775, y = 2.18, pos = 4, cex = 1.75)
    # text(labels = sign_filter_alpha, x = 3.91, y = 2.18, pos = 4, cex = 1.75, font = 2)
    
    for(sex in 1:2){
      # for(sex in 1){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
        # for(timepoint in 1){
        
        # #plot sex symbol
        # text(labels = c("\u2642", "\u2640")[sex], x = timepoint - 1.01, y = 1 + sex - ifelse(sex == 1, 1.0675, 1.0825), pos = 4, cex = 4, col = cols$Sex[sex])
        
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          next()
        }
        
        #plot KDE for ratios
        
        #add vertical line to mark zero
        # segments(y0 = sex - 1, y1 = sex, lwd = 2, col = "lightgray", lty = 2,
        #          x0 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb,
        #          x1 = (0 / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb)
        
        # #timepoint axes
        # if(sex == 1 | (tissue == "t64-ovaries" & sex == 2)){
        #   if(any(slope_rescale > 0) & any(slope_rescale < 0)){
        #     tick_vals <- round(seq(from = 0, to = max(abs(slope_rescale)) * sign(slope_rescale[which.max(abs(slope_rescale))]), length.out = 7), 2)  
        #     tick_vals <- c(rev(seq(from = 0, to = min(abs(slope_rescale)) * sign(slope_rescale[which.min(abs(slope_rescale))]), 
        #                            by = -diff(tick_vals)[1])), tick_vals)
        #     tick_vals[-which(diff(tick_vals) < 1E-6)]
        #   } else {
        #     tick_vals <- round(seq(from = slope_rescale[1], to = slope_rescale[2], length.out = 10), 2)  
        #   }
        #   tick_locs <- (tick_vals / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb
        #   axis(1, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -2)
        #   mtext(text = tick_vals, side = 1, at = tick_locs, cex = 1, line = -0.5)
        # }
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        trait_probs <- probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]]
        trait_probs <- trait_probs[order_traits]
        
        points(x = trait_probs / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = trait_locs, pch = c(18, 16)[sex], col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 0.75), cex = 2)
                 
        
        
        # #sex axes
        # if(timepoint == 1 & !(tissue == "t64-ovaries" & sex == 1) & !(tissue == "t63-testes" & sex == 2)){
        #   # tick_vals <- round(seq(from = logFC_rescale[1], to = logFC_rescale[2], length.out = 5), 1)
        #   tick_vals <- round(seq(from = 0, to = logFC_rescale[2], length.out = 4), 1)
        #   tick_locs <- (tick_vals / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb
        #   axis(2, at = tick_locs, labels = rep("", length(tick_vals)), lwd = 2, cex.axis = 2, tck = -0.015, line = -7.75)
        #   mtext(text = tick_vals, side = 2, at = tick_locs, cex = 1, line = -6.75)
        # }
        
        ### label genes
        # if(length(delsub_colocs$feature_ID) != 0){
        #   n_genes_to_label_by_PP4 <- 3
        #   ptl <- order((delsub_colocs$pp4), decreasing = T)[1:n_genes_to_label_by_PP4] #points to label
        #   text(x = (delsub_colocs$abs_slope[ptl] / diff(slope_rescale) - slope_rescale[1] / diff(slope_rescale)) * (1-2*tarb) + timepoint - 1 + tarb - 0.01,
        #        y = (-(log10(1-delsub_colocs$pp4[ptl])) / diff(logFC_rescale) - logFC_rescale[1] / diff(logFC_rescale)) * (1-2*tarb) + sex - 1 + tarb + 0.0125,
        #        labels = delsub$gene_name[ptl], col = grDevices::adjustcolor(p_val_cols[ceiling(log(delsub$pval_beta[ptl]) / min_logpval * 100)], alpha.f = 0.75),
        #        cex = 1.5, pos = 4, srt = 65)
        # }
        
        
      }
      
  }
  
  if(arnold_mods){
    addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
  }
  
  dev.off()
  }
  
  

}


#### Stephen's Figure -- the ratio of colocs on either side of the line ####
plot_ratio_colocs_DE_inequality <- T
change_names_in_plot = T

tissues <- names(deg_eqtl_list)
pp4_threshold <- 1e-4
if(!exists("coloc_list") & !exists("colocs")){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
  gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
  gonocs$motrpac_tissue <- "t1000-gonads"
  colocs <- rbind(colocs, gonocs)
  coloc_phenotypes <- sort(unique(colocs$gwas_trait))
  colocs <- lapply(coloc_phenotypes, function(coloc_phenotype) colocs[colocs$gwas_trait == coloc_phenotype,])
  names(colocs) <- coloc_phenotypes
  for(coloc_phenotype in coloc_phenotypes){
    print(coloc_phenotype)
    colocs[[coloc_phenotype]] <-  lapply(tissues, function(tissue) colocs[[coloc_phenotype]][colocs[[coloc_phenotype]]$motrpac_tissue == tissue,])
    names(colocs[[coloc_phenotype]]) <- tissues
  }
}

#get names and categories
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]


if(!exists(left_or_right_of_line) & !file.exists(paste0("~/Documents/DExEQTL/LoRoL_probs_", pp4_threshold, "pp4_", discr_coloc_PA_filter, "discrColocPAfilter"))){
left_or_right_of_line <- array(data = 0, dim = c(2,4,length(names(deg_eqtl_list)),length(coloc_phenotypes),4), 
                               dimnames = list(sex = c("male", "female"), 
                                               timepoint = paste0(c(1,2,4,8), "w"),
                                               tissue = names(deg_eqtl_list),
                                               coloc_phenotype = coloc_phenotypes,
                                               side = c("n_left", "n_right", "sum_log_1-pp4_left", "sum_log_1-pp4_right")))
sign_filter_alpha <- 0.1
discr_coloc_PA_filter <- 0.1
for(tissue in names(deg_eqtl_list)){

  for(sex in 1:2){

    sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
    
    for(timepoint in 1:4){

      timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
      
      if(length(deg_eqtl_list[[tissue]]$selection_fdr) > 0){
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
      } else {
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
      }
      
      delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
      delsub$abs_logFC_minus_abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) #* sign(delsub$logFC)
      
      if((sex == 1 & tissue == "t64-ovaries") | (sex == 2 & tissue == "t63-testes" | nrow(delsub) == 0)){
        next()
      }
      
      #snag colocalizing genes
      print(tissue)
      for(coloc_phenotype in coloc_phenotypes){
        cat(paste0(timepoint, " "))
        delsub$pp4 <- sapply(delsub$gene_id, function(cg) colocs[[coloc_phenotype]][[tissue]]$p4[colocs[[coloc_phenotype]][[tissue]]$gene_id == cg][1]) 
        delsub$pp4[delsub$pp4 > (1 - 1E-6)] <- (1 - 1E-6)
        delsub_colocs <- delsub[!is.na(delsub$pp4),]
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "n_left"] <- sum(delsub_colocs$abs_logFC_minus_abs_slope < 0 & delsub_colocs$pp4 > discr_coloc_PA_filter)
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "n_right"] <- sum(delsub_colocs$abs_logFC_minus_abs_slope > 0 & delsub_colocs$pp4 > discr_coloc_PA_filter)
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_left"] <- sum(-log10(1-delsub_colocs$pp4[(delsub_colocs$abs_logFC_minus_abs_slope < 0)]))
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_right"] <- sum(-log10(1-delsub_colocs$pp4[(delsub_colocs$abs_logFC_minus_abs_slope > 0)]))
        # print(left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_right"] - 
        #       left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_left"])
      
      }
      
      
    }
    
  }
  # dev.off()
}
save(left_or_right_of_line, file = paste0("~/Documents/DExEQTL/LoRoL_probs_", pp4_threshold, "pp4_", discr_coloc_PA_filter, "discrColocPAfilter"))
}

if(!exists("left_or_right_of_line")){
  load(paste0("~/Documents/DExEQTL/LoRoL_probs_", pp4_threshold, "pp4_", discr_coloc_PA_filter, "discrColocPAfilter"))
}
#compute difference in log-odds of at least one gene on left side vs right side
#or just diff in log prob
#can also do chi2 test for difference in incidence

sign_filter_alpha <- 0.1
if(!exists("n_deGenes_lorol") | sign_filter_alpha != sigfilt_used){
  n_deGenes_lorol <- array(data = 0, dim = c(2,4,length(names(deg_eqtl_list)),2), 
                     dimnames = list(sex = c("male", "female"), 
                                     timepoint = paste0(c(1,2,4,8), "w"),
                                     tissue = names(deg_eqtl_list),
                                     lorol = c("left", "right")
                     ))
  
  for(sex in c("male", "female")){
    for(timepoint in paste0(c(1,2,4,8), "w")){
      for(tissue in names(deg_eqtl_list)){
        # motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == sex)
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == timepoint)
        left_inds <- which( (abs(deg_eqtl_list[[tissue]]$logFC) - abs(deg_eqtl_list[[tissue]]$abs_slope)) < 0 )
        right_inds <- which( (abs(deg_eqtl_list[[tissue]]$logFC) - abs(deg_eqtl_list[[tissue]]$abs_slope)) > 0 )
        n_deGenes_lorol[sex, timepoint, tissue, "left"] <- (length(intersect(intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds), left_inds)))
        n_deGenes_lorol[sex, timepoint, tissue, "right"] <- (length(intersect(intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds), right_inds)))
        
        if(!(n_deGenes_lorol[sex, timepoint, tissue, "right"] + n_deGenes_lorol[sex, timepoint, tissue, "left"] == (length((intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds)))))){
          print("whoops")
        }
      }
    }
  }
  sigfilt_used <- sign_filter_alpha
  n_deGenes_lorol[n_deGenes_lorol == 0] <- 1 #make division play nice later
}

#specify colors
coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)
padded_nums <- as.character(1:length(coloc_phenotypes))
padded_nums <- (sapply(padded_nums, function(pn) paste0(paste0(rep("0", max(nchar(padded_nums)) - nchar(pn)), collapse = ""), pn)))
cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"



plot_phenotype_summary_figure_nicole = T
nicole_mods <- T
arnold_mods <- F
arnold_bayes <- png::readPNG(source = "~/Pictures/ArnoldBayes.png")
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))

reorder_vertical <- T
partition_by_category <- T
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))

use_geometric_mean <- T

if(plot_ratio_colocs_DE_inequality){
  
  max_horiz_axis <- 2.05
  
  for(timepoint in 1:4){
    
    if(reorder_vertical){
      order_traits <- order(sapply(coloc_phenotypes, function(coloc_phenotype) sum(left_or_right_of_line[,,,coloc_phenotype,4]) - sum(left_or_right_of_line[,,,coloc_phenotype,3])), decreasing = T)
      
    } else {
      order_traits <- 1:length(coloc_phenotypes)
    }
    
    if(partition_by_category){
      order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])], sort(order_traits))]
      trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])] ))
      trait_category_counts$cN <- cumsum(trait_category_counts$N)
    }
    
    if(use_geometric_mean){
      logprob_range <- range(as.vector(sapply(coloc_phenotypes, function(coloc_phenotype) sapply(1:2, function(sex) sapply(names(deg_eqtl_list), function(tissue) 
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_right"]  / n_deGenes_lorol[sex, timepoint, tissue, "right"] - 
          left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_left"]  / n_deGenes_lorol[sex, timepoint, tissue, "left"])))))
    } else {
      logprob_range <- range(as.vector(sapply(coloc_phenotypes, function(coloc_phenotype) sapply(1:2, function(sex) sapply(names(deg_eqtl_list), function(tissue) 
        left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_right"] - left_or_right_of_line[sex, timepoint, tissue, coloc_phenotype, "sum_log_1-pp4_left"])))))
    }
    max_prob <- diff(logprob_range) + ifelse(use_geometric_mean, 0.01, 1)
    

    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/LoRoL_", ifelse(use_geometric_mean, "geommean_", ""),
                                           length(coloc_phenotypes), "traits_", paste0(c(1,2,4,8), "w")[timepoint], ifelse(arnold_mods, "_AB", ""),".pdf"), 
                         width = 1500 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
    par(mfrow = c(1, 1), mar = c(4,3,3,3), xpd = NA)
    
    
    plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    ylocs <- seq(10, 0, length.out = length(coloc_phenotypes))
    
    #put boxes around categories, if we're doing that  
    if(partition_by_category){
      segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
      for(traitcat in 1:(nrow(trait_category_counts)-1)){
        segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
      }
      catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
      catnames <- as.vector(trait_category_counts$V1)
      catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
      catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
      catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
      catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
      catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/3
      text(x = -0.1, y = catylocs, labels = catnames, pos = 2, srt = 90, col = "grey35")
    }
    
    #plot names of traits
    if(change_names_in_plot){
      coloc_phenotypes_newnames <- traitname_map[match(coloc_phenotypes, traitname_map[,1]),2]
    } else {
      coloc_phenotypes_newnames <- coloc_phenotypes
    }
    
    if(nicole_mods){
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = 1, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    
    text(labels = "Phenotypes", x = nameloc, y = 10.15, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    if(use_geometric_mean){
      text(labels = "Ratio in Geometric Mean Posterior Probability of Colocalization", x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
    } else {
      text(labels = "Ratio in Posterior Probability of Colocalization in at Least One Gene", x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)  
    }
    
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5) 
    }
    
    #axes
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.2, y1 = -0.1, lwd = 2)
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
    
    #horiz axis ticks and nums
    if(use_geometric_mean){
      
      probs <- c(-0.2, -0.1, 0, 0.1)
      segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               y0 = -0.1, y1 = -0.15, lwd = 2)
      text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           y = -0.135, pos = 1, cex = 1.25,
           labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
      
    } else {
      
      segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               y0 = -0.1, y1 = -0.15, lwd = 2)
      text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
           y = -0.135, pos = 1, cex = 1.25,
           labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
    }
    
    
    #let's get minor tick marks in there too
    if(use_geometric_mean){
      minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
      minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
      minticks <- minticks[minticks < max_horiz_axis]
      #hmm being weird
    } else {
      minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
      minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
      minticks <- minticks[minticks < max_horiz_axis]
      segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    }
    
    
    #dashed lines for ticks
    if(use_geometric_mean){
      if(nicole_mods){
        segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
                 x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
      } else {
        segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
                 x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
      }
    } else {
      if(nicole_mods){
        segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
                 x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
      } else {
        segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
                 x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
      }  
    }
    
    
    #helpful hint line about direction
    
    if(use_geometric_mean){
      segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
      text(labels = "eQTL > DEG", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
      text(labels = "DEG > eQTL", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
      text(labels = "Stronger Coloc\nSignal When", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    } else {
      segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
      text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
      text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
      text(labels = "Stronger Coloc\nSignal When", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    }
    
    #plot legend for colors etc.
    if(nicole_mods){rect(xleft = 1.97, ybottom = 8.35, ytop = 10.1, xright = 2.2, border = NA, col = "white")}
    text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
         y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    
    text(labels = c("Male", "Female"),
         y = seq(8.5, 8.425, length.out = 2), x = 2, pos = 4)
    points(x = rep(1.9925, 2), y = seq(8.505, 8.43, length.out = 2), col = cols$Sex, pch = c(18, 16), cex = 2)
    
    for(tissue in names(deg_eqtl_list)){
      
      for(sex in 1:2){

        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          next()
        }
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        
        if(use_geometric_mean){
          trait_probs <- left_or_right_of_line[sex, timepoint, tissue,,"sum_log_1-pp4_right"] / n_deGenes_lorol[sex, timepoint, tissue, "right"] - 
            left_or_right_of_line[sex, timepoint, tissue,,"sum_log_1-pp4_left"] / n_deGenes_lorol[sex, timepoint, tissue, "left"]
        } else {
          trait_probs <- left_or_right_of_line[sex, timepoint, tissue,,"sum_log_1-pp4_right"]- left_or_right_of_line[sex, timepoint, tissue,,"sum_log_1-pp4_left"]
        }
        trait_probs <- trait_probs[order_traits]
        
        points(x = (trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = trait_locs, pch = c(18, 16)[sex], col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 0.75), cex = 2)
        
        
        
      }
      
    }
    
    if(arnold_mods){
      addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
    }
    
    dev.off()
  }
  
  
  
}


#### find number of DE genes per tissue ####

# system("gsutil cp gs://mawg-data/pass1b-06/transcript-rna-seq/dea/transcript_rna_seq_20210126.RData /tmp/")
# load("/tmp/transcript_rna_seq_20210126.RData")
# timewise_dea = data.table(transcript_rna_seq$timewise_dea)
# degs = timewise_dea[selection_fdr < 0.1]
# degs_8w = degs[comparison_group == '8w' & adj_p_value < 0.1]
# table(degs_8w[,tissue_abbreviation], degs_8w[,sex])

sign_filter_alpha <- 0.1
use_selection_fdr_too <- F
if(!exists("n_deGenes") | sign_filter_alpha != sigfilt_used | use_selection_fdr_too != selection_FDR_used){
n_deGenes <- array(data = 0, dim = c(2,4,length(names(deg_eqtl_list))), 
                               dimnames = list(sex = c("male", "female"), 
                                               timepoint = paste0(c(1,2,4,8), "w"),
                                               tissue = names(deg_eqtl_list)
                                               ))
for(sex in c("male", "female")){
  for(timepoint in paste0(c(1,2,4,8), "w")){
    for(tissue in names(deg_eqtl_list)){
      if(use_selection_fdr_too){
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)  
      } else {
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)  
      }
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == sex)
      timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == timepoint)
      n_deGenes[sex, timepoint, tissue] <- (length(intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds)))
    }
  }
}
sigfilt_used <- sign_filter_alpha
selection_FDR_used <- use_selection_fdr_too
}

t(n_deGenes[,"8w",])

grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/num_DE_genes_across_timepoints_",ifelse(use_selection_fdr_too, "selectionFDR+adjPVAL_used_", "just_adjPVAL_used_"),"atAlpha",sign_filter_alpha,".pdf"), 
                     width = 10, height = 4.5, family="Arial Unicode MS")

par(mar = c(5,6,3,8))
plot(1,1,xlim = range(n_deGenes), ylim = c(0.5,4.5), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#axis labels and title
text(labels = "Timepoint", x = -0.144  * max(n_deGenes), y = 4.25, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25", srt = 90)
text(labels = "# Genes", x = max(n_deGenes)/2, y = -0.45, pos = 1, cex = 2.5,  xpd = NA, family="Courier", col = "grey25")
text(labels = "# Genes DE at each Timepoint", x = max(n_deGenes)/2, y = 4.75, pos = 3, cex = 2.5,  xpd = NA)

#axes
box("plot", lwd = 2)
shadowtext(x = rep(-0.04 * max(n_deGenes),4), y = 1:4, labels = paste0(c(1,2,4,8), "w"), 
           pos = 2, xpd = NA, cex = 1.5, font = 2, col = cols$Time, r = 0.05)
segments(y0 = 1:4, y1 = 1:4, x0 = -1E6, x1 = 1E6, lwd = 0.5, lty = 3)
xlocs <- round(seq(range(n_deGenes)[1], range(n_deGenes)[2], by = 50))
segments(x0 = xlocs, x1 = xlocs, y0 = 0.325, y1 = 0.2, lwd = 2, xpd = NA)
text(x = xlocs, y = 0.2, pos = 1, cex = 1, labels = xlocs, xpd = NA)
#horiz axis ticks and nums

#plot points themselves
for(sex in 1:2){
  for(timepoint in 1:4){
    for(tissue in names(deg_eqtl_list)){
      points(x = n_deGenes[sex, timepoint, tissue], y = timepoint, cex =2,
             pch = c(18, 16)[sex], col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 0.95))
    }
  }
}

text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
     y = seq(4.5,1.175,length.out = length(names(deg_eqtl_list))), x = 1.075 * max(n_deGenes), pos = 4, xpd = NA, cex = 0.9)
points(x = rep(1.0674 * max(n_deGenes), length(names(deg_eqtl_list))), y = seq(4.55, 1.1755,length.out = length(names(deg_eqtl_list))), 
       col = cols$Tissue, pch = 15, cex = 1.85, xpd = NA)

text(labels = c("Male", "Female"),
     y = seq(0.75,0.5,length.out = 2), x = 1.075 * max(n_deGenes), pos = 4, xpd = NA, cex = 0.9)
points(x = rep(1.0674 * max(n_deGenes), 2), y = seq(0.75,0.5,length.out = 2), 
       col = cols$Sex, pch = c(18, 16), cex = 2, xpd = NA)

dev.off()

#make barplot of just the 8w
grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/num_DE_genes_across_timepoints_8w_",ifelse(use_selection_fdr_too, "selectionFDR+adjPVAL_used_", "just_adjPVAL_used_"),"atAlpha_",sign_filter_alpha,".pdf"), 
                     width = 15, height = 10, family="Arial Unicode MS")

par(mar = c(8,6,3,2))
plot(1,1,xlim = c(0,length(names(deg_eqtl_list))+1), ylim = range(n_deGenes), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")

#axis labels and title
text(labels = "# Genes", x = max(n_deGenes)/2, y = -0.45, pos = 1, cex = 2.5,  xpd = NA, family="Courier", col = "grey25")
text(labels = "# Genes DE in each Tissue", x = max(n_deGenes)/2, y = 4.75, pos = 3, cex = 2.5,  xpd = NA)

#axes
ylocs <- round(seq(range(n_deGenes)[1], range(n_deGenes)[2], by = 50))
xlocs <- 1:length(names(deg_eqtl_list)); names(xlocs) = names(deg_eqtl_list)[order(apply(n_deGenes[,"8w",],2,sum), decreasing = T)]

segments(x0 = xlocs, x1 = xlocs, y0 = 0.325, y1 = 0.2, lwd = 2, xpd = NA)
text(x = xlocs + 0.25, y = -12.5 / 500 * max(n_deGenes), pos = 2, cex = 1.25, labels = sapply(names(xlocs), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))), 
     xpd = NA, srt = 45, col = cols$Tissue[names(xlocs)])
text(labels = c("\u2642", "\u2640")[1], x = xlocs - 0.25, y = 4, pos = 1, cex = 1, col = cols$Sex[1], font = 2)
text(labels = c("\u2642", "\u2640")[2], x = xlocs + 0.25, y = 4, pos = 1, cex = 1, col = cols$Sex[2], font = 2)
text(x = 10, y = max(n_deGenes), cex = 3, labels = "8w Timepoint", col = cols$Time["8w"])

#vertical axis
segments(x0 = 0, x1 = 0, y0 = 0, y1 = max(n_deGenes), lwd = 3, xpd = NA)
segments(x0 = 0, x1 = -0.25, y0 = ylocs, y1 = ylocs, lwd = 2, xpd = NA)
text(x = -0.125, y = ylocs, labels = ylocs, pos = 2, xpd = NA)
text(x = -1.5, y = max(n_deGenes) * 0.95, labels = paste0("# genes differentially expressed at FDR = ", sign_filter_alpha), pos = 2, xpd = NA, srt = 90, cex = 2)
#horiz axis ticks and nums

#plot rects themselves
width_room <- 0.035
for(sex in 1:2){
  for(timepoint in 4){
    for(tissue in names(deg_eqtl_list)){
      rect(xleft = xlocs[tissue] - ifelse(sex == 1, 0.5 - width_room/2, - width_room/2), xright  = xlocs[tissue] + ifelse(sex == 1, - width_room/2, 0.5 - width_room/2), border = grDevices::adjustcolor(cols$Sex[sex], alpha.f = 1), lwd = 2,
           ytop = n_deGenes[sex, timepoint, tissue], ybottom = 0, col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 1))
      text(x = xlocs[tissue] - ifelse(sex == 1, (0.5-width_room)/2, -(0.5-width_room)/2), y = n_deGenes[sex, timepoint, tissue] - max(n_deGenes) / 125, labels = n_deGenes[sex, timepoint, tissue], cex = 0.75, pos = 3)
    }
  }
}

text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
     y = seq(4.5,1.175,length.out = length(names(deg_eqtl_list))), x = 1.075 * max(n_deGenes), pos = 4, xpd = NA, cex = 0.9)
points(x = rep(1.0674 * max(n_deGenes), length(names(deg_eqtl_list))), y = seq(4.55, 1.1755,length.out = length(names(deg_eqtl_list))), 
       col = cols$Tissue, pch = 15, cex = 1.85, xpd = NA)

text(labels = c("Male", "Female"),
     y = seq(0.75,0.5,length.out = 2), x = 1.075 * max(n_deGenes), pos = 4, xpd = NA, cex = 0.9)
points(x = rep(1.0674 * max(n_deGenes), 2), y = seq(0.75,0.5,length.out = 2), 
       col = cols$Sex, pch = c(18, 16), cex = 2, xpd = NA)

dev.off()




#### summarizing by trait figure -- geometric mean edition ####

sign_filter_alpha <- 0.1
if(!exists("n_deGenes") | sign_filter_alpha != sigfilt_used){
  n_deGenes <- array(data = 0, dim = c(2,4,length(names(deg_eqtl_list))), 
                     dimnames = list(sex = c("male", "female"), 
                                     timepoint = paste0(c(1,2,4,8), "w"),
                                     tissue = names(deg_eqtl_list)
                     ))
  
  for(sex in c("male", "female")){
    for(timepoint in paste0(c(1,2,4,8), "w")){
      for(tissue in names(deg_eqtl_list)){
        # motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == sex)
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == timepoint)
        n_deGenes[sex, timepoint, tissue] <- (length(intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds)))
      }
    }
  }
  sigfilt_used <- sign_filter_alpha
  n_deGenes[n_deGenes == 0] <- 1 #make division play nice later
}

change_names <- F
change_names_in_plot = T
pp4_threshold <- 0.01
if(!exists("coloc_list") & !exists("colocs") | last_import_pp4t != pp4_threshold | change_names != names_changed){
  load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
  colocs <- coloc_list
  rm(coloc_list)
  gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
  gonocs$motrpac_tissue <- "t1000-gonads"
  colocs <- rbind(colocs, gonocs)
  coloc_phenotypes <- sort(unique(colocs$gwas_trait))
  last_import_pp4t <- pp4_threshold
  trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
  if(change_names){
    traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
    colocs$gwas_trait <- traitname_map[match(colocs$gwas_trait, traitname_map[,1]),2]
    coloc_phenotypes <- sort(unique(colocs$gwas_trait))
    names_changed = T
  } else {
    names_changed = F
  }
}

load(file = paste0("~/Documents/DExEQTL/probs_noColocs_byTrait_", pp4_threshold, "_threshold"))

coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)
padded_nums <- as.character(1:length(coloc_phenotypes))
padded_nums <- (sapply(padded_nums, function(pn) paste0(paste0(rep("0", max(nchar(padded_nums)) - nchar(pn)), collapse = ""), pn)))

plot_phenotype_summary_figure_geommean = T
nicole_mods <- T
arnold_mods <- F
partition_by_category <- T
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))
arnold_bayes <- png::readPNG(source = "~/Pictures/ArnoldBayes.png")
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))
reorder_vertical <- T

if(change_names){
  for(i1 in seq_along(probs_noColocs)){
    for(i2 in seq_along(probs_noColocs[[1]])){
      for(i3 in seq_along(probs_noColocs[[1]][[1]])){
        names(probs_noColocs[[i1]][[i2]][[i3]]) <- traitname_map[match(names(probs_noColocs[[i1]][[i2]][[i3]]), traitname_map[,1]),2]
      }
    }
  }
}


if(plot_phenotype_summary_figure_geommean){
  
  for(timepoint in 1:4){
    
    if(reorder_vertical){
      order_traits <- order(apply(sapply(coloc_phenotypes, (function(coloc_phenotype) sapply(1:2, function(sex) sapply(1:4, function(timepoint) sapply(names(deg_eqtl_list), function(tissue) 
        probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]][[coloc_phenotype]] / n_deGenes[sex, timepoint,tissue]
      ))))), 2, sum), decreasing = T)
    } else {
      order_traits <- 1:length(coloc_phenotypes)
    }
    
    if(partition_by_category){
      order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])], sort(order_traits))]
      trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])] ))
      trait_category_counts$cN <- cumsum(trait_category_counts$N)
    }
    
    max_prob <- max(unlist(lapply(names(deg_eqtl_list), function(tissue) lapply(1:2, function(sex) 
      probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]] / n_deGenes[sex, timepoint,tissue]))))
    max_prob <- max_prob + 10^round(log(max_prob))
    
    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/byPhenotype_colocalization_summary_geommean_", 
                                           length(coloc_phenotypes), "traits_", paste0(c(1,2,4,8), "w")[timepoint],"_nicole", ifelse(arnold_mods, "_AB", ""),".pdf"), 
                         width = 1500 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
    par(mfrow = c(1, 1), mar = c(4,3,2,3), xpd = NA)
    
    plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    ylocs <- seq(10, 0, length.out = length(coloc_phenotypes))
    
    #put boxes around categories, if we're doing that  
    if(partition_by_category){
      segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
      for(traitcat in 1:(nrow(trait_category_counts)-1)){
        segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
      }
      catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
      catnames <- as.vector(trait_category_counts$V1)
      catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
      catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
      catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
      catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
      catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/3
      text(x = -0.1, y = catylocs, labels = catnames, pos = 2, srt = 90, col = "grey35")
    }
    
    #plot names of traits
    if(change_names_in_plot){
      coloc_phenotypes_newnames <- traitname_map[match(coloc_phenotypes, traitname_map[,1]),2]
    } else {
      coloc_phenotypes_newnames <- coloc_phenotypes
    }
    
    if(nicole_mods){
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = 1, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    text(labels = "Phenotypes", x = nameloc, y = 10.15, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    text(labels = "Posterior Probability of Colocalization (Geometric Mean)", x = 1.225, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5) 
    }
    
    #axes
    max_horiz_axis <- 2.05
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.2, y1 = -0.1, lwd = 2)
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
    
    #horiz axis ticks and nums
    segments(x0 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
             x1 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
             y0 = -0.1, y1 = -0.15, lwd = 2)
    text(x = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
         y = -0.135, pos = 1, cex = 1.25,
         labels =  c(0, sapply(seq(10^round(log(max_prob)), max_prob, by = 10^round(log(max_prob))), function(num) as.expression(bquote(1-10^-.(num))))))
    
    #let's get minor tick marks in there too
    minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
    minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
    minticks <- minticks[minticks < max_horiz_axis]
    segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    
    #dashed lines for ticks
    if(nicole_mods){
      segments(x0 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
    } else {
      segments(x0 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = seq(0, max_prob, by = 10^round(log(max_prob))) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
    }
    #plot legend for colors etc.
    if(nicole_mods){rect(xleft = 1.97, ybottom = 8.35, ytop = 10.1, xright = 2.2, border = NA, col = "white")}
    text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
         y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    
    text(labels = c("Male", "Female"),
         y = seq(8.5, 8.425, length.out = 2), x = 2, pos = 4)
    points(x = rep(1.9925, 2), y = seq(8.505, 8.43, length.out = 2), col = cols$Sex, pch = c(18, 16), cex = 2)
    
    for(tissue in names(deg_eqtl_list)){
      
      
      for(sex in 1:2){
        # for(sex in 1){
        
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        #incompatible gonads message
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2)){
          next()
        }
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        trait_probs <- probs_noColocs[[tissue]][[c("male", "female")[sex]]][[paste0(c(1,2,4,8), "w")[timepoint]]] / n_deGenes[sex, timepoint, tissue]
        trait_probs <- trait_probs[order_traits]
        
        points(x = trait_probs / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = trait_locs, pch = c(18, 16)[sex], col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 0.85), cex = 2)
        
      }
      
    }
    
    if(arnold_mods){
      addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
    }
    
    dev.off()
  }
  
  
  
}


#### quick figure to show gastroc kde ####

gastroc_rats <- list(male = "", female = "")
for(sex in c("male", "female")){
  for(timepoint in paste0("8w")){
    for(tissue in names(deg_eqtl_list)[14]){
      # motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
      motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == sex)
      timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == timepoint)
      delsub <- deg_eqtl_list[[tissue]][intersect(intersect(motrpac_signif_inds, sex_inds), timepoint_inds),]
      gastroc_rats[[sex]] <- abs(delsub$logFC) - abs(delsub$abs_slope)
    }
  }
}

mdens <- density(gastroc_rats$male, bw = 0.2)
fdens <- density(gastroc_rats$female, bw = 0.2)
mdens$y <- mdens$y[mdens$x < 2.5]
mdens$x <- mdens$x[mdens$x < 2.5]
fdens$y <- fdens$y[fdens$x < 2.5]
fdens$x <- fdens$x[fdens$x < 2.5]


par(mar = c(4,4,2,2))
plot(1,1,col = "white", xlim = c(-1.3,3), ylim = c(0,1.3), 
     xlab = latex2exp::TeX("Abs. Differential Expression ÷ Abs. Lead eVariant eQTL Effect Size (log_{2}ratio)"),
     ylab = "Density", cex.lab = 1.25)
lines(mdens$x, mdens$y, lwd = 2, col = cols$Sex["male"])
lines(fdens$x, fdens$y, lwd = 2, col = cols$Sex["female"])
polygon(x = mdens$x, y = mdens$y, col = grDevices::adjustcolor(cols$Sex["male"], alpha.f = 0.35))
polygon(x = fdens$x, y = fdens$y, col = grDevices::adjustcolor(cols$Sex["female"], alpha.f = 0.35))
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 10, lwd = 2, lty = 2, col = "grey50")
shadowtext(x = 0, y = 0.2, labels = round(sum(gastroc_rats$female > 0) / length(gastroc_rats$female) * 100), 
     col = cols$Sex["female"], cex = 1.5, pos = 4, font = 2)
shadowtext(x = 0, y = 0.2, labels = round(sum(gastroc_rats$female < 0) / length(gastroc_rats$female) * 100), 
     col = cols$Sex["female"], cex = 1.5, pos = 2, font = 2)
shadowtext(x = 0, y = 1, labels = round(sum(gastroc_rats$male > 0) / length(gastroc_rats$male) * 100), 
     col = cols$Sex["male"], cex = 1.5, pos = 4, font = 2)
shadowtext(x = 0, y = 1, labels = round(sum(gastroc_rats$male < 0) / length(gastroc_rats$male) * 100), 
     col = cols$Sex["male"], cex = 1.5, pos = 2, font = 2)
legend(x = "topright", legend = c("males", "females"), col = cols$Sex, pch = 15, cex = 1.5, box.lty = 3, box.lwd = 2, ncol = 1)
box(which = "plot", lwd = 2)


#### comparing coloc probabilities in DE vs. non-DE genes ####

plot_delta_coloc_DEvsNotDE <- T
change_names = T
change_names_in_plot = T

tissues <- names(deg_eqtl_list)
pp4_threshold <- 0
if(pp4_threshold == 0 & !exists("colocs")){
  load(file = "~/data/smontgom/coloc_list_of_lists_all_p4-0") #loads 'colocs'
  coloc_phenotypes <- sort(names(colocs))
} else {
  if(!exists("coloc_list") & !exists("colocs")){
    load(paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
    colocs <- coloc_list
    rm(coloc_list)
    gonocs <- colocs[colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[1] | colocs$motrpac_tissue == c("t63-testes", "t64-ovaries")[2],]
    gonocs$motrpac_tissue <- "t1000-gonads"
    colocs <- rbind(colocs, gonocs)
    coloc_phenotypes <- sort(unique(colocs$gwas_trait))
    colocs <- lapply(coloc_phenotypes, function(coloc_phenotype) colocs[colocs$gwas_trait == coloc_phenotype,])
    names(colocs) <- coloc_phenotypes
    for(coloc_phenotype in coloc_phenotypes){
      print(coloc_phenotype)
      colocs[[coloc_phenotype]] <-  lapply(tissues, function(tissue) colocs[[coloc_phenotype]][colocs[[coloc_phenotype]]$motrpac_tissue == tissue,])
      names(colocs[[coloc_phenotype]]) <- tissues
    }
  }
}
#get names and categories
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]

sign_filter_alpha <- 0.1
if(!exists("pp4s_DE_notDE") & !file.exists(paste0("~/Documents/DExEQTL/pp4s_DE_notDE_", pp4_threshold, "pp4_alsoRelDEeffects"))){
  sexes <- c("male", "female"); names(sexes) <- sexes
  timepoints <- paste0(c(1,2,4,8), "w"); names(timepoints) <- timepoints
  tissues <- names(deg_eqtl_list); names(tissues) <- tissues
  names(coloc_phenotypes) <- coloc_phenotypes
  pp4s_DE_notDE <- lapply(coloc_phenotypes, function(phenotype) 
    lapply(tissues, function(tissue) 
    lapply(sexes, function(sex) 
    lapply(timepoints, function(timepoint) list(DE_coloc_probs = NA, nDE_coloc_probs = NA, DE_rel_effsiz = NA, nDE_rel_effsiz = NA)
      ))))
    
  for(tissue in names(deg_eqtl_list)){
    
    for(sex in 1:2){
      
      sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
      
      for(timepoint in 1:4){
        
        timepoint_inds <- which(deg_eqtl_list[[tissue]]$comparison_group == paste0(c(1,2,4,8), "w")[timepoint])
        

        motrpac_signif_inds <- which(deg_eqtl_list[[tissue]]$selection_fdr < sign_filter_alpha & deg_eqtl_list[[tissue]]$adj_p_value < sign_filter_alpha)
        motrpac_insignif_inds <- setdiff(1:length(deg_eqtl_list[[tissue]]$selection_fdr), motrpac_signif_inds)

        delsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_signif_inds),]
        antidelsub <- deg_eqtl_list[[tissue]][intersect(intersect(sex_inds, timepoint_inds), motrpac_insignif_inds),]
        delsub$abs_logFC_minus_abs_slope <- (abs(delsub$logFC) - delsub$abs_slope) #* sign(delsub$logFC)
        antidelsub$abs_logFC_minus_abs_slope <- (abs(antidelsub$logFC) - antidelsub$abs_slope) #* sign(delsub$logFC)
        
        if((sex == 1 & tissue == "t64-ovaries") | (sex == 2 & tissue == "t63-testes" | nrow(delsub) == 0)){
          next()
        }
        
        #snag colocalizing genes
        print(tissue)
        for(coloc_phenotype in coloc_phenotypes){
          
          delsub_coloc_inds <- sapply(delsub$gene_id, function(cg) which(colocs[[coloc_phenotype]][[tissue]]$gene_id == cg)[1])
          delsub_coloc_inds <- delsub_coloc_inds[!is.na(delsub_coloc_inds)]
          antidelsub_coloc_inds <- setdiff(1:(length(colocs[[coloc_phenotype]][[tissue]]$gene_id)), delsub_coloc_inds)
          # delsub$pp4[delsub$pp4 > (1 - 1E-6)] <- (1 - 1E-6)
          # delsub <- delsub[!is.na(delsub$pp4),]
          pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["DE_coloc_probs"]] <- colocs[[coloc_phenotype]][[tissue]]$p4[delsub_coloc_inds]
          pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["nDE_coloc_probs"]] <- colocs[[coloc_phenotype]][[tissue]]$p4[antidelsub_coloc_inds]
          pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["DE_rel_effsiz"]] <- delsub$abs_logFC_minus_abs_slope
          pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["nDE_rel_effsiz"]] <- antidelsub$abs_logFC_minus_abs_slope
          
          cat(paste0(timepoint, " (", " " , ")"))
          
        }
        
        
      }
      
    }
    # dev.off()
  }
  save(pp4s_DE_notDE, file = paste0("~/Documents/DExEQTL/pp4s_DE_notDE_", pp4_threshold, "pp4_alsoRelDEeffects"))
}

if(!exists("pp4s_DE_notDE")){
  load(paste0("~/Documents/DExEQTL/pp4s_DE_notDE_", pp4_threshold, "pp4_alsoRelDEeffects"))
}


hist(sapply(coloc_phenotypes, function(phenotype) 
 sapply(tissues, function(tissue) 
  sapply(sexes, function(sex) 
   sapply(timepoints, function(timepoint)
      length(pp4s_DE_notDE[[phenotype]][[tissue]][[sex]][[timepoint]][["DE_coloc_probs"]])
      )))))


#specify colors
coloc_cols = viridis::viridis(50, begin = 0.5, end = 1)
plot(1:length(coloc_cols), 1:length(coloc_cols), col = coloc_cols, pch = 19)
padded_nums <- as.character(1:length(coloc_phenotypes))
padded_nums <- (sapply(padded_nums, function(pn) paste0(paste0(rep("0", max(nchar(padded_nums)) - nchar(pn)), collapse = ""), pn)))
cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

plot_DEvNDE_Figure_colocs <- T
nicole_mods <- T
arnold_mods <- F
arnold_bayes <- png::readPNG(source = "~/Pictures/ArnoldBayes.png")
pheno_cols <- disco::disco("rainbow")
pheno_cols <- rev(colorRampPalette(pheno_cols)(length(coloc_phenotypes)))

reorder_vertical <- T
partition_by_category <- T
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))

#munge into array for efficiency
if(!exists("pp4s_DE_notDE_array")){
  pp4s_DE_notDE_array <- array(data = 0, dim = c(2,4,length(names(deg_eqtl_list)),length(coloc_phenotypes),4), 
                                 dimnames = list(sex = c("male", "female"), 
                                                 timepoint = paste0(c(1,2,4,8), "w"),
                                                 tissue = names(deg_eqtl_list),
                                                 coloc_phenotype = coloc_phenotypes,
                                                 value = c("DE_coloc_probs", "nDE_coloc_probs", "DE_rel_effsiz", "nDE_rel_effsiz")))
  
  for(tissue in names(deg_eqtl_list)){
    for(sex in 1:2){
      for(timepoint in 1:4){
        for(coloc_phenotype in coloc_phenotypes){
          pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "DE_coloc_probs"] <- mean(log(pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["DE_coloc_probs"]]))
          pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "nDE_coloc_probs"] <- mean(log(pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["nDE_coloc_probs"]]))
          pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "DE_rel_effsiz"] <- mean(pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["DE_rel_effsiz"]])
          pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "nDE_rel_effsiz"] <- mean(pp4s_DE_notDE[[coloc_phenotype]][[tissue]][[sex]][[timepoint]][["nDE_rel_effsiz"]])
        }
      }
    }
  }
}

use_geometric_mean <- T

if(plot_DEvNDE_Figure_colocs){
  
  max_horiz_axis <- 2.05
  
  for(timepoint in 1:4){
    
    if(reorder_vertical){
      
      order_traits <- order(sapply(coloc_phenotypes, function(coloc_phenotype) sum(pp4s_DE_notDE_array[,,,coloc_phenotype,1] - pp4s_DE_notDE_array[,,,coloc_phenotype,2], na.rm = T)), decreasing = T)
      
    } else {
      order_traits <- 1:length(coloc_phenotypes)
    }
    
    if(partition_by_category){
      order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])], sort(order_traits))]
      trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes[order_traits], traitwise_partitions[,1])] ))
      trait_category_counts$cN <- cumsum(trait_category_counts$N)
    }
    
    if(use_geometric_mean){
      logprob_range <- range(as.vector(sapply(coloc_phenotypes, function(coloc_phenotype) sapply(1:2, function(sex) sapply(names(deg_eqtl_list), function(tissue) 
        pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "DE_coloc_probs"]  - 
          pp4s_DE_notDE_array[sex, timepoint, tissue, coloc_phenotype, "nDE_coloc_probs"]  )))), na.rm = T)
    } 
    max_prob <- diff(logprob_range) + ifelse(use_geometric_mean, 0.01, 1)
    
    
    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/DEvsNDE_", ifelse(use_geometric_mean, "geommean_", ""),
                                           length(coloc_phenotypes), "traits_", paste0(c(1,2,4,8), "w")[timepoint], ifelse(arnold_mods, "_AB", ""),".pdf"), 
                         width = 1500 / 72, height = 2000 / 72 *18 / 17, family="Arial Unicode MS")
    par(mfrow = c(1, 1), mar = c(4,3,3,3), xpd = NA)
    
    
    plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    ylocs <- seq(10, 0, length.out = length(coloc_phenotypes))
    
    #put boxes around categories, if we're doing that  
    if(partition_by_category){
      segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
      for(traitcat in 1:(nrow(trait_category_counts)-1)){
        segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
      }
      catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
      catnames <- as.vector(trait_category_counts$V1)
      catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
      catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
      catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
      catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
      catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/3
      text(x = -0.1, y = catylocs, labels = catnames, pos = 2, srt = 90, col = "grey35")
    }
    
    #plot names of traits
    if(change_names_in_plot){
      coloc_phenotypes_newnames <- traitname_map[match(coloc_phenotypes, traitname_map[,1]),2]
    } else {
      coloc_phenotypes_newnames <- coloc_phenotypes
    }
    
    if(nicole_mods){
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = 1, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_newnames), ": ", coloc_phenotypes_newnames)[order_traits], col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    
    text(labels = "Phenotypes", x = nameloc, y = 10.15, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    if(use_geometric_mean){
      text(labels = "Ratio in Geometric Mean Posterior Probability of Colocalization (DE / nDE)", x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
    } else {
      text(labels = "Ratio in Posterior Probability of Colocalization in at Least One Gene", x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)  
    }
    
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5) 
    }
    
    #axes
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.2, y1 = -0.1, lwd = 2)
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
    
    #horiz axis ticks and nums
    if(use_geometric_mean){
      
      probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
      segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               y0 = -0.1, y1 = -0.15, lwd = 2)
      text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
           y = -0.135, pos = 1, cex = 1.25,
           labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
      
    } else {
      
      segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               y0 = -0.1, y1 = -0.15, lwd = 2)
      text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
           y = -0.135, pos = 1, cex = 1.25,
           labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
    }
    
    
    #let's get minor tick marks in there too
    if(use_geometric_mean){
      minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
      minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
      minticks <- minticks[minticks < max_horiz_axis]
      segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
      #hmm being weird
    } else {
      minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
      minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
      minticks <- minticks[minticks < max_horiz_axis]
      segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    }
    
    
    #dashed lines for ticks
    if(use_geometric_mean){
      if(nicole_mods){
        segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
                 x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
      } else {
        segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
                 x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
      }
    } else {
      if(nicole_mods){
        segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
                 x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
      } else {
        segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
                 x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = -0.1, lwd = 1, lty = 2, col = "grey75")
      }  
    }
    
    
    #helpful hint line about direction
    
    if(use_geometric_mean){
      segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
               x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
      text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
      text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
      text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    } else {
      segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
               x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
      text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
      text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
      text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    }
    
    #plot legend for colors etc.
    if(nicole_mods){rect(xleft = 1.97, ybottom = 8.35, ytop = 10.1, xright = 2.2, border = NA, col = "white")}
    # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
    #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    text(labels = sapply(names(deg_eqtl_list)[-18], function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
         y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list)[-18])), x = 2, pos = 4)
    points(x = rep(1.9925, length(names(deg_eqtl_list)[-18])), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list)[-18])), col = cols$Tissue, pch = 15, cex = 1.75)
    
    text(labels = c("Male", "Female"),
         y = seq(8.5, 8.425, length.out = 2), x = 2, pos = 4)
    points(x = rep(1.9925, 2), y = seq(8.505, 8.43, length.out = 2), col = cols$Sex, pch = c(18, 16), cex = 2)
    
    for(tissue in names(deg_eqtl_list)){
      
      for(sex in 1:2){
        
        sex_inds <- which(deg_eqtl_list[[tissue]]$sex == c("male", "female")[sex])
        
        if((tissue == "t64-ovaries" & sex == 1) | (tissue == "t63-testes" & sex == 2) | (tissue == "t1000-gonads")){
          next()
        }
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        
        if(use_geometric_mean){
          trait_probs <- pp4s_DE_notDE_array[sex, timepoint, tissue,,"DE_coloc_probs"] - 
            pp4s_DE_notDE_array[sex, timepoint, tissue,,"nDE_coloc_probs"]
        }
        
        trait_probs <- trait_probs[order_traits]
        
        points(x = (trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = trait_locs, pch = c(18, 16)[sex], col = grDevices::adjustcolor(cols$Tissue[tissue], alpha.f = 0.75), cex = 2)
        
        
        
      }
      
    }
    
    if(arnold_mods){
      addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
    }
    
    dev.off()
  }
  
  
  
}

#sanity check
pp4s_DE_notDE_array["female", "8w", tissues[grep(tissues, pattern = "lung")], coloc_phenotypes[grep(coloc_phenotypes, pattern = "height")], "DE_coloc_probs"] - 
  pp4s_DE_notDE_array["female", "8w", tissues[grep(tissues, pattern = "lung")], coloc_phenotypes[grep(coloc_phenotypes, pattern = "height")], "nDE_coloc_probs"]


#stephen requests genes that colocalize in male colon, male small intestine, female ovaries for alzheimers
tissues_to_focus_on <- c("t61-colon", "t64-ovaries", "t67-small-intestine")
DE_genes_in_tissues <- deg_eqtl_list[tissues_to_focus_on]
DE_genes_in_tissues <- lapply(names(DE_genes_in_tissues), function(x) DE_genes_in_tissues[[x]][(DE_genes_in_tissues[[x]]$selection_fdr < sign_filter_alpha & DE_genes_in_tissues[[x]]$adj_p_value < sign_filter_alpha),])
# DE_genes_in_tissues <- lapply(names(DE_genes_in_tissues), function(x) DE_genes_in_tissues[[x]][(DE_genes_in_tissues[[x]]$adj_p_value < sign_filter_alpha),])
names(DE_genes_in_tissues) <- tissues_to_focus_on
DE_genes_in_tissues <- list(
  DE_genes_in_tissues$`t61-colon`$gene_id[DE_genes_in_tissues$`t61-colon`$comparison_group == "8w" & DE_genes_in_tissues$`t61-colon`$sex == "male"],
  DE_genes_in_tissues$`t64-ovaries`$gene_id[DE_genes_in_tissues$`t64-ovaries`$comparison_group == "8w" & DE_genes_in_tissues$`t64-ovaries`$sex == "female"],
  DE_genes_in_tissues$`t67-small-intestine`$gene_id[DE_genes_in_tissues$`t67-small-intestine`$comparison_group == "8w" & DE_genes_in_tissues$`t67-small-intestine`$sex == "male"]
); names(DE_genes_in_tissues) <- tissues_to_focus_on

colocs_in_tissues <- colocs[[traitname_map$Tag[grep("Alzh", traitname_map$new_Phenotype)]]][tissues_to_focus_on]

colocs_in_tissues$`t61-colon`$gene_name <- deg_eqtl_list$`t61-colon`$gene_name[match(colocs_in_tissues$`t61-colon`$gene_id, deg_eqtl_list$`t61-colon`$gene_id)]
colocs_in_tissues$`t64-ovaries`$gene_name <- deg_eqtl_list$`t64-ovaries`$gene_name[match(colocs_in_tissues$`t64-ovaries`$gene_id, deg_eqtl_list$`t64-ovaries`$gene_id)]
colocs_in_tissues$`t67-small-intestine`$gene_name <- deg_eqtl_list$`t67-small-intestine`$gene_name[match(colocs_in_tissues$`t67-small-intestine`$gene_id, deg_eqtl_list$`t67-small-intestine`$gene_id)]

colocs_in_tissues$`t61-colon`[match(DE_genes_in_tissues$`t61-colon`, colocs_in_tissues$`t61-colon`$gene_id),]
colocs_in_tissues$`t64-ovaries`[match(DE_genes_in_tissues$`t64-ovaries`, colocs_in_tissues$`t64-ovaries`$gene_id),]
colocs_in_tissues$`t67-small-intestine`[match(DE_genes_in_tissues$`t67-small-intestine`, colocs_in_tissues$`t67-small-intestine`$gene_id),]

(mean(log(colocs_in_tissues$`t61-colon`[match(DE_genes_in_tissues$`t61-colon`, colocs_in_tissues$`t61-colon`$gene_id),"p4"])) -
  mean(log(colocs_in_tissues$`t61-colon`[-match(DE_genes_in_tissues$`t61-colon`, colocs_in_tissues$`t61-colon`$gene_id),"p4"])))
(mean(log(colocs_in_tissues$`t64-ovaries`[match(DE_genes_in_tissues$`t64-ovaries`, colocs_in_tissues$`t64-ovaries`$gene_id),"p4"])) - 
  mean(log(colocs_in_tissues$`t64-ovaries`[-match(DE_genes_in_tissues$`t64-ovaries`, colocs_in_tissues$`t64-ovaries`$gene_id),"p4"])))
(mean(log(colocs_in_tissues$`t67-small-intestine`[match(DE_genes_in_tissues$`t67-small-intestine`, colocs_in_tissues$`t67-small-intestine`$gene_id),"p4"])) - 
  mean(log(colocs_in_tissues$`t67-small-intestine`[-match(DE_genes_in_tissues$`t67-small-intestine`, colocs_in_tissues$`t67-small-intestine`$gene_id),"p4"])))

mean(mean(log(colocs_in_tissues$`t61-colon`[match(DE_genes_in_tissues$`t61-colon`, colocs_in_tissues$`t61-colon`$gene_id),"p4"])) > 
       log(colocs_in_tissues$`t61-colon`[-match(DE_genes_in_tissues$`t61-colon`, colocs_in_tissues$`t61-colon`$gene_id),"p4"]))
mean(mean(log(colocs_in_tissues$`t64-ovaries`[match(DE_genes_in_tissues$`t64-ovaries`, colocs_in_tissues$`t64-ovaries`$gene_id),"p4"])) > 
       log(colocs_in_tissues$`t64-ovaries`[-match(DE_genes_in_tissues$`t64-ovaries`, colocs_in_tissues$`t64-ovaries`$gene_id),"p4"]))
mean(mean(log(colocs_in_tissues$`t67-small-intestine`[match(DE_genes_in_tissues$`t67-small-intestine`, colocs_in_tissues$`t67-small-intestine`$gene_id),"p4"])) > 
       log(colocs_in_tissues$`t67-small-intestine`[-match(DE_genes_in_tissues$`t67-small-intestine`, colocs_in_tissues$`t67-small-intestine`$gene_id),"p4"]))



#### LDSC output ####

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


if(file.exists("~/data/smontgom/ldsc_cluster_results.txt")){
  ldsc_results <- as.data.frame(fread("~/data/smontgom/ldsc_cluster_results.txt"))
} else{
  ldsc_results <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results) <- colnames(fread(ldsc_results_paths[[1]][1]))
  ldsc_results$cluster <- rep(1:length(cluster_names), length(gwas_names))
  ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 1:nrow(ldsc_results)){
    output <- fread(ldsc_results_paths[[ldsc_results$gwas[i]]][grep(pattern = paste0(cluster_names[ldsc_results$cluster[i]], ".res"), ldsc_results_paths[[ldsc_results$gwas[i]]])][1])
    ldsc_results[i,1:ncol(output)] <- output[grep(pattern = cluster_names[ldsc_results$cluster[i]], output$Category),]
  }
  fwrite(ldsc_results, "~/data/smontgom/ldsc_cluster_results.txt")
}

ldsc_results$gwas <- stringr::str_replace_all(ldsc_results$gwas, "imputed_", "")

ldsc_results$logPVal_enrichment <- log10(ldsc_results$Enrichment_p)
hist(ldsc_results$logPVal_enrichment)
hist(ldsc_results$logPVal_enrichment + log10(nrow(ldsc_results)))
abline(v = log10(0.1), col = 2)
plot((ldsc_results$Prop._h2 - ldsc_results$Prop._SNPs) / ldsc_results$Prop._h2_std_error, ldsc_results$Enrichment_p)
plot(ldsc_results$Enrichment / ldsc_results$Enrichment_std_error, ldsc_results$Enrichment_p)
plot(ldsc_results$Prop._h2 / ldsc_results$Prop._SNPs, ldsc_results$Enrichment)

cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"
cols$cluster <- cols$Tissue[1:length(cluster_names)]
names(cols$cluster) <- cluster_names

#read in conditional analyses
if(file.exists("~/data/smontgom/ldsc_cluster_results_conditional.txt")){
  ldsc_results_conditional <- as.data.frame(fread("~/data/smontgom/ldsc_cluster_results_conditional.txt"))
} else{
  ldsc_results_conditional_paths <- list.files(ldsc_output_dir)
  ldsc_results_conditional_paths <- paste0(ldsc_output_dir, ldsc_results_conditional_paths[grep(ldsc_results_conditional_paths, pattern = "results")])
  cluster_names <- paste0("Cluster_", 1:15)
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) ldsc_results_conditional_paths[grep(pattern = gwas, ldsc_results_conditional_paths)])
  names(ldsc_results_conditional_paths) <-  gwas_names
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) 
    ldsc_results_conditional_paths[[gwas]][grep(pattern = "1-15", ldsc_results_conditional_paths[[gwas]])])
  names(ldsc_results_conditional_paths) <-  gwas_names
  
  ldsc_results_conditional <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_conditional_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_conditional) <- colnames(fread(ldsc_results_conditional_paths[[1]][1]))
  ldsc_results_conditional$cluster <- rep(1:length(cluster_names), length(gwas_names))
  ldsc_results_conditional$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    output <- fread(ldsc_results_conditional_paths[[ldsc_results_conditional$gwas[i]]][1])
    output$Category <- stringr::str_replace(pattern = "L2_0", replacement = "", output$Category)
    output <- output[output$Category %in% cluster_names,]
    reorder_output <- order(match(stringr::str_replace(output$Category, pattern = "Cluster_", replacement = ""), 
          ldsc_results_conditional[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_conditional[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_conditional$gwas <- stringr::str_remove(ldsc_results_conditional$gwas, "imputed_")
  ldsc_results_conditional$logPVal_enrichment <- log10(ldsc_results_conditional$Enrichment_p)
  fwrite(ldsc_results_conditional, "~/data/smontgom/ldsc_cluster_results_conditional.txt")
}

all(as.numeric(stringr::str_remove_all(ldsc_results_conditional$Category, pattern = "Cluster_")) == ldsc_results_conditional$cluster)





coloc_phenotypes_sub <- traits_with_satisfactory_heritaility


#specify graph parameters
max_point_cex <- 3.5
point_cex_power <- 0.5
n_points_for_legend = 7
buffer_min_and_max = 0.05
minimum_enrichment_logPval <- min(ldsc_results_sub$logPVal_enrichment)
opacity_insig_points <- 0.2
opacity_sig_points <- 0.8

plot_LDSC_comparison = T
reorder_vertical <- T
use_conditional_model = F
use_enrichment = F
partition_by_category <- T
total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])
ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_with_satisfactory_heritaility,]

if(use_conditional_model){
  ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_with_satisfactory_heritaility,]
}
use_heritability = !use_enrichment
fix_axes = F
fix_axes_bounds_enrichment = c(0,5)
fix_axes_bounds_heritability = c(0,1)

total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])
ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_with_satisfactory_heritaility,]

if(use_conditional_model){
  ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_with_satisfactory_heritaility,]
}

trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))

bad_boys <- sort(unique(ldsc_results_sub$gwas[ldsc_results_sub$Enrichment < 0 | ldsc_results_sub$Prop._h2 > 1 | ldsc_results_sub$Prop._h2 < 0]))
bad_boys_cols <- rep(1, length(gwas_names))
bad_boys_cols[match(bad_boys, coloc_phenotypes_sub)] <- 2

plot(density(log_files$h2))
segments(x0 = c(log_files$h2[log_files$gwas %in% bad_boys]), x1 = c(log_files$h2[log_files$gwas %in% bad_boys]), y0 = 0, y1 = 5, col = "red", lwd = 0.4)
hist(log_files$h2[log_files$gwas %in% bad_boys])
hist(sapply(log_files$h2[log_files$gwas %in% bad_boys], function(h) mean(h > log_files$h2)), cex.main = 1)
par(mar = c(4,4,4,4))
plot(log_files$h2[log_files$gwas %in% bad_boys], log_files$h2se[log_files$gwas %in% bad_boys], xlim = c(-0.003, 0.006), ylim = c(-0.003, 0.006), 
     xlab = "total heritability", ylab = "total heritability SE")
hist(log_files$h2[log_files$gwas %in% bad_boys] / log_files$h2se[log_files$gwas %in% bad_boys])
hist(log_files$h2 / log_files$h2se, breaks = seq(min(log_files$h2 / log_files$h2se), max(log_files$h2 / log_files$h2se), length.out = 100))
abline(a = 0, b = 1, xpd = F)
hist(log_files$h2, breaks = seq(min(log_files$h2), max(log_files$h2), length.out = 50))
hist(log_files$chi2, breaks = seq(min(log_files$chi2), max(log_files$chi2), length.out = 50))
hist(log_files$chi2[log_files$gwas %in% bad_boys])
log_files$gwas[log_files$h2 > 1]

change_names <- F
change_names_in_plot = T
nicole_mods = T
arnold_mods = F
if(plot_LDSC_comparison){
  
  max_horiz_axis <- 2.05
  
    
    if(reorder_vertical){
      
      if(use_enrichment){
        order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Enrichment)), decreasing = T)
      } else {
        order_traits <- order(sapply(coloc_phenotypes_sub, function(coloc_phenotype) mean(ldsc_results_sub[ldsc_results_sub$gwas == coloc_phenotype,]$Prop._h2)), decreasing = T)
      }
      
    } else {
      order_traits <- 1:length(coloc_phenotypes_sub)
    }
    
    if(partition_by_category){
      order_traits <- order_traits[order(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])], sort(order_traits))]
      trait_category_counts <- as.data.table(table(traitwise_partitions[,2][match(coloc_phenotypes_sub[order_traits], traitwise_partitions[,1])] ))
      trait_category_counts$cN <- cumsum(trait_category_counts$N)
    }
    
    if(use_enrichment){
      if(fix_axes){
        logprob_range <- fix_axes_bounds_enrichment
      } else {
        logprob_range <- range(ldsc_results_sub$Enrichment)
      }
    } else {
      if(fix_axes){
        logprob_range <- fix_axes_bounds_heritability
      } else {
        logprob_range <- range(ldsc_results_sub$Prop._h2)
      }
    }
    
    logprob_range <- c(logprob_range[1] - diff(logprob_range) * buffer_min_and_max / 2, logprob_range[2] + diff(logprob_range) * buffer_min_and_max / 2)
    
    max_prob <- diff(logprob_range)
    
    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/LDSC_comparison_heritabilities", ifelse(use_enrichment, "_enrichment", "_heritability"), 
                                           ifelse(fix_axes, "_fixedAxes", ""), ifelse(use_conditional_model, "_conditionalModel", "_unconditionalModel"), ".pdf"), 
                         width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub) , family="Arial Unicode MS")
    par(mfrow = c(1, 1), mar = c(4,3,3,3), xpd = NA)
    
    
    plot(1,1,xlim = c(0,2.1), ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
    
    ylocs <- seq(10, 0, length.out = length(coloc_phenotypes_sub))
    
    #put boxes around categories, if we're doing that  
    if(partition_by_category){
      segments(x0 = -0.1, x1 = -0.1, y0 = 0+diff(ylocs)[1]/2, y1 = 10-diff(ylocs)[1]/2, col = "grey35")
      for(traitcat in 1:(nrow(trait_category_counts)-1)){
        segments(y0 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, y1 = ylocs[trait_category_counts$cN[traitcat]] + diff(ylocs)[1]/2, x0 = -0.1, x1 = nameloc, col = "grey35")
      }
      
      catylocs <- (ylocs[trait_category_counts$cN] + c(ylocs[1], ylocs[trait_category_counts$cN][-length(trait_category_counts$cN)])) / 2 - diff(ylocs)[1]/1.5
      catnames <- as.vector(trait_category_counts$V1)
      catnames <- sapply(catnames, function(x) strsplit(x, split = c("-"))[[1]][1])
      catnames <- sapply(catnames, function(x) strsplit(x, split = c(" "))[[1]][1])
      catylocs[catnames == "Hair"] <- catylocs[catnames == "Hair"] + diff(ylocs)[1]/1.5
      catylocs[catnames == "Aging"] <- catylocs[catnames == "Aging"] - diff(ylocs)[1]/3
      catylocs[catnames == "Endocrine"] <- catylocs[catnames == "Endocrine"] - diff(ylocs)[1]/3
      for(catname in 1:length(catnames)){
        text(x = -0.1, y = catylocs[catname] + ifelse(catnames[catname] == "Endocrine", 0.1, 0), labels = catnames[catname], pos = 2, srt = 90, col = "grey35")
      }
      
    }
    
    #plot names of traits
    if(change_names_in_plot){
      coloc_phenotypes_sub_newnames <- traitname_map[match(coloc_phenotypes_sub, traitname_map[,1]),2]
    } else {
      coloc_phenotypes_sub_newnames <- coloc_phenotypes_sub
    }
    
    if(nicole_mods){
      text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", coloc_phenotypes_sub_newnames)[order_traits], col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", coloc_phenotypes_sub_newnames)[order_traits], col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    #vert axis label
    text(labels = "Phenotypes", x = nameloc, y = 10.15, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    #horiz axis label
    text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Clusters"), x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
    
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
    }
    
    #axes
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.2, y1 = -0.1, lwd = 2)
    segments(x0 = nameloc, x1 = max_horiz_axis, y0 = -0.1, y1 = -0.1, lwd = 2)
    
    #horiz axis ticks and nums
    # if(use_geometric_mean){
    #   
    #   probs <- ceiling(logprob_range[1]):ceiling(logprob_range[2])
    #   segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            y0 = -0.1, y1 = -0.15, lwd = 2)
    #   text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #        y = -0.135, pos = 1, cex = 1.25,
    #        labels =  c(sapply(probs, function(num) as.expression(bquote(10^.(num))))))
    #   
    # } else {
    #   
    #   segments(x0 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            x1 = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            y0 = -0.1, y1 = -0.15, lwd = 2)
    #   text(x = seq(0, max_prob, by = 1) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #        y = -0.135, pos = 1, cex = 1.25,
    #        labels =  c(sapply(seq(ceiling(logprob_range[1]), ceiling(logprob_range[2]), by = 1), function(num) as.expression(bquote(10^.(num))))))
    # }
    
    if(use_enrichment){
      probs <- seq(ceiling(logprob_range[1]), floor(logprob_range[2]), by = ifelse(fix_axes, 0.5, 1))  
    } else {
      probs <- seq(ceiling(logprob_range[1]), (logprob_range[2]), by = ifelse(fix_axes, 0.1, 0.1))
    }
    
    segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
             x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
             y0 = -0.1, y1 = -0.15, lwd = 2)
    text(x = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc,
         y = -0.135, pos = 1, cex = 1.25,
         labels =  probs)
  
    
    # #let's get minor tick marks in there too
    # if(use_geometric_mean){
    #   minticks <- log10(as.vector(sapply(1:(length(probs)-1), function(p) seq(10^probs[p], 10^probs[p+1], length.out = 10))))
    #   minticks <- (minticks - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc
    #   minticks <- minticks[minticks < max_horiz_axis]
    #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    #   #hmm being weird
    # } else {
    #   minticks <- log10(as.vector(t(t(1:9)) %*% t(10^seq(0, max_prob, by = 10^round(log(max_prob))))))
    #   minticks <- minticks / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025
    #   minticks <- minticks[minticks < max_horiz_axis]
    #   segments(x0 = minticks, x1 =  minticks, y0 = -0.1, y1 = -0.125, lwd = 1)
    # }
    
    
    #dashed lines for ticks
    segments(x0 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             x1 = (probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = -0.075, lwd = 2, lty = 2, col = "grey75")
    
    
    #helpful hint line about direction
    
    # if(use_geometric_mean){
    #   segments(x0 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
    #            x1 = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
    #   text(labels = "non-DE Genes", srt = 90, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
    #   text(labels = "DE Genes", srt = 270, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
    #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = (0 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    # } else {
    #   segments(x0 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, 
    #            x1 = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y0 = 10, y1 = 10.4, lwd = 2, lty = 2, col = "grey75")
    #   text(labels = "eQTL > DEG", srt = 90, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 2, col = "grey75", cex = 0.75)
    #   text(labels = "DEG > eQTL", srt = 270, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.4, pos = 4, col = "grey75", cex = 0.75)
    #   text(labels = "Stronger Coloc\nSignal In", srt = 0, x = floor(-logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc + 0.025, y = 10.3875, pos = 3, col = "grey75", cex = 0.75)
    # }
    
    #plot legend for colors etc.
    rect(xleft = 1.97, ybottom = 8.15 - 0.05*n_points_for_legend, ytop = 10.1, xright = 2.2, border = NA, col = "white")
    # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
    #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    text(labels = stringr::str_replace_all(cluster_names, "_", " "),
         y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
    points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
    
    #plot legend for points
    pt_loc_expander <- 2
    point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
    points_legend_pchs <- rep(19, n_points_for_legend)
    point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
    points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    
    points(y = 8.55 - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
           x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
    text(labels = latex2exp::TeX("log_{10}(enrichment p-val)"), y = 8.6 , x = 1.975, cex = 1.1, pos = 4,  font = 2)
    text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
    text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
         y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
    points(pch = c(19,18), y = c(8.325, 8.225) , x = rep(2.175,2), cex = c(1.25, 1.75), 
           col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
    
    
    text(labels = point_legend_cexes_log10_pvals,y = 8.545 - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
           x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
    
    
    #figure out cex params
    for(cluster in 1:length(cluster_names)){
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        if(use_enrichment){
          trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
          trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
        } else {
          trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
          trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
        }
        trait_probs <- trait_probs[order_traits]
        
        point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
        point_cex <- point_cex * max_point_cex
        
        good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
        
        pchs <- rep(19, length(trait_probs))
        pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
        opacities <- rep(opacity_insig_points, length(trait_probs))
        opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
        point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[cluster], alpha.f = opcty))
        
        points(x = ((trait_probs - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc)[good_points], y = (trait_locs)[good_points], 
               pch = pchs[good_points], col = point_cols[good_points], cex = point_cex[good_points])
        
        
      }
      
    
    
    if(arnold_mods){
      addImg(obj = arnold_bayes, x = 1.9, y = 1.6775, width = 0.5)
    }
    
    dev.off()
}
  
#quick comparison of conditional vs unconditional models
par(mfrow = c(1,2), mar = c(4,4,4,2))
plot(ldsc_results_conditional$Enrichment,  ldsc_results$Enrichment, main = "Heritability Enrichment",
     xlab = "Conditional Model Enrichments", ylab = "Unconditional Model Enrichments", pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)
plot(ldsc_results_conditional$Prop._h2,  ldsc_results$Prop._h2, main = "Proportion Heritability", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))
abline(0,1, col = "red", lty = 2)

#### cluster-based marginal enrichment ####

trait_type_specific = T
grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/LDSC_cluster-specific_enrichments", ifelse(trait_type_specific, "_categorized-by-categories"),".pdf"), 
                     width = 1500 / 72, height = 2000 / 72, family="Arial Unicode MS")
par(mfrow = c(1, 1), mar = c(4,3,1,3), xpd = NA)

range_xs <- c(0, max(ldsc_results_sub$Enrichment[!(ldsc_results_sub$gwas %in% bad_boys)]) * 1.2)
plot(1,1,xlim = range_xs, ylim = c(0,10), col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
cis <- 1:15

ylocs_bottoms <- seq(0, 10, length.out = length(cis) + 1)[-(length(cis) + 1)]
yheight = diff(ylocs_bottoms)[1] * 0.9

order_cs <- order(sapply(cis, function(ci) mean(ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == ci & !(ldsc_results_sub$gwas %in% bad_boys)])))
unique_trait_categories <- unique(traitwise_partitions$Category)
cols$category <- cols$Tissue[1:length(unique_trait_categories)+1]
names(cols$category) <- unique_trait_categories
for(ci in cis){
  if(trait_type_specific){
    for(ti in unique_trait_categories){
      
      enrichments <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == ci & 
                                                !(ldsc_results_sub$gwas %in% bad_boys) & 
                                                  ldsc_results_sub$gwas %in% traitwise_partitions$Tag[traitwise_partitions$Category == ti]]
      if(length(enrichments) == 0){
        next()
      } else if(length(enrichments) == 1){
        dens_enrich <- density(enrichments, from = range_xs[1], to = range_xs[2], bw = 0.1)
      } else {
        dens_enrich <- density(enrichments, from = range_xs[1], to = range_xs[2])
      }
      ylocs <- dens_enrich$y / max(dens_enrich$y) * yheight + ylocs_bottoms[order_cs[ci]]
      xlocs <- dens_enrich$x
      ylocs[1] <- ylocs[length(ylocs)] <- ylocs_bottoms[order_cs[ci]]
      polygon(x = xlocs, y = ylocs, col = adjustcolor(cols$category[ti], 0.5))
      text(labels = paste0("Cluster ", ci), x = range_xs[1], y = ylocs_bottoms[order_cs[ci]] + 0.55, pos = 2, col = "black", cex = 1.5, srt = 90)
      }
    
  } else {
    enrichments <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == ci & !(ldsc_results_sub$gwas %in% bad_boys)]
    dens_enrich <- density(enrichments, from = range_xs[1], to = range_xs[2])
    ylocs <- dens_enrich$y / max(dens_enrich$y) * yheight + ylocs_bottoms[order_cs[ci]]
    xlocs <- dens_enrich$x
    ylocs[1] <- ylocs[length(ylocs)] <- ylocs_bottoms[order_cs[ci]]
    polygon(x = xlocs, y = ylocs, col = "lightgrey")
    text(labels = paste0("Cluster ", ci), x = xlocs[min(which(cumsum(dens_enrich$y) / sum(dens_enrich$y) > 0.5))], y = ylocs_bottoms[order_cs[ci]], pos = 3, col = "white", cex = 2)  
  }
  
}
segments(x0 = range_xs[1], x1 = range_xs[2], lwd = 3, y0 = -0.1)
xticks <- seq(range_xs[1], range_xs[2], by = 1)
segments(x0 = xticks, x1 = xticks, y0 = -0.1, y1 = -0.15, lwd = 3)
text(labels = xticks, x = xticks, y = -0.15, pos = 1, cex = 2)
text(labels = "Enrichment Proportion", x = diff(range_xs) / 2, y = -0.35, pos = 1, cex = 3)
text(labels = "LDSC Enrichment Across Clusters", x = diff(range_xs) / 2, y = 10.1, pos = 3, cex = 4)

if(trait_type_specific){
  legend(legend = names(cols$category), x = range_xs[2] - diff(range_xs) / 5, y = 10.45, col = adjustcolor(cols$category, 0.75), pch = 15, pt.cex = 3, ncol = 2, cex = 1.1)
}

dev.off()

#### pairwise genetic correlation matrix ####
ldsc_directory <- "~/repos/ldsc/"

#could maybe have genetic correlations in the upper right, significance in the bottom left

#let's get the gwas summary stats into a format that munge_sumstats.py can process
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
traitnames <- gwas_summary_files
traitnames <- sapply(traitnames, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
traitnames <- sapply(traitnames, function(trait) strsplit(trait, "imputed_")[[1]][2])
traitnames <- as.character(traitnames)
gcors <- lapply(traitnames, function(x) integer())
names(gcors) <- traitnames  
for(gi in 1:length(gwas_summary_files)){
 logfile <- readLines(paste0(ldsc_directory, "output/signed_gcor/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_pairwise-Gcorrs.log"))
 logfile <- logfile[(grep(x = logfile, pattern = "Summary of Genetic Correlation Results")+1):(length(logfile)-3)]
 writeLines(logfile, con = paste0(ldsc_directory, "temp.txt"))
 gcors[[gi]] <- fread(paste0(ldsc_directory, "temp.txt"))
 file.remove(paste0(ldsc_directory, "temp.txt"))
 
 gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
 gcors[[gi]]$p1 <- sapply(gcors[[gi]]$p1, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
 gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
 gcors[[gi]]$p2 <- sapply(gcors[[gi]]$p2, function(trait) strsplit(trait, "gwas_sumstats/proper_format/signed/imputed_")[[1]][2])
}

gcor_mat <- diag(length(traitnames))
colnames(gcor_mat) <- rownames(gcor_mat) <- traitnames
for(ri in rownames(gcor_mat)){
  for(ci in colnames(gcor_mat)){
    gcor_mat[ri, ci] <- gcors[[ri]]$rg[match(ci, gcors[[ri]]$p2)]  
  }
}
diag(gcor_mat) <- rep(1, length(traitnames))

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

#### actually plot genetic correlations ####

grDevices::cairo_pdf(filename = paste0("~/Documents/DExEQTL/genetic_correlations.pdf"), 
                     width = 1200 / 72, height = 600 / 72, family="Arial Unicode MS")
par(mar = c(8,6,3,4), mfrow = c(1,2), xpd = NA)


#plotting params
traits <- rownames(gcor_mat)
rate = 0.001
exp_dist_cols <- round(cumsum(c(1, dexp(1:100, rate = rate) / min(dexp(1:100, rate = rate)))))
heatcols <- viridis::magma(max(exp_dist_cols))[exp_dist_cols]
heatcols <- RColorBrewer::brewer.pal(11, "RdBu")[-c(1,11)]
heatcols <- rev(colorRampPalette(heatcols)(101))
# plot(1:length(heatcols), 1:length(heatcols), cex = 2, pch = 19, col = heatcols)

par(mar = c(8,6,3,4), xpd = NA)
plot(1,1,xlim = c(0,length(traits)), ylim = c(0,length(traits)), xpd = NA,
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
sorted_row_inds <- order(cmdscale(1-gcor_mat, k = 1))
for(rowi in 1:length(traits)){
  text(labels = traits[sorted_row_inds[rowi]], x = length(traits), y = length(traits) - rowi + 1, 
       col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
         match(traits[sorted_row_inds[rowi]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
       , pos = 4, xpd = NA, cex = 0.5, font = 1)
  for(colj in 1:length(traits)){
    rect(yb = length(traits) - rowi + 1 - 1/2, yt = length(traits) - rowi + 1 + 1/2, 
         xl = colj - 1/2, xr = colj + 1/2, pch = 15, cex = 1, border = NA,
           col = heatcols[round((gcor_mat[sorted_row_inds[rowi], sorted_row_inds[colj]] + 1) / 2 * 100) + 1])
    
    if(gcor_mat_sig[sorted_row_inds[rowi], sorted_row_inds[colj]]){
      points(y = length(traits) - rowi + 1, x = colj, pch = 19, col = "white", cex = 0.2)
    }
    
    if(rowi == 1){
      text(labels = traits[sorted_row_inds[colj]], x = colj-1.25, y = 0, 
           col = cols$category[traitwise_partitions$Category[match(traitname_map$Tag[
             match(traits[sorted_row_inds[colj]], traitname_map$new_Phenotype)],traitwise_partitions$Tag)]]
           , pos = 4, srt = 270, xpd = NA, cex = 0.5, font = 1)
    }
    
    
    
    
    
  }
}

xl = -3; xr = -1; yb = 27; yt = 57
rect(xleft = xl, xright = xr, col = heatcols, border = NA,
     yb = seq(yb, yt, length.out = length(heatcols)+1)[-(length(heatcols)+1)], yt= seq(yb, yt, length.out = length(heatcols)+1)[-1])
rect(xleft = xl, xright = xr, ybottom = yb, ytop = yt)
text(labels = -5:5/5, x = xl+1, pos = 2, y = seq(yb, yt, length.out = 11), cex = 1)
text(x = mean(c(xl, xr)), y = yt - 1, labels = latex2exp::TeX("r_g"), pos = 3, cex = 2, font = 2)
text("Genetic Correlation Matrix", x = length(traits) / 2 + 0.5, y = length(traits) + 0.5, pos = 3, cex = 2.5, font = 2)

points(rep((xl+xr)/2, length(cols$category)), 1:length(cols$category)*2, pch = 15, col = cols$category, cex = 1.75)
text(rep((xl+xr)/2, length(cols$category)), 1:length(cols$category)*2, pos = 2, col = 1, 
     labels = sapply(sapply(names(cols$category), function(x) strsplit(x, " ")[[1]][1]), function(y) strsplit(y, "-")[[1]][1]), cex = 0.75)

text(labels = latex2exp::TeX(paste0("• mark p-val < 10^{", round(log10(0.025  / choose(dim(gcor_mat)[1], 2)), 2), "}")), 
     x = -1, y = yt + nrow(gcor_mat)/5, cex = 0.75, srt = 90)
text(labels = "bonferroni", x = -2.5, y = yt + nrow(gcor_mat)/5, cex = 0.75, srt = 90, font = 3, family = "Arial")

eigcor <- eigen(Matrix::nearPD(gcor_mat)$mat)
par(mar = c(8,7,3,4), xpd = NA)
dcor <- dim(gcor_mat)[1]
cols_scree <- c("darkblue", "darkorange")
plot(eigcor$values / sum(eigcor$values) * 100 / max(eigcor$values / sum(eigcor$values) * 100), type = "l", lwd = 2, 
     ylab = "", xlab = "", cex.axis = 1.25, cex.lab = 1.5, xaxt = "n", yaxt = "n", frame = F, col = cols_scree[1], ylim = c(0,1))
text("Proportion Variance Explained", y = 0.5, srt = 90, pos = 3, x = -dcor/10, xpd = NA, cex = 2, col = cols_scree[1])
text("Cum. Proportion Variance Explained", y = 0.5, srt = 270, pos = 3, x = (dcor)*1.1, xpd = NA, cex = 2, col = cols_scree[2])

lines(1:length(eigcor$values), cumsum(eigcor$values) / sum(eigcor$values) * 100 / max(cumsum(eigcor$values) / sum(eigcor$values) * 100), lwd = 2,
      col = cols_scree[2])
rect(xleft = 0, xright = dcor, ybottom = 0, ytop = 1, lwd =2)
segments(x0 = 0, x1 = -1, y0 = seq(0,1,0.1), y1 = seq(0,1,0.1), lwd= 2, col = cols_scree[1])
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, lwd= 2, col = cols_scree[1])
text(x = -1, y = seq(0,1,0.1), pos = 2, col = cols_scree[1], 
     labels = round(seq(0, max(eigcor$values / sum(eigcor$values) * 100), length.out= 11), 1))
segments(x0 = dcor, x1 = dcor+1, y0 = seq(0,1,0.1), y1 = seq(0,1,0.1), lwd= 2, col = cols_scree[2])
segments(x0 = dcor, x1 = dcor, y0 = 0, y1 = 1, lwd= 2, col = cols_scree[2])
text(x = dcor + 1, y = seq(0,1,0.1), pos = 4, col = cols_scree[2], 
     labels = round(seq(0, max(cumsum(eigcor$values) / sum(eigcor$values) * 100), length.out= 11), 1))
text(labels = "Scree Plot", cex = 3, x = dcor/2, y = 1, pos = 3)
segments(x0 = seq(1,dcor,10), x1 = seq(1,dcor,10), y0 = 0, y1 = -0.02, lwd= 2, col = 1)
text(x = seq(1,dcor,10), y = -0.02, lwd= 2, col = 1, pos = 1, labels = seq(1,dcor,10))
text(labels = "Index", cex = 2, x = dcor/2, y = -0.075, pos = 1)



dev.off()

# 
# #plot eigenvector loadings and scree plot
# wiou_eigen <- eigen(gcor_mat)
# 
# par(mar = c(3.5,0.5,0.5,0.5), xpd = NA)
# for(i in 1:length(traits)){
#   plot(100,100,xlim = c(0,length(traits)), ylim = c(0,max(abs(wiou_eigen$vectors[,i]))), 
#        col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
#   text(paste0("PC ", i, " Loadings"), x = length(traits) / 2, y = 0.95*max(abs(wiou_eigen$vectors[,i])), cex = 1.25)
#   rect(xleft = 1:length(traits)-1/2, xright = 1:length(traits)+1/2, 
#        ybottom = 0, ytop = abs(wiou_eigen$vectors[,i])[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], 
#        col = cols$Tissue[order(abs(wiou_eigen$vectors[,i]), decreasing = T)])
#   text(labels = traits[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], x = 1:length(traits) - 0.75, 
#        y = -max(abs(wiou_eigen$vectors[,i]))*0.01, col = cols$Tissue[order(abs(wiou_eigen$vectors[,i]), decreasing = T)], 
#        pos = 4, srt = 270)
# }
# 
# par(mar = c(5,3,1,2))
# plot(cex.lab = 1.5, x = 1:length(traits), y = wiou_eigen$values, type = "l", xlab = "ordered eigenvalue indices", ylab = "eigenvalue", lwd = 2)
# dev.off()


# for(i in 1:nrow(log_files)){
#   logfile <- readLines(ldsc_log_paths[[log_files$gwas[i]]][grep(pattern = paste0(cluster_names[log_files$cluster[i]], ".log"), ldsc_log_paths[[log_files$gwas[i]]])][1])
#   h2 <- logfile[grep(logfile, pattern = "Total Observed scale h2")]
#   log_files$h2[i] <- as.numeric(strsplit(h2, " ")[[1]][5])
#   log_files$h2se[i] <- as.numeric(substr(x = strsplit(h2, " ")[[1]][6], start = 2, stop = nchar(strsplit(h2, " ")[[1]][6]) - 1))
#   chi2 <- as.numeric(strsplit(logfile[grep(logfile, pattern = "Mean Chi")], split = " ")[[1]][3])
#   log_files$chi2[i] <- chi2
# }
# 


#### ldsc cts results ####
ldsc_results_dir <- "~/repos/ldsc/output/cts/"
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
orig_clusters <- paste0("Cluster_", 1:15)
new_clusters <- c("t53-cortex",  "t56-vastus-lateralis",  "t58-heart",  "t61-colon",  "t62-spleen",  "t64-ovaries",  "t65-aorta",  "t66-lung",  "t67-small-intestine",  "t69-brown-adipose",  "t70-white-adipose",  "t30-blood-rna",  "t52-hippocampus",  "t54-hypothalamus",  "t63-testes",  "t55-gastrocnemius",  "t60-adrenal",  "t59-kidney",  "t68-liver",  "1w",  "2w",  "4w",  "8w",  "female-1w",  "male-2w",  "female-4w",  "male-8w",  "all_DE")

geneset_info <- expand.grid(c(7,15), c("t55-gastrocnemius", "t68-liver", "t58-heart"))
new_clusters <- paste0("Cluster_", trimws(apply(geneset_info, 1, paste0, collapse = "-")), "")

output_files <- c(sapply(c("_Cluster_Addenda.cell_type_results.txt", "_Cluster_1-15.cell_type_results.txt"), function(x) paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), x)))
output_files <- paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), c("_Cluster_1-15.cell_type_results.txt"))
output_files <- paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), c("_Cluster_Addenda.cell_type_results.txt"))
ldsc_results <- do.call(rbind, lapply(1:length(output_files), function(i) cbind(fread(paste0(ldsc_results_dir, output_files[i])),  trait = strsplit(output_files[i], "_Cluster")[[1]][1])))

ldsc_results$twoTailedPVal <- ldsc_results$Coefficient_P_value
ldsc_results$twoTailedPVal[ldsc_results$twoTailedPVal > 0.5] <- 1 - ldsc_results$twoTailedPVal[ldsc_results$twoTailedPVal > 0.5]
ldsc_results$twoTailedPVal <- ldsc_results$twoTailedPVal * 2

#filter by category
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
ldsc_results$trait_name <- gsub(ldsc_results$trait, pattern = "imputed_", replacement = "")
ldsc_results$trait_category <- traitwise_partitions$Category[match(ldsc_results$trait_name, traitwise_partitions$Tag)]
ldsc_results <- ldsc_results[ldsc_results$trait_category %in% c("Cardiometabolic", "Aging", "Anthropometric", "Immune")]

hist(log(ldsc_results$Coefficient_P_value), xlab = "log 1-tailed p-value", main = "", breaks = -200:0/10)
abline(v = log(0.05 / nrow(ldsc_results)), col = 2, lwd = 4)
text(x = log(0.05 / nrow(ldsc_results)), y = 200, 
     labels = latex2exp::TeX("Bonferroni-adjusted $\\alpha$"), srt = 90, pos = 2)
hist(ldsc_results$Coefficient_P_value)
hist(log(ldsc_results$twoTailedPVal), breaks = -200:0/10)
hist((ldsc_results$twoTailedPVal))
abline(v = log(0.05 / nrow(ldsc_results)), col = 2, lwd = 4)
qvalue::pi0est(ldsc_results$Coefficient_P_value)
qvalue::pi0est(ldsc_results$twoTailedPVal)
sum(log(ldsc_results$Coefficient_P_value) < log(0.05 / nrow(ldsc_results)))
head(ldsc_results[order(ldsc_results$Coefficient_P_value, decreasing = F),], 20)
head(ldsc_results[order(ldsc_results$twoTailedPVal, decreasing = F),], 20)
head(ldsc_results[order(ldsc_results$Coefficient, decreasing = T),], 20)

#### TWAS results ####
twas_tissues <- fread(file = "~/repos/fusion_twas-master/output/all_results.txt")
metaxscan_results <- fread(file = "~/repos/MetaXcan/software/results/all_results.txt")
motrpac_gtex_map = c('t30-blood-rna'='Whole_Blood',
                     't52-hippocampus'='Brain_Hippocampus',
                     't53-cortex'='Brain_Cortex',
                     't54-hypothalamus'='Brain_Hypothalamus',
                     't55-gastrocnemius'='Muscle_Skeletal',
                     't56-vastus-lateralis'='Muscle_Skeletal',
                     't58-heart'='Heart_Left_Ventricle',
                     't59-kidney'='Kidney_Cortex',
                     't60-adrenal'='Adrenal_Gland',
                     't61-colon'='Colon_Transverse',
                     't62-spleen'='Spleen',
                     't63-testes'='Testis',
                     't64-ovaries'='Ovary',
                     't66-lung'='Lung',
                     't67-small-intestine'='Small_Intestine_Terminal_Ileum',
                     't68-liver'='Liver',
                     't70-white-adipose'='Adipose_Subcutaneous')

#compare twas to coloc results
load('~/data/smontgom/coloc_list_pp4threshold_1e-04.RData')
geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")

#whoops there is redundancy here
tiss_gene_twas <- paste0(twas_tissues$tissue, "_", twas_tissues$ID) 
tiss_gene_coloc <- paste0(coloc_list$motrpac_tissue, "_", coloc_list$gene_id) 

for(i in 1:length(deg_eqtl_list)){
  cat(paste0(i, " "))
  temp <- deg_eqtl_list[[i]]
  tiss_gene_twas_equiv <- paste0(motrpac_gtex_map[match(names(deg_eqtl_list)[i], names(motrpac_gtex_map))], "_", temp$gene_name)
  tiss_gene_coloc_equiv <- paste0(names(deg_eqtl_list)[i], "_", temp$gene_id)
  temp$TWAS_PVAL <- twas_tissues$TWAS.P[match(tiss_gene_twas_equiv, tiss_gene_twas)]
  temp$colocPP4 <- coloc_list$p4[match(tiss_gene_coloc_equiv, tiss_gene_coloc)]
  deg_eqtl_list[[i]] <- temp
}

twas_pvals <- (do.call(rbind,deg_eqtl_list)$TWAS_PVAL)
coloc_pp4s <- (do.call(rbind,deg_eqtl_list)$colocPP4)
plot(-log(twas_pvals), -log(1-coloc_pp4s))

all_gene_names <- unique(map$human_gene_symbol)
gene_names_DE <- unique(do.call(rbind, deg_eqtl_list)$gene_name)
mean(!is.na(match(gene_names_DE, all_gene_names)))
gene_names_TWAS <- unique(twas_tissues$ID)
gene_names_metaxscan <- unique(metaxscan_results$gene_name)
mean(!is.na(match(gene_names_TWAS, all_gene_names)))
length(setdiff(gene_names_DE, gene_names_TWAS)) / length(gene_names_DE)
sapply(names(deg_eqtl_list), function(i) mean(!is.na(deg_eqtl_list[[i]]$TWAS_PVAL))*100)
twas_tissuegenes <- sapply(unique(twas_tissues$tissue), function(tiss) unique(twas_tissues$ID[twas_tissues$tissue == tiss]))
sapply(twas_tissuegenes, function(x) length(x))
length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
gene_IDs_metaxscan <- unique(metaxscan_results$gene)
gene_IDs_DE <- unique(do.call(rbind, deg_eqtl_list)$gene_id)
length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
length(setdiff(gene_IDs_DE, gene_IDs_metaxscan)) / length(gene_IDs_DE)

#now let's use the TWAS results from Barbeira et al. 2021
twas_results_directory <- "~/data/smontgom/eqtl/"
tissue = motrpac_gtex_map[6]
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
gwas_names <- stringr::str_replace_all(gwas_names, "imputed_", "")
trait <- gwas_names[grep("hyperte", gwas_names, T)]
twas <- fread(file = paste0(twas_results_directory, "spredixcan_igwas_gtexmashrv8_", trait, "__PM__", tissue, ".csv"))

jac <- function(x1, x2){length(intersect(x1, x2)) / length(union(x1, x2))}
jac(unique(deg_eqtl_list[[names(tissue)]]$gene_name), twas$gene_name)
jac(unique(deg_eqtl_list[[names(tissue)]]$gene_id), twas$gene)
1 - (length(setdiff(unique(deg_eqtl_list[[names(tissue)]]$gene_id[deg_eqtl_list[[names(tissue)]]$selection_fdr < 0.05]), twas$gene)) / 
       length(unique(deg_eqtl_list[[names(tissue)]]$gene_id[deg_eqtl_list[[names(tissue)]]$selection_fdr < 0.05])))

available_twas_traits <- list.files(path = twas_results_directory)
available_twas_tissues <- table(sapply(available_twas_traits, function(x) strsplit(x, "__PM__")[[1]][2]))
available_twas_traits <- table(sapply(gsub(x = available_twas_traits, "spredixcan_igwas_gtexmashrv8_", ""), function(x) strsplit(x, "__P")[[1]][1]))

# check cluster membership vs. FDR
# load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
# load('~/data/smontgom/dea/transcript_rna_seq_20210804.RData')
# cluster_membership <- cluster_membership[cluster_membership$cluster %in% c(1,3,7,15) & cluster_membership$assay == "TRNSCRPT",]
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
# rna_dea_list <- lapply(setNames(unique(rna_dea$timewise_dea$tissue_abbreviation), unique(rna_dea$timewise_dea$tissue_abbreviation)), 
#                        function(tissue_abbrev) rna_dea$timewise_dea[rna_dea$timewise_dea$tissue_abbreviation == tissue_abbrev,])
# 
# selection_FDRs <- sapply(1:nrow(cluster_membership), function(ri)
#   rna_dea_list[[cluster_membership$tissue[ri]]]$selection_fdr[rna_dea_list[[cluster_membership$tissue[ri]]]$feature_ID == cluster_membership$feature_ID[ri]][1])
# hist(selection_FDRs, breaks = 0:100/100)
# mean(selection_FDRs > 0.1)

#compile the twas into one big data frame?
all_twas <- lapply(setNames(gwas_names, gwas_names), function(trait) data.table());
for(trait in gwas_names){
  print(trait)
  all_tissues <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), function(tissue) data.table());
  for(tissue in names(motrpac_gtex_map)){
    twas <- fread(file = paste0(twas_results_directory, "spredixcan_igwas_gtexmashrv8_", trait, "__PM__", motrpac_gtex_map[names(motrpac_gtex_map) == tissue], ".csv"))
    twas <- twas[,c("gene", "gene_name", "zscore", "pvalue")]
    twas <- twas[-which(is.na(twas$zscore) & is.na(twas$pvalue)),]
    twas$tissue <- tissue
    twas$trait <- trait
    all_tissues[[tissue]] <- twas
  }
  all_twas[[trait]] <- do.call(rbind, all_tissues)
}



data <- data.frame(trait = names(sort(sapply(gwas_names, function(trait) 1-qvalue::pi0est(all_twas[[trait]]$pvalue)$pi0))),
                   pi1 = sort(sapply(gwas_names, function(trait) 1-qvalue::pi0est(all_twas[[trait]]$pvalue)$pi0)))
ggplot(data, aes(x=trait, y=pi1)) + 
  geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
  coord_flip() + scale_x_discrete(limits=data$trait)

all.twas <- do.call(rbind, all_twas)
1-qvalue::pi0est(all.twas$pvalue)$pi0
hist(all.twas$pvalue)

#filter by category
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
all.twas$trait_category <- traitwise_partitions$Category[match(all.twas$trait, traitwise_partitions$Tag)]
some.twas <- all.twas[all.twas$trait_category %in% c("Cardiometabolic", "Aging", "Anthropometric", "Immune")]
hist(some.twas$pvalue)
hist(log(all.twas$pvalue), breaks = -700:1)
1-qvalue::pi0est(some.twas$pvalue)$pi0
sum(some.twas$pvalue < 0.05 / nrow(some.twas))

data <- data.frame(Trait = names(sort(table(some.twas[p.adjust(some.twas$pvalue, method = "bonf") < 0.05,]$trait))),
                   n_sig_bonferroni = as.integer(sort(table(some.twas[p.adjust(some.twas$pvalue, method = "bonf") < 0.05,]$trait))))
ggplot(data, aes(x=Trait, y=n_sig_bonferroni)) + 
  geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
  coord_flip() + scale_x_discrete(limits=data$Trait)                      

data <- data.frame(trait = names(sort(table(all.twas[all.twas$pvalue < 0.05 / nrow(all.twas),]$trait))),
                   n_sig = as.integer(sort(table(all.twas[all.twas$pvalue < 0.05 / nrow(all.twas),]$trait))))
ggplot(data, aes(x=trait, y=n_sig)) + 
  geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
  coord_flip() + scale_x_discrete(limits=data$trait)                      

data <- data.frame(trait = names(sort(table(some.twas[some.twas$pvalue < 0.05 / nrow(some.twas),]$gene_name))),
                   n_sig = as.integer(sort(table(some.twas[some.twas$pvalue < 0.05 / nrow(some.twas),]$gene_name))))
ggplot(data, aes(x=trait, y=n_sig)) + 
  geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
  coord_flip() + scale_x_discrete(limits=data$trait)    

#ok now let's try to x-reference it to the differential expression results!
DEsub <- deg_eqtl_list[[tissue]]
hist(DEsub$selection_fdr)
sum(p.adjust(some.twas$pvalue, method = "BH") < 0.05) / length(some.twas$pvalue)
some.twas$pvalue_BH <- p.adjust(some.twas$pvalue, method = "BH")

salient_twas <- unique(some.twas$trait)
deg_eqtl_list_TWAS <- deg_eqtl_list
for(tissue_i in names(deg_eqtl_list_TWAS)){
  
  print(tissue_i)
  # twas_addition <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas)))
  # twas_addition_nominalPValue <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas))) 
  # twas_addition_zscore <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas)))
  # colnames(twas_addition) <- colnames(twas_addition_zscore) <- salient_twas
  twas_sub <- some.twas[tissue == tissue_i]
  
  gtex_motrpac <- rna_dea$timewise_dea[tissue == tissue_i]
  for(trait_j in salient_twas){
    # twas_sub_sub <- twas_sub[trait == trait_j]
    # twas_addition[,trait_j] <- twas_sub_sub$pvalue_BH[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
     
    # subset in TWAS
    twas_sub_sub <- twas_sub[trait == trait_j]
    twas_sub_sub[, human_ensembl_gene := gsub('\\..*','',gene)]
    twas_sub_sub$human_gene_symbol <- twas_sub_sub$gene_name
    p_twas_in_map <- round(mean(unique(twas_sub_sub$human_gene_symbol) %in% unique(map$human_gene_symbol)), 3)
    
    # match human genes with rat ensembl genes
    twas_sub_sub = merge(twas_sub_sub, map, by='human_gene_symbol')
    twas_sub_sub <- twas_sub_sub[,c("feature_ID", "gene_name", "pvalue", "zscore")]
    colnames(twas_sub_sub)[-match(c("feature_ID", "gene_name"), colnames(twas_sub_sub))] <- paste0(trait_j, ".", 
                        colnames(twas_sub_sub)[-match(c("feature_ID", "gene_name"), colnames(twas_sub_sub))])
    
    matched_twas <- twas_sub_sub[match(gtex_motrpac$feature_ID, twas_sub_sub$feature_ID),-c("feature_ID")]
    if(trait_j != salient_twas[1]){matched_twas <- matched_twas[,-"gene_name"]}
    gtex_motrpac <- cbind(gtex_motrpac, matched_twas)

    cat(paste0("\nprop of twas genes in map: ", p_twas_in_map,
                 "\nprop of resulting feature_IDs matched to GTEx: ", round(mean(!is.na(matched_twas[,1])), 3)))

    # twas_addition_nominalPValue[,trait_j] <- twas_sub_sub$pvalue[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
    # twas_addition_zscore[,trait_j] <- twas_sub_sub$zscore[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
  }
  deg_eqtl_list_TWAS[[tissue_i]] <- gtex_motrpac
  
  # colnames(twas_addition) <- paste0(colnames(twas_addition), "_BH_PValue")
  # colnames(twas_addition_nominalPValue) <- paste0(colnames(twas_addition_nominalPValue), "_Nominal_PValue")
  # colnames(twas_addition_zscore) <- paste0(colnames(twas_addition_zscore), "_zscore")
  # deg_eqtl_list_TWAS[[tissue_i]] <- cbind(deg_eqtl_list_TWAS[[tissue_i]], twas_addition, twas_addition_nominalPValue, twas_addition_zscore)

}

#get number of p-values
# deg_eqtl_list_TWAS_all <- do.call(rbind, deg_eqtl_list_TWAS)
# deg_eqtl_list_TWAS_all <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$sex == "female"]
# deg_eqtl_list_TWAS_all <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$comparison_group == "8w"]
# 
# all_twas_pvalues <- unlist(as.data.frame(deg_eqtl_list_TWAS_all)[,grep(colnames(deg_eqtl_list_TWAS_all), pattern = ".pvalue")])
# all_twas_pvalues <- all_twas_pvalues[!is.na(all_twas_pvalues)]
# n_twas_comparisons <- length(all_twas_pvalues)

#first do simple gene intersect
if(!exists("cluster_membership")){
  load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
  cluster_membership <- cluster_membership[cluster_membership$cluster %in% c(1,3,7,15) & cluster_membership$assay == "TRNSCRPT",]
  map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
  cluster_genes <- map$human_gene_symbol[match(cluster_membership$feature_ID, map$feature_ID)]
  cluster_genes <- cluster_genes[!is.na(cluster_genes)]
}

# deg_eqtl_list_TWAS_all_1.3.7.15 <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$gene_name %in% cluster_genes,]
# all_twas_pvalues_clusters_1.3.7.15 <- unlist(as.data.frame(deg_eqtl_list_TWAS_all_1.3.7.15)[,grep(colnames(deg_eqtl_list_TWAS_all_1.3.7.15_subsetSexTime), pattern = ".pvalue")])
# all_twas_pvalues_clusters_1.3.7.15 <- all_twas_pvalues_clusters_1.3.7.15[!is.na(all_twas_pvalues_clusters_1.3.7.15)]
# n_twas_comparisons_clusters_1.3.7.15 <- length(all_twas_pvalues_clusters_1.3.7.15)

#hmm wait that is not quite right... let's subset it to the cluster-tissue subset first, and then further subset it to the FDR genes
deg_eqtl_list_TWAS_cluster_subset <- lapply(deg_eqtl_list_TWAS[-which(names(deg_eqtl_list_TWAS) %in% c("t1000-gonads", "t63-testes", "t64-ovaries"))], function(tissue_dea) 
    tissue_dea[tissue_dea$feature_ID %in% cluster_membership$feature_ID[cluster_membership$tissue == tissue_dea$tissue_abbreviation[1] & 
                                                                        cluster_membership$ome == "TRNSCRPT" & 
                                                                        cluster_membership$cluster %in% c(1,3,7,15)],] 
    #& tissue_dea$selection_fdr < 0.1,]
  )

deg_eqtl_list_TWAS_cluster_subset_ntests <- do.call(rbind, deg_eqtl_list_TWAS_cluster_subset)
deg_eqtl_list_TWAS_cluster_subset_ntests <- deg_eqtl_list_TWAS_cluster_subset_ntests[deg_eqtl_list_TWAS_cluster_subset_ntests$sex == "female" &
                                                                                       deg_eqtl_list_TWAS_cluster_subset_ntests$comparison_group == "8w",]
twas_pvalues_clusters <- unlist(as.data.frame(deg_eqtl_list_TWAS_cluster_subset_ntests)[,grep(colnames(deg_eqtl_list_TWAS_cluster_subset_ntests), pattern = ".pvalue")])
twas_pvalues_clusters <- twas_pvalues_clusters[!is.na(twas_pvalues_clusters)]
n_twas_comparisons_clusters <- length(twas_pvalues_clusters)

n_deg_sigtwas_intersect <- as.data.frame(matrix(0, nrow = length(deg_eqtl_list_TWAS_cluster_subset), ncol = length(salient_twas), dimnames = list(names(deg_eqtl_list_TWAS_cluster_subset), salient_twas)))
for(tissue_i in names(deg_eqtl_list_TWAS_cluster_subset)){
  
  print(tissue_i)
  DELT <- as.data.frame(deg_eqtl_list_TWAS_cluster_subset[[tissue_i]])
  
  DE_inds <- which(DELT$selection_fdr <= 1) #since we're subsetting to monotonic clusters already, but why not
  TWAS_inds <- apply(log(DELT[,grep(colnames(DELT), pattern = ".pvalue")]) <= (log(0.05) - log(n_twas_comparisons_clusters)), 2, which)
  
  intersect_inds <- lapply(TWAS_inds, function(twi) intersect(DE_inds, twi))
  intersect_genes <- lapply(intersect_inds, function(ii) unique(DELT$gene_name[ii]))
  # intersect_genes <- lapply(intersect_genes, function(ii) intersect(ii, cluster_genes)) #subset to just monotonic sex-homogenous clusters
  
  n_intersect <- sapply(intersect_genes, function(ig) length(ig))
  names(n_intersect) <- gsub(names(n_intersect), pattern = ".pvalue", replacement = "")
  n_deg_sigtwas_intersect[tissue_i,names(n_intersect)] <- n_intersect

}

n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,order(apply(n_deg_sigtwas_intersect, 2, sum), decreasing = T)]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,apply(n_deg_sigtwas_intersect, 2, sum) != 0]

#now plot it
par(xpd = T)
plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", xlim= c(-5,ncol(n_deg_sigtwas_intersect)), ylim = c(-5,nrow(n_deg_sigtwas_intersect)))
heatmap_cols <- viridisLite::viridis(n = max(n_deg_sigtwas_intersect)*100+1)
heatmap_cols <- heatmap_cols[round(log(1:max(n_deg_sigtwas_intersect)) / log(max(n_deg_sigtwas_intersect)) * max(n_deg_sigtwas_intersect) * 100 + 1)]
for(ri in 1:nrow(n_deg_sigtwas_intersect)){
  text(x = 0.5, y = ri, pos = 2, labels = rownames(n_deg_sigtwas_intersect)[ri])
  for(ci in 1:ncol(n_deg_sigtwas_intersect)){
    if(ri == 1){text(x = ci+0.75, y = 0.2, pos = 2, srt = 45,
                     labels = trait_categories$new_Phenotype[match(colnames(n_deg_sigtwas_intersect)[ci], trait_categories$Tag)])}
    rect(xleft = ci + 1/2,
         xright = ci - 1/2,
         ybottom = ri - 1/2,
         ytop =  ri + 1/2,
         col = heatmap_cols[n_deg_sigtwas_intersect[ri, ci]])
    text(n_deg_sigtwas_intersect[ri, ci], x = ci, y = ri, col = "white", cex = 0.85)
  }
}
#legend
for(i in 0:max(n_deg_sigtwas_intersect)){
  if(i %in% round(seq(0, max(n_deg_sigtwas_intersect), length.out = 10))){
    text(labels = i, x = ncol(n_deg_sigtwas_intersect) + 2.25, pos = 4, cex = 0.75,
       y = (i) * nrow(n_deg_sigtwas_intersect) / (max(n_deg_sigtwas_intersect)+1) + 0.675)
  }
  rect(xleft = ncol(n_deg_sigtwas_intersect) + 2 + 1/2,
       xright = ncol(n_deg_sigtwas_intersect) + 2 - 1/2,
       ybottom = (i - 1/2) * nrow(n_deg_sigtwas_intersect) / (max(n_deg_sigtwas_intersect)+1) + 0.75,
       ytop =  (i + 1/2) * nrow(n_deg_sigtwas_intersect) / (max(n_deg_sigtwas_intersect)+1) + 0.75,
       col = heatmap_cols[i]) 
}
text(labels = latex2exp::TeX("n_{intersect}"), pos = 3, font = 2,
     x = ncol(n_deg_sigtwas_intersect) + 2, y = nrow(n_deg_sigtwas_intersect) + 0.5)
text(latex2exp::TeX("number of genes in monotonic, sex-homogenous clusters and w/ Bonferroni significant TWAS at $\\alpha$ = 0.05"), 
     x = 22, y = nrow(n_deg_sigtwas_intersect) + 0.5, pos = 3, cex = 1.25, font = 2)








#now do the proportion
deg_sigtwas_proportion <- array(NA, dim = c(length(deg_eqtl_list_TWAS_cluster_subset), length(salient_twas), 4, 3, 2), 
                                dimnames = list(names(deg_eqtl_list_TWAS_cluster_subset), salient_twas, paste0(2^(0:3), "w"), c("p", "n", "genes"), c("male", "female")))
nuniq <- function(x) length(unique(x))
for(tissue_i in names(deg_eqtl_list_TWAS_cluster_subset)){
  
  print(tissue_i)
  
  for(time in paste0(2^(0:3), "w")){
    
    for(sex in c("male", "female")){
      
      DELT <- as.data.frame(deg_eqtl_list_TWAS_cluster_subset[[tissue_i]])
      
      #subset to time, sex, and unique gene entries
      DELT <- DELT[DELT$comparison_group == time,]
      DELT <- DELT[DELT$sex == sex,]
      if(nrow(DELT) == 0){next()}
      DELT <- DELT[match(unique(DELT$gene_name), DELT$gene_name),] #pull out only first gene entry
      
      DE_inds <- which(DELT$adj_p_value <= 1)
      # TWAS_inds <- apply(DELT[,grep(colnames(DELT), pattern = "BH_PValue")] < 0.05, 2, which)
      TWAS_inds <- apply(log(DELT[,grep(colnames(DELT), pattern = ".pvalue")]) <= (log(0.05) - log(n_twas_comparisons_clusters)), 2, which)
      
      intersect_inds <- lapply(TWAS_inds, function(twi) intersect(DE_inds, twi))
      intersect_inds <- lapply(intersect_inds, function(ii) ii[DELT$gene_name[ii] %in% cluster_genes]) #subset to just monotonic sex-homogenous clusters
      intersect_genes <- lapply(intersect_inds, function(ii) unique(DELT$gene_name[ii]))
      
      
      intersect_sign_logFC <- lapply(intersect_inds, function(ii) sign(DELT$logFC[ii]))
      names(intersect_sign_logFC) <- gsub(names(intersect_sign_logFC), pattern = ".pvalue", replacement = "")
      intersect_sign_TWAS <- lapply(gsub(names(intersect_inds), pattern = ".pvalue", replacement = ""), function(ii) 
                             sign(DELT[intersect_inds[paste0(ii, ".pvalue")][[1]], paste0(ii, ".zscore")]))
      intersect_sign_match <- lapply(1:length(intersect_sign_logFC), function(i) intersect_sign_logFC[[i]]*intersect_sign_TWAS[[i]])
      
      intersect_genes_strings <- lapply(1:length(intersect_sign_logFC), function(i) 
        paste0(intersect_genes[[i]], " (", c("-", "", "+")[intersect_sign_match[[i]] + 2], ")", collapse = " ~ "))
      intersect_genes_strings <- unlist(intersect_genes_strings)
      intersect_genes_strings[intersect_genes_strings == " ()"] <- ""
      
      intersect_sign_match <- sapply(intersect_sign_match, function(x) mean(x == 1))
      intersect_sign_match_n <- sapply(intersect_inds, function(x) length(x))
      
      deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "p", sex] <- intersect_sign_match
      deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "n", sex] <- intersect_sign_match_n
      deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "genes", sex] <- intersect_genes_strings
    
    }
    
  }
  
}

deg_sigtwas_proportion[,salient_twas[grep(salient_twas, pattern = "type_1_diabetes", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,salient_twas[grep(salient_twas, pattern = "reported_hypertension", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,salient_twas[grep(salient_twas, pattern = "heart_problems", ignore.case = T)],,"n",sex]
deg_sigtwas_proportion[,salient_twas[grep(salient_twas, pattern = "Body_mass_index_BMI", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,salient_twas[grep(salient_twas, pattern = "Body_fat_percentage", ignore.case = T)],,"p",sex]

#plot lines for proportions
trait <- salient_twas[grep(salient_twas, pattern = "reported_hypertension", ignore.case = T)][1]
trait <- salient_twas[grep(salient_twas, pattern = "Body_fat_percentage", ignore.case = T)][1]

cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))




par(mfrow = c(1,2), xpd = NA, mar = c(5,5,5,17))
lwd <- 3
for(sex in c("male", "female")){
  
  #retrieve data
  d <- apply(deg_sigtwas_proportion[, trait,,"p",sex], 2, as.numeric)
  dn <-  apply(deg_sigtwas_proportion[, trait,,"n",sex], 2, as.numeric)
  colnames(d) <- colnames(dn) <- colnames(deg_sigtwas_proportion[, trait,,"p",sex])
  rownames(d) <- rownames(dn) <- rownames(deg_sigtwas_proportion[, trait,,"p",sex])
  dn <- dn[complete.cases(d),]
  d <- d[complete.cases(d),]
  dg <- deg_sigtwas_proportion[, trait,,"genes",sex]
  dg <- lapply(dg[,"8w"], function(x) strsplit(x, " ~ ")[[1]])
  
  #iterate through d to make identical lines parallel
  line_thickness <- lwd / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
  need_to_increment <- matrix(T, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
  while(any(need_to_increment)){
    need_to_increment <- matrix(F, nrow = nrow(d), ncol = ncol(d), dimnames = list(rownames(d), colnames(d)))
    for(coi in 1:(ncol(d)-1)){
      unchanging_tissues <- rownames(unique(d[,c(coi,coi+1)]))
      need_to_increment[setdiff(rownames(d), unchanging_tissues),c(coi,coi+1)] <- T
    }
    d[need_to_increment] <- d[need_to_increment] + line_thickness * 1.05
  }
  
  #find coordinates to plot tissue names
  xylocs_tissue_names <- FField::FFieldPtRep(coords = cbind(rep(4 + 0.25, nrow(d)), 
                                                            c(d[,"8w"] * 100 + rnorm(nrow(d), 0, 1E-2))),
                                             rep.fact = 20, adj.max = 1)
  xylocs_tissue_names$y <- xylocs_tissue_names$y / 100
  
  #find coordinates to plot gene names
  xylocs_tissues_genes <- data.frame(gene = unlist(dg[rownames(xylocs_tissue_names)[order(xylocs_tissue_names$y, decreasing = T)]]))
  xylocs_tissues_genes$tissue <- rownames(xylocs_tissue_names)[sapply(rownames(xylocs_tissues_genes), function(x) 
    grep(strsplit(x, "-")[[1]][1], rownames(xylocs_tissue_names)))]
  xylocs_tissues_genes$x <- 5.75
  # xylocs_tissues_genes$y <- xylocs_tissue_names$y[match(xylocs_tissues_genes$tissue, rownames(xylocs_tissue_names))]
  # xylocs_tissues_genes$y <- xylocs_tissues_genes$y + 1:length(xylocs_tissues_genes$y)/1E4
  xylocs_tissues_genes$y <- seq(max(xylocs_tissue_names$y) + 0.05, min(xylocs_tissue_names$y) - 0.1, length.out = nrow(xylocs_tissues_genes))
  # xylocs_tissues_genes_locs <- FField::FFieldPtRep(coords = cbind(xylocs_tissues_genes$x, 
  #                                                                 xylocs_tissues_genes$y * 100), 
  #                                                  rep.fact = 20, adj.max = 1, iter.max = 5E3)
  # xylocs_tissues_genes$y <- xylocs_tissues_genes_locs$y / 100
  
  plot(100,100,xlim = c(1,4), ylim = c(0,1), xpd=NA, 
       ylab = "Proportion Positive Effects on GWAS Trait", xlab = "Timepoint", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)
  
  
  #plot faded positive and negative regions
  rect(xl = 1, xr = 4, yb = 0.5, ytop = 1,
       col = grDevices::adjustcolor("red", 0.1), border = NA)
  rect(xl = 1, xr = 4, yb = 0, ytop = 0.5,
       col = grDevices::adjustcolor("blue", 0.1), border = NA)
  
  text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 3, pos = 3,
       x = par("usr")[2] * 0.5 + par("usr")[1] * 0.5, y = par("usr")[3] * 0 + par("usr")[4] * 1)
  
  #horizontal axis
  segments(x0 = 1:4, x1 = 1:4, y0 = - 0.02, y1 = - 0.04)
  segments(x0 = 1, x1 = 4, y0 = - 0.02, y1 = - 0.02)
  text(x = 1:4, y = - 0.04, labels = paste0(2^(0:3), "w"), pos = 1)
  
  #vertical axis
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5)
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1)
  text(x = 1 - 3 * 0.02, y = 0:5/5, labels = 0:5/5, pos = 2)
  
  
  for(tissue in rownames(d)){
    
    if(all(is.na(d[tissue,]))){
      next()
    }
    
    lines(1:4, d[tissue,], lwd = lwd, col = cols$Tissue[tissue])
    
    #plot tissue names
    text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = 1,
         y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
         labels = paste0(tissue, " (", dn[tissue,"8w"], ")"), col = cols$Tissue[tissue], pos = 4)
    segments(x0 = 4, y0 = d[tissue,"8w"], 
             x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.075, 
             y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
             col = cols$Tissue[tissue], lty = 3)
    
    #plot gene names
    tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
    text(labels = tissue_genes$gene,
         x = tissue_genes$x,
         y = tissue_genes$y,
         cex = 0.875, col = cols$Tissue[tissue], pos = 4)
    for(gene_i in 1:nrow(tissue_genes)){
      segments(x0 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 
                 strwidth(paste0(tissue, " (", dn[tissue,"8w"], ")   "), units = "user"), 
               y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
               x1 = tissue_genes$x[gene_i] + 0.75*strwidth("  ", units = "user"),
               y1 = tissue_genes$y[gene_i],
               col = cols$Tissue[tissue], lty = 3)
    }
    
  }
  # legend(x = 1, y = 1.5, legend = tissue_names, 
  #        col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
  # segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")
  
}
text(trait_categories$new_Phenotype[trait_categories$Tag == trait], x =0, y = 1.175, cex = 2, font = 2)


#make data table for nicole 
# save(deg_eqtl_list_TWAS, file = "~/data/smontgom/deg_eqtl_twas.RData")

# plot(p.adjust(do.call(rbind, deg_eqtl_list_TWAS)$UKB_20002_1223_self_reported_type_2_diabetes_Nominal_PValue, method = "bonf"),
# do.call(rbind, deg_eqtl_list_TWAS)$UKB_20002_1223_self_reported_type_2_diabetes_BH_PValue)

#### GREX, or just expression data works too ####
# load libraries
library(foreach)
library(doParallel)
library(parallel)

#initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

library(data.table)

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
tissues <- names(motrpac_gtex_map)
sds_expression <- list()
for(tissue in tissues){
  print(tissue)
  d <- fread(file = paste0("~/repos/gtex-pipeline/log2-normalized-expression/log2-normalized-expression_",tissue,".expression.bed.gz"))
  d$gene_id <- gsub('\\..*','',d$gene_id)
  
  # ec <- fread(file = paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_",tissue,"_expected_count.gct.gz"))
  # ec$gene_id <- gsub('\\..*','',ec$gene_id)
  # 
  # tpm <- fread(file = paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_",tissue,"_tpm.gct.gz"))
  # tpm$gene_id <- gsub('\\..*','',tpm$gene_id)
  # 
  # diff(as.integer(sapply(colnames(d), function(name) grep(name, colnames(ec)))))
  # diff(as.integer(sapply(colnames(d), function(name) grep(name, colnames(tpm)))))
  # 
  # plot(as.numeric(d[1,-c(1:4)]), log2(as.numeric(ec[match(d$gene_id[1], ec$gene_id),-1])+1))
  # plot(as.numeric(d[3,-c(1:4)]), log2(as.numeric(tpm[match(d$gene_id[3], tpm$gene_id),-1])+1))
  # plot(as.numeric(tpm[1,-1]), as.numeric(ec[1,-1]))

  
  indivs <- colnames(d)[grep("GTEX", colnames(d))]
  sds_expression[[tissue]] <- setNames(apply(d[,..indivs], 1, sd), d$gene_id)
}
# sapply(sds_expression, function(x) length(x))
if(all(sapply(tissues, function(ti) all(names(sds_expression[[1]]) == names(sds_expression[[ti]]))))){
  sds_expression <- do.call(rbind, sds_expression)
}
rownames(sds_expression) <- tissues
GTEx_SampleSize <- read.csv("~/data/smontgom/GTEx_SampleSizes.csv", header = T)
colnames(GTEx_SampleSize) <- gsub("X..", "", colnames(GTEx_SampleSize))
GTEx_SampleSize <- data.frame(tissue = names(motrpac_gtex_map), sample_size = GTEx_SampleSize$RNASeq.and.Genotyped.samples[match(motrpac_gtex_map, GTEx_SampleSize$Tissue)])
hist(do.call(rbind, sds_expression))
library(arrow)
ldsc_directory <- "~/repos/ldsc/"

process_GWASy_sumstats <- F
if(process_GWASy_sumstats){
  foreach(tissue=tissues, .packages = c("data.table", "arrow")) %dopar% {
  # for(tissue in tissues){
    print(tissue)
    if(!dir.exists(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue))){dir.create(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue))}
    
    for(cri in c(1:22,"X")){
      
      cat(paste0(" ", cri))
      if(!dir.exists(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri))){dir.create(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri))}
      
      RSID_POS_MAP <- fread(paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
      
      eQTL_sumstats <- read_parquet(paste0("~/repos/gtex-pipeline/tensorQTL_output/log2-normalized-expression_", tissue,".cis_qtl_pairs.chr", cri, ".parquet"))
      eQTL_sumstats$ENSG <- gsub('\\..*','',eQTL_sumstats$phenotype_id)
      
      vids <- eQTL_sumstats$variant_id
      vids <- substr(vids, start = 4, nchar(vids))
      vids <- strsplit(vids, "_", T)
      vids <- as.data.table(data.frame(data.table::transpose(vids)))
      colnames(vids) <- c("CHROM", "POS", "REF", "ALT", "BUILD")
      vids$RSID <- RSID_POS_MAP$ID[match(vids$POS, RSID_POS_MAP$POS)]
      eQTL_sumstats <- cbind(eQTL_sumstats, vids)
      eQTL_sumstats <- eQTL_sumstats[!is.na(eQTL_sumstats$RSID),]
      eQTL_sumstats$Z <- round(eQTL_sumstats$slope / eQTL_sumstats$slope_se, 5)
      eQTL_sumstats$N <-  paste0(GTEx_SampleSize$sample_size[GTEx_SampleSize$tissue == tissue], ".00")
      
      genes <- unique(eQTL_sumstats$ENSG)
      eQTL_sumstats <- split(eQTL_sumstats , f = eQTL_sumstats$ENSG)
      eQTL_sumstats <- lapply(eQTL_sumstats, function(gid){
                              gene_info <- gid[,c("RSID", "REF", "ALT", "Z", "N")]
                              colnames(gene_info) <- c("SNP", "A1", "A2", "Z", "N")
                              return(gene_info)
      })
      
      for(gene_name in names(eQTL_sumstats)){
        
        fwrite(eQTL_sumstats[[gene_name]], file = paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri, "/", gene_name, ".sumstats.gz"), sep = "\t", append = F, col.names = T)
      }
      
    }
  }
}

#compare to heritability estimates from https://github.com/WheelerLab/GenArchDB
library(data.table)
h2dir <- "~/repos/GenArchDB/GenArch_reml-no-constrain_h2/"
h2_files <- list.files(h2dir)
h2_files <- h2_files[intersect(grep("GTEx", h2_files), grep(".TS.", h2_files))]
intersect_rec <- function(x) {
  if(length(x) <= 1){
    return(x)
  }else if(length(x) > 2){
    x[[1]] <- intersect(x[[1]], x[[2]])
    intersect_rec(x[-2])
  } else {
    intersect(x[[1]], x[[2]])
  }
}
h2_tissue_to_file <- sapply(motrpac_gtex_map, function(tissue) h2_files[intersect_rec(sapply(strsplit(tissue, "_")[[1]], function(x) grep(x, h2_files)))][1])
tissue = names(motrpac_gtex_map)[7]
h2 <- fread(paste0(h2dir, h2_tissue_to_file[tissue]))
h2$gene_id <- gsub('\\..*','',h2$ensid)

#compare to gcta heritability estimates
gcta_directory <- "~/repos/gcta_1.93.2beta_mac/"
load(paste0(gcta_directory, "gcta_output_GTEx_allTissues.RData")) #gcta_output

par(mfrow = c(4,4))
for(tissue in setdiff(names(motrpac_gtex_map), "t59-kidney")){
  h2 <- fread(paste0(h2dir, h2_tissue_to_file[tissue]))
  h2$gene_id <- gsub('\\..*','',h2$ensid)
  h2s <- h2$global.h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, h2$gene_id)]
  prop_pos <- round(sum(h2s > 0, na.rm = T) / sum(!is.na(h2s)) * 100, 1)
  hist(h2s, main= paste0(tissue, " (", prop_pos, "% positive)"), xlab = "heritability")
}
for(tissue in setdiff(names(motrpac_gtex_map), "t59-kidney")){
  h2 <- as.data.frame(fread(paste0(h2dir, h2_tissue_to_file[tissue])))
  h2$gene_id <- gsub('\\..*','',h2$ensid)
  h2$myh2 <- gcta_output[[tissue]]$h2[match(h2$gene_id, gcta_output[[tissue]]$ENSG)]
  h2s <- h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x[deg_eqtl_list[[tissue]]$comparison_group == "1w" & deg_eqtl_list[[tissue]]$sex == deg_eqtl_list[[tissue]]$sex[1]], 
                  h2$gene_id), c("local.h2", "myh2")]
  plot(h2s[complete.cases(h2s),], main= paste0(tissue, " (r = ", round(cor(h2s[complete.cases(h2s),])[1,2], 3), ")"), 
       xlab = "wheeler h2", ylab = "my gcta h2", pch = 19, col = adjustcolor(1, 0.1))
}

#compute relative exercise effect z-scores

for(tissue in tissues){
  print(tissue)
  deg_eqtl_list[[tissue]]$phenotypic_expression_Z <- deg_eqtl_list[[tissue]]$logFC / sds_expression[tissue, match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, colnames(sds_expression))]
  deg_eqtl_list[[tissue]]$genetic_expression_Z <- deg_eqtl_list[[tissue]]$logFC / 
    (sds_expression[tissue, match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, colnames(sds_expression))] * 
       sqrt(gcta_output[[tissue]]$h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, gcta_output[[tissue]]$ENSG)]))
  deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z <- deg_eqtl_list[[tissue]]$logFC / 
    (sds_expression[tissue, match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, colnames(sds_expression))] * 
       sqrt(gcta_output[[tissue]]$h2[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, gcta_output[[tissue]]$ENSG)] + 
              2 * gcta_output[[tissue]]$SE[match(deg_eqtl_list[[tissue]]$human_ensembl_gene.x, gcta_output[[tissue]]$ENSG)]))
}

#compute quantile plots
qs2use <- 1:9999/10000
EZ_PZ <- lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
  sapply(tissues, function(tissue) quantile(
    x = deg_eqtl_list[[tissue]]$phenotypic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
    probs = qs2use, na.rm = T)))
EZ_PZ <- lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
  sapply(tissues, function(tissue) quantile(
    x = deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
    probs = qs2use, na.rm = T)))
EZ_PZ <- lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
  sapply(tissues, function(tissue) quantile(
    x = (deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z[deg_eqtl_list[[tissue]]$comparison_group == ti]), 
    probs = qs2use, na.rm = T)))

tissue <- tissues[3]
plot(deg_eqtl_list[[tissue]]$logFC, deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z)

#### just the plotting ####
#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}

ti = "8w"
f_p <- 0.4
f_x <- 2
ylims = c(min(sort(EZ_PZ[[ti]])[sort(EZ_PZ[[ti]]) != -Inf]),max(sort(EZ_PZ[[ti]],T)[sort(EZ_PZ[[ti]],T) != Inf]))
ylims <- squish_middle_x(ylims, f_x)

plot(100,100,xlim = c(0,1.25), ylim = ylims, xpd=NA, ylab = "Z-Score",
     main = latex2exp::TeX("Ratio of Exercise DE to \\sqrt{Variance in log_2(Gene Expression)}"), xlab = "Quantile", xaxt = "n", yaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25)

xylocs_tissue_names <- cbind(rep(1.05, ncol(EZ_PZ[[ti]])), redistribute(as.numeric(squish_middle_x(tail(EZ_PZ[[ti]], 1), f_x)), diff(ylims) / 30))
rownames(xylocs_tissue_names) <- colnames(EZ_PZ[[ti]]); colnames(xylocs_tissue_names) <- c("x", "y")

#horiz axis
segments(x0 = 0, y0 = ylims[1] - diff(ylims) / 100, x1 = 1, y1 = ylims[1] - diff(ylims) / 100, lwd = 2)
segments(x0 = 0:10/10, y0 = ylims[1] - diff(ylims) / 100, x1 = 0:10/10, y1 = ylims[1] - diff(ylims) / 50, lwd = 2, xpd = NA)
horiz_axis_labels <- round(unsquish_middle_p(0:10/10, f_p), 3)
horiz_axis_labels[1] <- 0; horiz_axis_labels[length(horiz_axis_labels)] <- 1;
text(labels = horiz_axis_labels, x = 0:10/10, y = rep(ylims[1] - diff(ylims) / 50, 10), pos = 1, xpd = NA)
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
  text(tissue, x = xylocs_tissue_names[tissue,"x"], y = xylocs_tissue_names[tissue,"y"], pos = 4, xpd = NA,
       col = cols$Tissue[tissue])
  segments(x0 = squish_middle_p(tail(qs2use, 1), f_p), x1 = xylocs_tissue_names[tissue,"x"]+0.01,
           y0 = squish_middle_x(tail(EZ_PZ[[ti]][,tissue], 1), f_x), y1 = xylocs_tissue_names[tissue,"y"],
           lty = 3, col = cols$Tissue[tissue])
}
shadowtext(x = xylocs_tissue_names[1,"x"], y = min(xylocs_tissue_names[,"y"]) - diff(ylims) / 7.5, labels = ti, 
           cex = 5, col = cols$Time[ti], pos = 4, r = 0.2) 
