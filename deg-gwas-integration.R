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
library(pracma)


#### define functions ####
polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
}
polarp <- function (t, r, col = 1, pch = 1, cex = 1, center = c(0,0)) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  points(xy[, 1] + center[1], xy[, 2] + center[2], cex = cex, col = col, pch = pch)
}
polarl <- function (t, r, col = 1, lwd = 1, center = c(0,0), adjx = 1) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  lines(xy[, 1] * adjx + center[1], xy[, 2] + center[2], lwd = lwd, col = col)
}
arc <- function(t1,t2,r1,r2,res=50,lwd = 1,col=1, mindist = T, self_adjust = 2 * pi / 45, random_selfing = T, clockwise_selfing = T, pointy_selfing = F,
                center = c(0,0), adjx = 1){
  if(mindist){
    if(abs(t1-t2) > pi){
      if(t1 > t2){
        t1 <- t1-2*pi
      } else {
        t2 <- t2-2*pi
      }
    }
  }
  if(abs(t1-t2) < 1E-6){
    
    if(random_selfing){
      lor <- sample(c(-1,1), 1)
    } else {
      lor <- ifelse(clockwise_selfing, -1, 1)
    }
    
    if(pointy_selfing){
      ts <- c(seq(t1, t1 + lor * self_adjust, length.out = res/2), seq(t1 + lor * self_adjust, t1, length.out = res/2))
      rs <- seq(r1, r2, length.out = res)
    } else {
      center <- center + pol2cart(cbind(t2, mean(c(r1, r2))))
      if(random_selfing){
        ts <- seq(t1, t1 + sample(c(-pi, pi), 1), length.out = res)
      } else {
        ts <- seq(t1, t1 + ifelse(clockwise_selfing, -pi, pi), length.out = res)
      }
      rs <- seq(max(c(r1,r2)) - mean(c(r1,r2)), mean(c(r1,r2)) - min(c(r1,r2)), length.out = res)
    }
    
  } else {
    ts <- seq(t1, t2, length.out = res)
    rs <- seq(r1, r2, length.out = res)
  }
  polarl(ts, rs, lwd = lwd,col=col, center = center, adjx = adjx)
}

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

GTEx_logo <- readJPEG("~/Documents/Documents - nikolai/GTEx_logo.jpg")
eqtl = '~/data/smontgom/GTEx_Analysis_v8_eQTL'
#deg = 'gs://mawg-data/pass1b-06/transcript-rna-seq/dea/'
deg = '~/data/smontgom/dea/'
#map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t')
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t', header=T)
# map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")
map <- unique(gene_map[,c("RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_ENSEMBL_ID", "HUMAN_ORTHOLOG_SYMBOL")])
colnames(map) <- c("feature_ID", "human_ensembl_gene", "human_gene_symbol")

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
tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
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
  ldsc_results$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 1:nrow(ldsc_results)){
    output <- fread(ldsc_results_paths[[ldsc_results$gwas[i]]][grep(pattern = paste0(ldsc_results$cluster[i], ".res"), 
                                                                    ldsc_results_paths[[ldsc_results$gwas[i]]])][1])
    ldsc_results[i,1:ncol(output)] <- output[grep(pattern = ldsc_results$cluster[i], output$Category),]
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

cols = list(Tissue=tissue_cols[tissue_abbrev], 
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
  # cluster_names <- paste0("Cluster_", 1:15)
  tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) ldsc_results_conditional_paths[grep(pattern = gwas, ldsc_results_conditional_paths)])
  names(ldsc_results_conditional_paths) <-  gwas_names
  ldsc_results_conditional_paths <- lapply(gwas_names, function(gwas) 
    ldsc_results_conditional_paths[[gwas]][grep(pattern = 
    paste0(sort(cluster_names)[c(1, length(cluster_names))], collapse = "-"), ldsc_results_conditional_paths[[gwas]])])
  names(ldsc_results_conditional_paths) <-  gwas_names
  
  ldsc_results_conditional <- as.data.frame(matrix(data = NA, ncol = ncol(fread(ldsc_results_conditional_paths[[1]][1])), nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_conditional) <- colnames(fread(ldsc_results_conditional_paths[[1]][1]))
  ldsc_results_conditional$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_conditional$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    output <- fread(ldsc_results_conditional_paths[[ldsc_results_conditional$gwas[i]]][1])
    output$Category <- stringr::str_replace(pattern = "L2_0", replacement = "", output$Category)
    output <- output[output$Category %in% paste0("Cluster_", cluster_names),]
    reorder_output <- order(match(stringr::str_replace(output$Category, pattern = "Cluster_", replacement = ""), 
          ldsc_results_conditional[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_conditional[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_conditional$gwas <- stringr::str_remove(ldsc_results_conditional$gwas, "imputed_")
  ldsc_results_conditional$logPVal_enrichment <- log10(ldsc_results_conditional$Enrichment_p)
  fwrite(ldsc_results_conditional, "~/data/smontgom/ldsc_cluster_results_conditional.txt")
}


#read in ldsc cts results outputted from *old* command
if(file.exists("~/data/smontgom/ldsc_cluster_results_cts_alt.txt")){
  ldsc_results_cts_alt <- as.data.frame(fread("~/data/smontgom/ldsc_cluster_results_cts_alt.txt"))
} else{
  ldsc_output_dir <- "~/repos/ldsc/output/original_command/"
  ldsc_results_cts_alt_paths <- list.files(ldsc_output_dir)
  ldsc_results_cts_alt_paths <- paste0(ldsc_output_dir, ldsc_results_cts_alt_paths[grep(ldsc_results_cts_alt_paths, pattern = "results")])
  # cluster_names <- paste0("Cluster_", 1:15)
  
  tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  
  ldsc_results_cts_alt_paths <- lapply(gwas_names, function(gwas) ldsc_results_cts_alt_paths[grep(pattern = gwas, ldsc_results_cts_alt_paths)])
  names(ldsc_results_cts_alt_paths) <-  gwas_names
  
  ldsc_results_cts_alt <- as.data.frame(matrix(data = NA, 
                                               ncol = ncol(fread(ldsc_results_cts_alt_paths[[1]][[1]][1])), 
                                               nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_cts_alt) <- colnames(fread(ldsc_results_cts_alt_paths[[1]][[1]][1]))
  ldsc_results_cts_alt$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_cts_alt$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    paths <- ldsc_results_cts_alt_paths[[ldsc_results_cts_alt$gwas[i]]]
    output <- do.call(rbind, lapply(setNames(cluster_names, cluster_names), function(cli){
      path <- paths[grep(cli, paths)]
      temp_output <- fread(path)
      temp_output <- temp_output[temp_output$Category == "L2_0"]
      temp_output$Category <- cli
      return(temp_output)
    }))
    reorder_output <- order(match(output$Category, ldsc_results_cts_alt[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_cts_alt[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_cts_alt$gwas <- stringr::str_remove(ldsc_results_cts_alt$gwas, "imputed_")
  ldsc_results_cts_alt$logPVal_enrichment <- log10(ldsc_results_cts_alt$Enrichment_p)
  fwrite(ldsc_results_cts_alt, "~/data/smontgom/ldsc_cluster_results_cts_alt.txt")
}
# all(as.numeric(stringr::str_remove_all(ldsc_results_conditional$Category, pattern = "Cluster_")) == ldsc_results_conditional$cluster)


#read in ldsc cts results outputted from *old* command
if(file.exists("~/data/smontgom/ldsc_cluster_results_cts_alt_fullyconditional.txt")){
  ldsc_results_cts_alt_fullyconditional <- as.data.frame(fread("~/data/smontgom/ldsc_cluster_results_cts_alt_fullyconditional.txt"))
} else{
  ldsc_output_dir <- "~/repos/ldsc/output/original_command/"
  ldsc_results_cts_alt_fullyconditional_paths <- list.files(ldsc_output_dir)
  ldsc_results_cts_alt_fullyconditional_paths <- paste0(ldsc_output_dir, ldsc_results_cts_alt_fullyconditional_paths[grep(ldsc_results_cts_alt_fullyconditional_paths, pattern = "results")])
  # cluster_names <- paste0("Cluster_", 1:15)
  
  tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
  tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
  tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
  cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
  
  ldsc_results_cts_alt_fullyconditional_paths <- lapply(gwas_names, function(gwas) ldsc_results_cts_alt_fullyconditional_paths[grep(pattern = gwas, ldsc_results_cts_alt_fullyconditional_paths)])
  names(ldsc_results_cts_alt_fullyconditional_paths) <-  gwas_names
  
  ldsc_results_cts_alt_fullyconditional <- as.data.frame(matrix(data = NA, 
                                               ncol = ncol(fread(ldsc_results_cts_alt_fullyconditional_paths[[1]][[1]][1])), 
                                               nrow = length(gwas_names) * length(cluster_names)))
  colnames(ldsc_results_cts_alt_fullyconditional) <- colnames(fread(ldsc_results_cts_alt_fullyconditional_paths[[1]][[1]][1]))
  ldsc_results_cts_alt_fullyconditional$cluster <- rep(cluster_names, length(gwas_names))
  ldsc_results_cts_alt_fullyconditional$gwas <- c(sapply(gwas_names, function(gwas_name) rep(gwas_name, length(cluster_names))))
  for(i in 0:(length(gwas_names)-1) * length(cluster_names) + 1){
    paths <- ldsc_results_cts_alt_fullyconditional_paths[[ldsc_results_cts_alt_fullyconditional$gwas[i]]]
    output <- do.call(rbind, lapply(setNames(cluster_names, cluster_names), function(cli){
      path <- paths[grep(cli, paths)]
      path <- path[grep("baseline-plus", path)]
      temp_output <- fread(path)
      temp_output <- temp_output[temp_output$Category == "L2_0"]
      temp_output$Category <- cli
      return(temp_output)
    }))
    reorder_output <- order(match(output$Category, ldsc_results_cts_alt_fullyconditional[i:(i+length(cluster_names)-1),"cluster"]))
    ldsc_results_cts_alt_fullyconditional[i:(i+length(cluster_names)-1),1:ncol(output)] <- output[reorder_output,]
  }
  ldsc_results_cts_alt_fullyconditional$gwas <- stringr::str_remove(ldsc_results_cts_alt_fullyconditional$gwas, "imputed_")
  ldsc_results_cts_alt_fullyconditional$logPVal_enrichment <- log10(ldsc_results_cts_alt_fullyconditional$Enrichment_p)
  fwrite(ldsc_results_cts_alt_fullyconditional, "~/data/smontgom/ldsc_cluster_results_cts_alt_fullyconditional.txt")
}

#specify graph parameters
max_point_cex <- 3.5
point_cex_power <- 0.35
n_points_for_legend = 7
buffer_min_and_max = 0.05
# minimum_enrichment_logPval <- min(ldsc_results_sub$logPVal_enrichment)
opacity_insig_points <- 0.2
opacity_sig_points <- 0.8

plot_LDSC_comparison = T
reorder_vertical <- T
use_conditional_model = F
use_cts_alt_model = T
use_enrichment = T
partition_by_category <- T
use_heritability = !use_enrichment
fix_axes = F
fix_axes_bounds_enrichment = c(0,5)
fix_axes_bounds_heritability = c(0,1)

#filter by h2 sigma
total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])

# ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_with_satisfactory_heritaility,]
# if(use_conditional_model){
#   ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_with_satisfactory_heritaility,]
# }

trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitname_map <- trait_categories[,c("Tag", "new_Phenotype")]
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
categories <- sort(unique(traitwise_partitions$Category))

#filter by category
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
traits_in_focal_categories <- traitwise_partitions$Tag[traitwise_partitions$Category %in% 
                                                         c("Cardiometabolic", "Aging", "Anthropometric", 
                                                           "Immune", "Psychiatric-neurologic")]

#get final trait list
traits_to_keep <- intersect(traits_with_satisfactory_heritaility, traits_in_focal_categories)
coloc_phenotypes_sub <- traits_to_keep


ldsc_results_sub <- ldsc_results[ldsc_results$gwas %in% traits_to_keep,]

if(use_conditional_model){
  ldsc_results_sub <- ldsc_results_conditional[ldsc_results_conditional$gwas %in% traits_to_keep,]
}

if(use_cts_alt_model){
  ldsc_results_sub <- ldsc_results_cts_alt[ldsc_results_cts_alt$gwas %in% traits_to_keep,]
}

use_fullyconditional_cts_alt_model = T
if(use_fullyconditional_cts_alt_model){
  ldsc_results_sub <- ldsc_results_cts_alt_fullyconditional[ldsc_results_cts_alt_fullyconditional$gwas %in% traits_to_keep,]
}


#BH correction significant
ldsc_results_sub$Enrichment_p_BH <- p.adjust(ldsc_results_sub$Enrichment_p, "BH")
minimum_enrichment_logPval <- log10(min(ldsc_results_sub$Enrichment_p_BH))


bad_boys <- sort(unique(ldsc_results_sub$gwas[ldsc_results_sub$Enrichment < 0 | ldsc_results_sub$Prop._h2 < 0]))
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
    logprob_range <- c(0,5)
    
    max_prob <- diff(logprob_range)
    
    nameloc <- 0.3 + ifelse(change_names | change_names_in_plot, 0, 0.3)
    
    grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/DExEQTL/LDSC_comparison_heritabilities", 
                                           ifelse(use_enrichment, "_enrichment", "_heritability"), 
                                           ifelse(fix_axes, "_fixedAxes", ""), 
                                           ifelse(use_fullyconditional_cts_alt_model, "_fullyConditional", ""), 
                                           # ifelse(use_conditional_model, "_conditionalModel", "_unconditionalModel"), 
                                           ifelse(use_cts_alt_model, "_cell-type-specific-Model", ""), ".pdf"), 
                         width = 1500 / 72, height = 2000 / 72 * 18 / 17 / 114 * length(coloc_phenotypes_sub), family="Arial Unicode MS")
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
      text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                  coloc_phenotypes_sub_newnames)[order_traits], 
           col = bad_boys_cols[order_traits], x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    } else {
      text(paste0(1:length(coloc_phenotypes_sub_newnames), ": ", 
                  coloc_phenotypes_sub_newnames)[order_traits], 
           col = pheno_cols, x = nameloc, y = ylocs, pos = 2, cex = 1, xpd = NA)
    }
    
    #vert axis label
    text(labels = "Phenotypes", x = nameloc, y = 10.25, pos = 2, cex = 2.5, xpd = NA, family="Courier", col = "grey25")
    #horiz axis label
    text(labels = paste0("Heritability ", ifelse(use_enrichment, "Enrichment", "Proportion"), " Across Clusters"), x = 1.2, y = -0.35, cex = 2.25, pos = 1, xpd = NA)
    
    
    #guiding lines for traits
    if(nicole_mods){
      segments(x0 = nameloc, x1 = max_horiz_axis, y0 = ylocs, y1 = ylocs, col = "grey80", lty = 3, lwd = 0.5)
    }
    
    #axes
    segments(x0 = nameloc, x1 = nameloc, y0 = 10.35, y1 = -0.1, lwd = 2)
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
    segments(x0 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             x1 = (1 - logprob_range[1]) / max_prob * (2 - nameloc - 0.025) + nameloc, 
             y0 = 10, y1 = -0.075, lwd = 2, lty = 1, col = "grey75")
    
    
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
    lower_legend_by <- 0.6
    rect(xleft = 1.97, ybottom = 8.15 - lower_legend_by - 0.05*n_points_for_legend, 
         ytop = 10.1 - lower_legend_by, xright = 2.2, border = NA, col = "white")
    # text(labels = sapply(names(deg_eqtl_list), function(ts) stringr::str_to_title(paste0(strsplit(ts, "-")[[1]][-1], collapse = " "))),
    #      y = seq(9.95,8.75,length.out = length(names(deg_eqtl_list))), x = 2, pos = 4)
    # points(x = rep(1.9925, length(names(deg_eqtl_list))), y = seq(9.955,8.755,length.out = length(names(deg_eqtl_list))), col = cols$Tissue, pch = 15, cex = 1.75)
    # text(labels = stringr::str_replace_all(cluster_names, "_", " "),
    #      y = seq(9.95,8.75,length.out = length(cluster_names)), x = 2, pos = 4)
    text(labels = sapply(cluster_names, function(cli) strsplit(cli, "-")[[1]][1]),
         y = seq(9.95,8.75 - lower_legend_by,length.out = length(cluster_names)), x = 2, pos = 4)
    points(x = rep(1.9925, length(cluster_names)), y = seq(9.955,8.755 - lower_legend_by,length.out = length(cluster_names)), col = cols$cluster, pch = 19, cex = 1.75)
    
    #plot legend for points
    lower_legend_by <- lower_legend_by + 0.2
    pt_loc_expander <- 4
    point_legend_cexes <- seq(from = 0.4, to = max_point_cex, length.out = n_points_for_legend)
    points_legend_pchs <- rep(19, n_points_for_legend)
    point_legend_cexes_log10_pvals <- round((point_legend_cexes / max_point_cex)^(1/point_cex_power) * minimum_enrichment_logPval, 2)
    # points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
    points_legend_pchs[point_legend_cexes_log10_pvals < (log10(0.05))] <- 18
    
    points(y = 8.55 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
           x = rep(1.9925, n_points_for_legend), cex = point_legend_cexes, col = "grey50", pch = points_legend_pchs)
    text(labels = latex2exp::TeX("log_{10}(enrichment p-val)"), y = 8.6 - lower_legend_by, x = 1.975, cex = 1.1, pos = 4,  font = 2)
    # text(labels = paste0("0.05 FDR @ ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n           (Bonferroni)"), y = 8.45 , x = 2.1, cex = 0.8, pos = 4,  font = 2)
    text(labels = latex2exp::TeX(paste0("0.05 FDR @ $\\alpha = 0.05$")), 
         y = 8.45 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
    text(labels = latex2exp::TeX(paste0("(Benjamini-Hochberg)")), 
         y = 8.35 - lower_legend_by, x = 2.07, cex = 0.8, pos = 4,  font = 2)
    # text(labels = paste0("< ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2), "\n> ", round(log10(0.05) - log10(nrow(ldsc_results_sub)),2)), 
    #      y = 8.25 , x = 2.175, cex = 1, pos = 4,  font = 2)
    text(labels = paste0("< ", round(log10(0.05),2), "\n> ", round(log10(0.05),2)), 
         y = 8.25 - lower_legend_by - 0.1, x = 2.145, cex = 1, pos = 4,  font = 2)
    points(pch = c(19,18), y = c(8.325, 8.225) - lower_legend_by - 0.1, x = rep(2.145,2), cex = c(1.25, 1.75), 
           col = sapply(c(opacity_insig_points, opacity_sig_points), function(opcty) adjustcolor("grey50", alpha.f = opcty)))
    
    
    text(labels = point_legend_cexes_log10_pvals,y = 8.545 - lower_legend_by - cumsum(point_legend_cexes + pt_loc_expander) / sum(point_legend_cexes) * 0.05*n_points_for_legend, 
           x = rep(1.99, n_points_for_legend) + point_legend_cexes / 200, cex = 1, pos = 4, pch = 19)
    
    
    #figure out cex params
    for(cluster in cluster_names){
        
        # horizontal lines for colocalizing traits
        trait_locs <- ylocs
        if(use_enrichment){
          trait_probs <- ldsc_results_sub$Enrichment[ldsc_results_sub$cluster == cluster]
          #for bonferroni
          # trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
          #for BH adjustment
          trait_logPvals <- log10((ldsc_results_sub$Enrichment_p_BH[ldsc_results_sub$cluster == cluster])[order_traits])
        } else {
          trait_probs <- ldsc_results_sub$Prop._h2[ldsc_results_sub$cluster == cluster]
          trait_logPvals <- (ldsc_results_sub$logPVal_enrichment[ldsc_results_sub$cluster == cluster])[order_traits]
        }
        trait_probs <- trait_probs[order_traits]
        
        point_cex <- (-trait_logPvals / -minimum_enrichment_logPval)^point_cex_power
        point_cex <- point_cex * max_point_cex
        
        good_points <- trait_probs >= logprob_range[1] & trait_probs <= logprob_range[2]
        
        pchs <- rep(19, length(trait_probs))
        opacities <- rep(opacity_insig_points, length(trait_probs))
        
        #mark "significant" points
        # pchs[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- 18
        # opacities[trait_logPvals < (log10(0.05) - log10(nrow(ldsc_results_sub)))] <- opacity_sig_points
        pchs[trait_logPvals < log10(0.05)] <- 18
        opacities[trait_logPvals < log10(0.05)] <- opacity_sig_points
        
        point_cols <- sapply(opacities, function(opcty) grDevices::adjustcolor(cols$Tissue[match(cluster, cluster_names)], alpha.f = opcty))
        
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

plot(ldsc_results$Prop._h2,  ldsc_results_cts_alt$Prop._h2, main = "Proportion Heritability", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))

plot(ldsc_results$Enrichment,  ldsc_results_cts_alt$Enrichment, main = "Proportion Heritability", xpd = NA,
     xlab = latex2exp::TeX("Conditional Model Prop. h^2"), ylab =  latex2exp::TeX("Uncnditional Model Prop. h^2"), pch = 19, col = adjustcolor(1, 0.5))

ldsc_results_conditional
ldsc_results_cts_alt

#### cluster-based marginal enrichment ####

trait_type_specific = T
grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/DExEQTL/LDSC_cluster-specific_enrichments", ifelse(trait_type_specific, "_categorized-by-categories"),".pdf"), 
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

#snag gwas log files
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
total_h2_sigma_thresh <- 7
traits_with_satisfactory_heritaility <- unique(log_files$gwas[log_files$h2 / log_files$h2se > total_h2_sigma_thresh])



#take 2
new_clusters <- c("t53-cortex",  "t56-vastus-lateralis",  "t58-heart",  "t61-colon",  "t62-spleen",  "t64-ovaries",  "t65-aorta",  "t66-lung",  "t67-small-intestine",  "t69-brown-adipose",  "t70-white-adipose",  "t30-blood-rna",  "t52-hippocampus",  "t54-hypothalamus",  "t63-testes",  "t55-gastrocnemius",  "t60-adrenal",  "t59-kidney",  "t68-liver",  "1w",  "2w",  "4w",  "8w",  "female-1w",  "male-2w",  "female-4w",  "male-8w",  "all_DE")

#take 3
geneset_info <- expand.grid(c(7,15), c("t55-gastrocnemius", "t68-liver", "t58-heart"))
new_clusters <- paste0("Cluster_", trimws(apply(geneset_info, 1, paste0, collapse = "-")), "")

#take 4
tissue_abbrev <- unique(MotrpacBicQC::bic_animal_tissue_code$abbreviation)
tissue_abbrev <- tissue_abbrev[!is.na(tissue_abbrev)]
tissue_abbrev <- setdiff(tissue_abbrev, c("PLASMA", "HYPOTH", "TESTES", "OVARY", "VENACV"))
cluster_names <- paste0(tissue_abbrev, "-sex_homogeneous_changing")
new_clusters <- sort(paste0(tissue_abbrev, "-sex_homogeneous_changing"))

output_files <- c(sapply(c("_Cluster_Addenda.cell_type_results.txt", "_Cluster_1-15.cell_type_results.txt"), 
                         function(x) paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), x)))
output_files <- paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), c("_Cluster_1-15.cell_type_results.txt"))
output_files <- paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), c("_Cluster_Addenda.cell_type_results.txt"))
output_files <- paste0(stringr::str_remove_all(gwas_summary_files, ".txt.gz"), "_cluster_", 
                       paste0(new_clusters[1], "-", new_clusters[length(new_clusters)]), 
                       ".cell_type_results.txt")

ldsc_results <- do.call(rbind, lapply(1:length(output_files), function(i) 
  cbind(fread(paste0(ldsc_results_dir, output_files[i])),  trait = strsplit(output_files[i], "_Cluster")[[1]][1])))

ldsc_results$twoTailedPVal <- ldsc_results$Coefficient_P_value
ldsc_results$twoTailedPVal[ldsc_results$twoTailedPVal > 0.5] <- 1 - ldsc_results$twoTailedPVal[ldsc_results$twoTailedPVal > 0.5]
ldsc_results$twoTailedPVal <- ldsc_results$twoTailedPVal * 2

#filter by category
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
ldsc_results$trait_name <- gsub(ldsc_results$trait, pattern = "imputed_", replacement = "")
ldsc_results$trait_name <- gsub(ldsc_results$trait_name, pattern = "_cluster.*", replacement = "")

ldsc_results$trait_category <- traitwise_partitions$Category[match(ldsc_results$trait_name, traitwise_partitions$Tag)]
ldsc_results <- ldsc_results[ldsc_results$trait_category %in% 
                               c("Cardiometabolic", "Aging", "Anthropometric", "Immune", "Psychiatric-neurologic")]

#filter by heritability threshold
ldsc_results <- ldsc_results[ldsc_results$trait_name %in% traits_with_satisfactory_heritaility,]


hist(log(ldsc_results$Coefficient_P_value), xlab = "log 1-tailed p-value", main = "", breaks = -200:0/10)
abline(v = log(0.05 / nrow(ldsc_results)), col = 2, lwd = 4)
text(x = log(0.05 / nrow(ldsc_results)), y = 200, 
     labels = latex2exp::TeX("Bonferroni-adjusted $\\alpha$"), srt = 90, pos = 2)
hist(ldsc_results$Coefficient_P_value)
hist(log(ldsc_results$twoTailedPVal), breaks = -200:0/10)
abline(v = log(0.05 / nrow(ldsc_results)), col = 2, lwd = 4)
hist((ldsc_results$twoTailedPVal))
qvalue::pi0est(ldsc_results$Coefficient_P_value)
qvalue::pi0est(ldsc_results$twoTailedPVal)
sum(log(ldsc_results$Coefficient_P_value) < log(0.05 / nrow(ldsc_results)))
head(ldsc_results[order(ldsc_results$Coefficient_P_value, decreasing = F),], 20)
head(ldsc_results[order(ldsc_results$twoTailedPVal, decreasing = F),], 20)
head(ldsc_results[order(ldsc_results$Coefficient, decreasing = T),], 20)

#### TWAS results ####
# twas_tissues <- fread(file = "~/repos/fusion_twas-master/output/all_results.txt")
# metaxscan_results <- fread(file = "~/repos/MetaXcan/software/results/all_results.txt")
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
# load('~/data/smontgom/coloc_list_pp4threshold_1e-04.RData')
# geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")

#whoops there is redundancy here
# tiss_gene_twas <- paste0(twas_tissues$tissue, "_", twas_tissues$ID) 
# tiss_gene_coloc <- paste0(coloc_list$motrpac_tissue, "_", coloc_list$gene_id) 
# 
# for(i in 1:length(deg_eqtl_list)){
#   cat(paste0(i, " "))
#   temp <- deg_eqtl_list[[i]]
#   tiss_gene_twas_equiv <- paste0(motrpac_gtex_map[match(names(deg_eqtl_list)[i], names(motrpac_gtex_map))], "_", temp$gene_name)
#   tiss_gene_coloc_equiv <- paste0(names(deg_eqtl_list)[i], "_", temp$gene_id)
#   temp$TWAS_PVAL <- twas_tissues$TWAS.P[match(tiss_gene_twas_equiv, tiss_gene_twas)]
#   temp$colocPP4 <- coloc_list$p4[match(tiss_gene_coloc_equiv, tiss_gene_coloc)]
#   deg_eqtl_list[[i]] <- temp
# }
# 
# twas_pvals <- (do.call(rbind,deg_eqtl_list)$TWAS_PVAL)
# coloc_pp4s <- (do.call(rbind,deg_eqtl_list)$colocPP4)
# plot(-log(twas_pvals), -log(1-coloc_pp4s))
# 
# all_gene_names <- unique(map$human_gene_symbol)
# gene_names_DE <- unique(do.call(rbind, deg_eqtl_list)$gene_name)
# mean(!is.na(match(gene_names_DE, all_gene_names)))
# gene_names_TWAS <- unique(twas_tissues$ID)
# gene_names_metaxscan <- unique(metaxscan_results$gene_name)
# mean(!is.na(match(gene_names_TWAS, all_gene_names)))
# length(setdiff(gene_names_DE, gene_names_TWAS)) / length(gene_names_DE)
# sapply(names(deg_eqtl_list), function(i) mean(!is.na(deg_eqtl_list[[i]]$TWAS_PVAL))*100)
# twas_tissuegenes <- sapply(unique(twas_tissues$tissue), function(tiss) unique(twas_tissues$ID[twas_tissues$tissue == tiss]))
# sapply(twas_tissuegenes, function(x) length(x))
# length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
# length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
# gene_IDs_metaxscan <- unique(metaxscan_results$gene)
# gene_IDs_DE <- unique(do.call(rbind, deg_eqtl_list)$gene_id)
# length(setdiff(gene_names_DE, gene_names_metaxscan)) / length(gene_names_DE)
# length(setdiff(gene_IDs_DE, gene_IDs_metaxscan)) / length(gene_IDs_DE)

#now let's use the TWAS results from Barbeira et al. 2021
twas_results_directory <- "~/data/smontgom/eqtl/"
# tissue = motrpac_gtex_map[6]
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
gwas_names <- stringr::str_replace_all(gwas_summary_files, ".txt.gz", "")
gwas_names <- stringr::str_replace_all(gwas_names, "imputed_", "")
# trait <- gwas_names[grep("hyperte", gwas_names, T)]
# twas <- fread(file = paste0(twas_results_directory, "spredixcan_igwas_gtexmashrv8_", trait, "__PM__", tissue, ".csv"))

# jac <- function(x1, x2){length(intersect(x1, x2)) / length(union(x1, x2))}
# jac(unique(deg_eqtl_list[[names(tissue)]]$gene_name), twas$gene_name)
# jac(unique(deg_eqtl_list[[names(tissue)]]$gene_id), twas$gene)
# 1 - (length(setdiff(unique(deg_eqtl_list[[names(tissue)]]$gene_id[deg_eqtl_list[[names(tissue)]]$selection_fdr < 0.05]), twas$gene)) / 
#        length(unique(deg_eqtl_list[[names(tissue)]]$gene_id[deg_eqtl_list[[names(tissue)]]$selection_fdr < 0.05])))

available_twas_traits <- list.files(path = twas_results_directory)
available_twas_tissues <- table(sapply(available_twas_traits, function(x) strsplit(x, "__PM__")[[1]][2]))
available_twas_traits <- table(sapply(gsub(x = available_twas_traits, "spredixcan_igwas_gtexmashrv8_", ""), 
                                      function(x) strsplit(x, "__P")[[1]][1]))

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

# data <- data.frame(trait = names(sort(sapply(gwas_names, function(trait) 1-qvalue::pi0est(all_twas[[trait]]$pvalue)$pi0))),
#                    pi1 = sort(sapply(gwas_names, function(trait) 1-qvalue::pi0est(all_twas[[trait]]$pvalue)$pi0)))
# ggplot(data, aes(x=trait, y=pi1)) + 
#   geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
#   coord_flip() + scale_x_discrete(limits=data$trait)

all.twas <- do.call(rbind, all_twas)
# 1-qvalue::pi0est(all.twas$pvalue)$pi0
# hist(all.twas$pvalue)

#filter by category
trait_categories <- read.csv("~/data/smontgom/gwas_metadata.csv", header = T)
traitwise_partitions <- trait_categories[,c("Tag", "Category")]
all.twas$trait_category <- traitwise_partitions$Category[match(all.twas$trait, traitwise_partitions$Tag)]
salient.categories <- c("Cardiometabolic", "Aging", "Anthropometric", 
                        "Immune", "Psychiatric-neurologic")
some.twas <- all.twas[all.twas$trait_category %in% salient.categories]

#cleanup memory
rm("all.twas")
rm("all_twas")

hist(some.twas$pvalue)
# hist(log(all.twas$pvalue), breaks = -700:1)
# 1-qvalue::pi0est(some.twas$pvalue)$pi0
# sum(some.twas$pvalue < 0.05 / nrow(some.twas))

# data <- data.frame(Trait = names(sort(table(some.twas[p.adjust(some.twas$pvalue, method = "bonf") < 0.05,]$trait))),
#                    n_sig_bonferroni = as.integer(sort(table(some.twas[p.adjust(some.twas$pvalue, method = "bonf") < 0.05,]$trait))))
# ggplot(data, aes(x=Trait, y=n_sig_bonferroni)) + 
#   geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
#   coord_flip() + scale_x_discrete(limits=data$Trait)                      
# 
# data <- data.frame(trait = names(sort(table(all.twas[all.twas$pvalue < 0.05 / nrow(all.twas),]$trait))),
#                    n_sig = as.integer(sort(table(all.twas[all.twas$pvalue < 0.05 / nrow(all.twas),]$trait))))
# ggplot(data, aes(x=trait, y=n_sig)) + 
#   geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
#   coord_flip() + scale_x_discrete(limits=data$trait)                      
# 
# data <- data.frame(trait = names(sort(table(some.twas[some.twas$pvalue < 0.05 / nrow(some.twas),]$gene_name))),
#                    n_sig = as.integer(sort(table(some.twas[some.twas$pvalue < 0.05 / nrow(some.twas),]$gene_name))))
# ggplot(data, aes(x=trait, y=n_sig)) + 
#   geom_bar(stat = "identity") + theme_bw(base_size = 10) + 
#   coord_flip() + scale_x_discrete(limits=data$trait)    

#ok now let's try to x-reference it to the differential expression results!
# DEsub <- deg_eqtl_list[[tissue]]
# hist(DEsub$selection_fdr)
# sum(p.adjust(some.twas$pvalue, method = "BH") < 0.05) / length(some.twas$pvalue)
# some.twas$adj_pvalue <- p.adjust(some.twas$pvalue, method = "BH")

if(!file.exists("~/data/smontgom/ihw_results_some.twas.RData")){
  ihw_twas <- as.data.frame(cbind(trait_tissue = as.factor(paste0(some.twas$trait, "~", some.twas$tissue)), pvalue = some.twas$pvalue))
  ihw_results <- IHW::ihw(pvalue ~ trait_tissue, data = ihw_twas, alpha = 0.05)  
  save(ihw_results, file = "~/data/smontgom/ihw_results_some.twas.RData")
} else {
  load("~/data/smontgom/ihw_results_some.twas.RData")
  some.twas$adj_pvalue <- IHW::adj_pvalues(ihw_results)
  rm(ihw_results)
}

salient_twas <- unique(some.twas$trait)
# deg_eqtl_list_TWAS <- deg_eqtl_list
# for(tissue_i in names(motrpac_gtex_map)){
# 
#   print(tissue_i)
#   # twas_addition <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas)))
#   # twas_addition_nominalPValue <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas)))
#   # twas_addition_zscore <- data.frame(matrix(nrow = nrow(deg_eqtl_list_TWAS[[tissue_i]]), ncol = length(salient_twas)))
#   # colnames(twas_addition) <- colnames(twas_addition_zscore) <- salient_twas
#   twas_sub <- some.twas[tissue == tissue_i]
# 
#   gtex_motrpac <- rna_dea$timewise_dea[tissue == tissue_i]
#   for(trait_j in salient_twas){
#     # twas_sub_sub <- twas_sub[trait == trait_j]
#     # twas_addition[,trait_j] <- twas_sub_sub$adj_pvalue[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
# 
#     # subset in TWAS
#     twas_sub_sub <- twas_sub[trait == trait_j]
#     twas_sub_sub[, human_ensembl_gene := gsub('\\..*','',gene)]
#     twas_sub_sub$human_gene_symbol <- twas_sub_sub$gene_name
#     p_twas_in_map <- round(mean(unique(twas_sub_sub$human_gene_symbol) %in% unique(map$human_gene_symbol)), 3)
# 
#     # match human genes with rat ensembl genes
#     twas_sub_sub = merge(twas_sub_sub, map, by='human_gene_symbol')
#     twas_sub_sub <- twas_sub_sub[,c("feature_ID", "gene_name", "pvalue", "zscore")]
#     colnames(twas_sub_sub)[-match(c("feature_ID", "gene_name"), colnames(twas_sub_sub))] <- paste0(trait_j, ".",
#                         colnames(twas_sub_sub)[-match(c("feature_ID", "gene_name"), colnames(twas_sub_sub))])
# 
#     matched_twas <- twas_sub_sub[match(gtex_motrpac$feature_ID, twas_sub_sub$feature_ID),-c("feature_ID")]
#     if(trait_j != salient_twas[1]){matched_twas <- matched_twas[,-"gene_name"]}
#     gtex_motrpac <- cbind(gtex_motrpac, matched_twas)
# 
#     cat(paste0("\nPr(twas in map): ", p_twas_in_map,
#                  "; Pr(feature_IDs in GTEx): ", round(mean(!is.na(matched_twas[,1])), 3)))
# 
#     # twas_addition_nominalPValue[,trait_j] <- twas_sub_sub$pvalue[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
#     # twas_addition_zscore[,trait_j] <- twas_sub_sub$zscore[match(deg_eqtl_list_TWAS[[tissue_i]]$gene_id, twas_sub_sub$gene)]
#   }
#   deg_eqtl_list_TWAS[[tissue_i]] <- gtex_motrpac
# 
#   # colnames(twas_addition) <- paste0(colnames(twas_addition), "_BH_PValue")
#   # colnames(twas_addition_nominalPValue) <- paste0(colnames(twas_addition_nominalPValue), "_Nominal_PValue")
#   # colnames(twas_addition_zscore) <- paste0(colnames(twas_addition_zscore), "_zscore")
#   # deg_eqtl_list_TWAS[[tissue_i]] <- cbind(deg_eqtl_list_TWAS[[tissue_i]], twas_addition, twas_addition_nominalPValue, twas_addition_zscore)
# 
# }

#get number of p-values
# deg_eqtl_list_TWAS_all <- do.call(rbind, deg_eqtl_list_TWAS)
# deg_eqtl_list_TWAS_all <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$sex == "female"]
# deg_eqtl_list_TWAS_all <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$comparison_group == "8w"]
# 
# all_twas_pvalues <- unlist(as.data.frame(deg_eqtl_list_TWAS_all)[,grep(colnames(deg_eqtl_list_TWAS_all), pattern = ".pvalue")])
# all_twas_pvalues <- all_twas_pvalues[!is.na(all_twas_pvalues)]
# n_twas_comparisons <- length(all_twas_pvalues)

#first do simple gene intersect
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")

if(!exists("cluster_membership") | !exists("node_metadata")){
  # load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
  # cluster_membership <- cluster_membership[cluster_membership$cluster %in% c(1,3,7,15) & cluster_membership$assay == "TRNSCRPT",]
  # map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
  # cluster_genes <- map$human_gene_symbol[match(cluster_membership$feature_ID, map$feature_ID)]
  # cluster_genes <- cluster_genes[!is.na(cluster_genes)]
  
  load("~/data/smontgom/graphical_analysis_results_20211220.RData")
  nodes_to_look_at_list <- list(c("1w_F1_M1", "1w_F-1_M-1"),
                             c("2w_F1_M1", "2w_F-1_M-1"),
                             c("4w_F1_M1", "4w_F-1_M-1"),
                             c("8w_F1_M1", "8w_F-1_M-1"))
  
  node_metadata_list <- lapply(setNames(nodes_to_look_at_list, paste0(2^(0:3), "w")), function(nodes_to_look_at){
    node_metadata <- lapply(setNames(nodes_to_look_at, nodes_to_look_at), function(node_to_look_at)
      cbind(do.call(rbind, strsplit(node_sets[node_to_look_at][[node_to_look_at]][grep(
        node_sets[node_to_look_at][[node_to_look_at]], pattern = "TRNSCRPT")], ";")), node_to_look_at))
    node_metadata <- as.data.table(do.call(rbind, node_metadata))
    colnames(node_metadata) <- c("ome","tissue","ensembl_gene", "node")
    
    
    node_metadata$rat_gene_symbol <- gene_map$RAT_SYMBOL[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_gene_symbol <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_ensembl_gene <- gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(node_metadata$ensembl_gene, gene_map$RAT_ENSEMBL_ID)]
    node_metadata$human_ensembl_gene <- gsub(node_metadata$human_ensembl_gene, pattern = "\\..*", replacement = "")
    node_metadata$cluster <- paste0(node_metadata$tissue, "-", node_metadata$node)
    node_metadata$cluster <- paste0(node_metadata$tissue, "-", "sex_homogeneous_changing")
    node_metadata
  })
  
}

# deg_eqtl_list_TWAS_all_1.3.7.15 <- deg_eqtl_list_TWAS_all[deg_eqtl_list_TWAS_all$gene_name %in% cluster_genes,]
# all_twas_pvalues_clusters_1.3.7.15 <- unlist(as.data.frame(deg_eqtl_list_TWAS_all_1.3.7.15)[,grep(colnames(deg_eqtl_list_TWAS_all_1.3.7.15_subsetSexTime), pattern = ".pvalue")])
# all_twas_pvalues_clusters_1.3.7.15 <- all_twas_pvalues_clusters_1.3.7.15[!is.na(all_twas_pvalues_clusters_1.3.7.15)]
# n_twas_comparisons_clusters_1.3.7.15 <- length(all_twas_pvalues_clusters_1.3.7.15)

#hmm wait that is not quite right... let's subset it to the cluster-tissue subset first, and then further subset it to the FDR genes
# deg_eqtl_list_TWAS_cluster_subset_list <- lapply(setNames(node_metadata_list, paste0(2^(0:3), "w")), function(node_metadata){
#   lapply(deg_eqtl_list_TWAS[-which(names(deg_eqtl_list_TWAS) %in% c("t1000-gonads", "t63-testes", "t64-ovaries"))], function(tissue_dea) 
#     tissue_dea[tissue_dea$feature_ID %in% node_metadata$ensembl_gene[node_metadata$tissue == tissue_dea$tissue_abbreviation[1]],]) 
# })

# deg_eqtl_list_TWAS_cluster_subset <- deg_eqtl_list_TWAS_cluster_subset_list[[2]]
# 
# deg_eqtl_list_TWAS_cluster_subset_ntests <- do.call(rbind, deg_eqtl_list_TWAS_cluster_subset)
# deg_eqtl_list_TWAS_cluster_subset_ntests <- deg_eqtl_list_TWAS_cluster_subset_ntests[deg_eqtl_list_TWAS_cluster_subset_ntests$sex == "female" &
#                                                                                        deg_eqtl_list_TWAS_cluster_subset_ntests$comparison_group == "8w",]
# twas_pvalues_clusters <- unlist(as.data.frame(deg_eqtl_list_TWAS_cluster_subset_ntests)[,grep(colnames(deg_eqtl_list_TWAS_cluster_subset_ntests), pattern = ".pvalue")])
# twas_pvalues_clusters <- twas_pvalues_clusters[!is.na(twas_pvalues_clusters)]
# n_twas_comparisons_clusters <- length(twas_pvalues_clusters)

n_deg_sigtwas_intersect <- as.data.frame(matrix(0, nrow = length(names(motrpac_gtex_map)), 
                                                   ncol = length(salient_twas), 
                                                   dimnames = list(names(motrpac_gtex_map), salient_twas)))
timepoint = "8w"
adj_pvalue_alpha <- 0.05
sig_twas_by_tissue <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), 
                         function(tissue) some.twas[some.twas$tissue == tissue & some.twas$adj_pvalue < adj_pvalue_alpha,])
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
for(tissue_i in names(motrpac_gtex_map)){
  
  print(tissue_i)
  # DELT <- as.data.frame(deg_eqtl_list_TWAS_cluster_subset[[tissue_i]])
  # 
  # DE_inds <- which(DELT$selection_fdr <= 1) #since we're subsetting to monotonic clusters already, but why not
  # TWAS_inds <- apply(log(DELT[,grep(colnames(DELT), pattern = ".pvalue")]) <= (log(0.05) - log(n_twas_comparisons_clusters)), 2, which)
  
  # intersect_inds <- lapply(TWAS_inds, function(twi) intersect(DE_inds, twi))
  # intersect_genes <- lapply(intersect_inds, function(ii) unique(DELT$gene_name[ii]))
  # intersect_genes <- lapply(intersect_genes, function(ii) intersect(ii, cluster_genes)) #subset to just monotonic sex-homogenous clusters
  

  sig_twas_by_trait <- lapply(setNames(salient_twas, salient_twas), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(salient_twas, salient_twas), 
                              function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  
  
  n_intersect <- sapply(sig_twas_by_trait_genes, function(twas_genes) 
    length(intersect(twas_genes, node_metadata_list[[timepoint]]$human_gene_symbol[
      node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]])))
  
  
  
  ##WAS HERE
  # names(n_intersect) <- gsub(names(n_intersect), pattern = ".pvalue", replacement = "")
  n_deg_sigtwas_intersect[tissue_i,names(n_intersect)] <- n_intersect

}

#tidy this up
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,order(apply(n_deg_sigtwas_intersect, 2, sum), decreasing = T)]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[,apply(n_deg_sigtwas_intersect, 2, sum) != 0]
rownames(n_deg_sigtwas_intersect) <- MotrpacBicQC::tissue_abbr[rownames(n_deg_sigtwas_intersect)]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[order(match(rownames(n_deg_sigtwas_intersect), MotrpacBicQC::tissue_order)),]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[nrow(n_deg_sigtwas_intersect):1,]
n_deg_sigtwas_intersect <- n_deg_sigtwas_intersect[!(rownames(n_deg_sigtwas_intersect) %in% c("OVARY", "TESTES")),]

#add in an extra row and column for totals
#needs to be ones that can be mapped to rat orthologs
#colsums
load("~/data/smontgom/genes_tested_in_transcriptome_DEA.RData")
all_orthologs_tested <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(genes_tested_in_transcriptome_DEA, gene_map$RAT_ENSEMBL_ID)]
all_orthologs_tested <- all_orthologs_tested[!is.na(all_orthologs_tested)]
sig_twas_by_trait_genes <- sapply(setNames(salient_twas, salient_twas), 
                                  function(trait_i) length(intersect(all_orthologs_tested, sig_twas_by_trait[[trait_i]]$gene_name)))
sig_twas_by_trait_genes <- sig_twas_by_trait_genes[colnames(n_deg_sigtwas_intersect)]
prop_twas_are_degs <- t(sapply(rownames(n_deg_sigtwas_intersect), function(tissue) n_deg_sigtwas_intersect[tissue,] / sig_twas_by_trait_genes))
prop_twas_are_degs <- apply(prop_twas_are_degs, 2, unlist)
prop_twas_are_degs <- prop_twas_are_degs[,order(apply(prop_twas_are_degs, 2, mean), decreasing = T)]
# hist(apply(prop_twas_are_degs, 2, mean))
#rowsums
all_twas_genes_tested <- unique(some.twas$gene_name)
n_genes_in_nodes <- sapply(setNames(MotrpacBicQC::tissue_abbr[names(motrpac_gtex_map)], MotrpacBicQC::tissue_abbr[names(motrpac_gtex_map)]), function(tissue) 
  length(intersect(all_twas_genes_tested, node_metadata_list[[timepoint]]$human_gene_symbol[node_metadata_list[[timepoint]]$tissue == tissue])))
n_genes_in_nodes <- n_genes_in_nodes[rownames(n_deg_sigtwas_intersect)]

prop_degs_are_twas <- (sapply(colnames(n_deg_sigtwas_intersect), function(trait) n_deg_sigtwas_intersect[,trait] / n_genes_in_nodes[rownames(n_deg_sigtwas_intersect)]))
prop_degs_are_twas <- apply(prop_degs_are_twas, 2, unlist)
prop_degs_are_twas <- prop_degs_are_twas[,order(apply(prop_degs_are_twas, 2, mean), decreasing = T)]

#### obtain a bayesian estimate of the twas proportion and perform an 'enrichment' analysis ####

data1 <- data.frame(count = as.integer(unlist(c(n_deg_sigtwas_intersect))))
data1$tissue <- rep(rownames(n_deg_sigtwas_intersect), ncol(n_deg_sigtwas_intersect))
data1$trait <- unlist(lapply(colnames(n_deg_sigtwas_intersect), function(tri) rep(tri, nrow(n_deg_sigtwas_intersect))))
data1$total <- n_genes_in_nodes[data1$tissue]
# data$total <- unlist(lapply(sig_twas_by_trait_genes, function(trait_total) rep(trait_total, nrow(n_deg_sigtwas_intersect))))
data1$TWAS_Hit <- "YES"

total_number_of_possible_hits <- length(intersect(all_orthologs_tested, all_twas_genes_tested))
data2 <- data1
data2$TWAS_Hit <- "NO"
# data2$count <- n_genes_in_nodes[data$tissue] - data$count
# data2$total <- unlist(lapply(total_number_of_possible_hits - sig_twas_by_trait_genes, function(trait_total) rep(trait_total, nrow(n_deg_sigtwas_intersect))))
data2$count <- sig_twas_by_trait_genes[data1$trait] - data1$count
data2$total <- (total_number_of_possible_hits - n_genes_in_nodes)[data1$tissue]

# plot((data2$count / data2$total)[data$trait == "UKB_50_Standing_height"], 
#      (data1$count / data1$total)[data$trait == "UKB_50_Standing_height"])
par(mfrow = c(3,5))
for(i in (1:15)[-4]){
  plot(logit((data2$count / data2$total)[data2$tissue == names(n_genes_in_nodes)[i]]), 
       logit((data1$count / data1$total)[data1$tissue == names(n_genes_in_nodes)[i]]),
       main = names(n_genes_in_nodes)[i], 
       col = adjustcolor(1, 0.5),
       pch = 19, cex = 2)
  abline(0,1)
}


# for(i in 1:30){
#   plot(logit((data2$count / data2$total)[data$trait == names(sig_twas_by_trait_genes)[i]]), 
#        logit((data$count / data$total)[data$trait == names(sig_twas_by_trait_genes)[i]]),
#        main = names(sig_twas_by_trait_genes)[i], 
#        col = tissue_cols[(data$tissue)[data$trait == names(sig_twas_by_trait_genes)[i]]],
#        pch = 19, cex = 2)
#   print(round(mean((data$count / data$total)[data$trait == names(sig_twas_by_trait_genes)[i]]),4))
#   print(round(mean((data2$count / data2$total)[data$trait == names(sig_twas_by_trait_genes)[i]]),4))
#   abline(0,-1)
# }
#check residuals

# data <- rbind(data1, data2)
# data <- data[!(data$tissue %in% unique(data1$tissue[data1$total == 0])),]
# 
# 
# 
# #average enrichment? quick check
# (mean((data$count / data$total)[data$TWAS_Hit == "YES"]) - 
#     mean((data$count / data$total)[data$TWAS_Hit == "NO"])) / 
#   mean((data$count / data$total)[data$TWAS_Hit == "NO"])
# 
# mean((((data$count / data$total)[data$TWAS_Hit == "YES"]) - 
#     ((data$count / data$total)[data$TWAS_Hit == "NO"])) > 0)
# 
# 
# d <- list(count = data$count,
#           total = data$total,
#           trait = as.integer(as.factor(data$trait)),
#           tissue = as.integer(as.factor(data$tissue)),
#           TWAS_hit = as.integer(data$TWAS_Hit == "YES") + 1,
#           n = nrow(data),
#           n_trait = ncol(n_deg_sigtwas_intersect),
#           n_tissue = nrow(n_deg_sigtwas_intersect),
#           n_cond = 2
#           )
# 
# hist(sapply(1:80, function(i) mean((d$count / d$total)[d$TWAS_hit == 2 & d$trait == i])), breaks = 0:500/1000)
# sapply(1:20, function(i) mean((d$count / d$total)[d$TWAS_hit == 2 & d$trait == i]))
# hist(apply(n_deg_sigtwas_intersect, 2, mean) / sig_twas_by_trait_genes, breaks = 0:150/1000)
# 
# hist(sapply(colnames(n_deg_sigtwas_intersect), function(i) mean((data2$count / data2$total)[data2$TWAS_Hit == "NO" & data2$trait == i])), 
#      breaks = 0:50/100)
# hist(sapply(colnames(n_deg_sigtwas_intersect), function(i) mean((data$count / data$total)[data$TWAS_Hit == "YES" & data$trait == i])), 
#      breaks = 0:50/100)
# 
# sapply(1:80, function(i) ((d$count / d$total)[d$TWAS_hit == 1 & d$trait == i]))
# sapply(1:80, function(i) ((d$tissue)[d$TWAS_hit == 1 & d$trait == i]))
# sapply(1:80, function(i) ((d$trait)[d$TWAS_hit == 1 & d$trait == i]))
# sapply(1:80, function(i) ((d$total)[d$TWAS_hit == 2 & d$trait == i]))
# sapply(1:80, function(i) ((d$count)[d$TWAS_hit == 2 & d$trait == i]))

# 
# 
# #subset to some tissues & traits for troubleshooting
# trait_subset <- names(sig_twas_by_trait_genes)[order(sig_twas_by_trait_genes, decreasing = T)][1:20]
# tissue_subset <- names(n_genes_in_nodes)[order(n_genes_in_nodes, decreasing = T)][1:5]
# d <- list(count = data$count[data$tissue %in% tissue_subset & data$trait %in% trait_subset],
#           total = data$total[data$tissue %in% tissue_subset & data$trait %in% trait_subset],
#           trait = as.integer(as.factor(data$trait[data$tissue %in% tissue_subset & data$trait %in% trait_subset])),
#           tissue = as.integer(as.factor(data$tissue[data$tissue %in% tissue_subset & data$trait %in% trait_subset])),
#           TWAS_hit = as.integer(data$TWAS_Hit[data$tissue %in% tissue_subset & data$trait %in% trait_subset] == "YES") + 1,
#           n = nrow(data[data$tissue %in% tissue_subset & data$trait %in% trait_subset,]),
#           n_trait = length(trait_subset),
#           n_tissue = length(tissue_subset),
#           n_cond = 2
# )
# 
# d <- list(count = data$count[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES"],
#           total = data$total[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES"],
#           trait = as.integer(as.factor(data$trait[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES"])),
#           tissue = as.integer(as.factor(data$tissue[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES"])),
#           TWAS_hit = as.integer(data$TWAS_Hit[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES"] == "YES") + 1,
#           n = nrow(data[data$tissue %in% tissue_subset & data$trait %in% trait_subset & data$TWAS_Hit == "YES",]),
#           n_trait = length(trait_subset),
#           n_tissue = length(tissue_subset),
#           n_cond = 1
# )

#are cell counts just the product of their row and column counts?
expected_prob <- t(t(n_genes_in_nodes / total_number_of_possible_hits)) %*% 
  t(sig_twas_by_trait_genes / total_number_of_possible_hits)
expected_prob <- expected_prob[,colnames(n_deg_sigtwas_intersect)]
expected_count <- expected_prob * total_number_of_possible_hits
par(mfrow = c(2,1), mar = c(4,4,1,1))
plot(as.matrix(n_deg_sigtwas_intersect), (expected_count))
sample_count <- t(sapply(rownames(expected_prob), function(ri) sapply(colnames(expected_prob), function(ci) 
  rbinom(size = total_number_of_possible_hits, n = 1, prob = expected_prob[ri, ci]))))
plot(as.matrix(n_deg_sigtwas_intersect), (sample_count))
sample_prob <- sample_count / total_number_of_possible_hits


#observed
par(mfrow = c(2,2), mar = c(4.5,4.25,1,1))
remove_inf <- function(x) x[x != Inf & x != -Inf]
partway <- function(x, p = 0.35) x[1]*p + x[2]*(1-p)
plot(y = range(remove_inf(c(logit(as.matrix(n_deg_sigtwas_intersect) / total_number_of_possible_hits)))), 
     x = range(remove_inf(c(logit(expected_prob)))),
     col = "white", 
     pch = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))),
     ylab = "true sample logodds", xlab = "expected logodds")
text(y = c(logit(as.matrix(n_deg_sigtwas_intersect) / total_number_of_possible_hits)), x = c(logit(expected_prob)),
     col = rep(tissue_cols[rownames(expected_prob)], ncol(expected_prob)), 
     labels = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))))
text(partway(par("usr")[1:2]), par("usr")[4], labels = "numbers represent trait indices", pos = 1)
legend(x = "bottomleft", legend = rownames(expected_prob), col = tissue_cols[rownames(expected_prob)], 
       pch = "X", ncol = 1, pt.cex = 1, cex = 0.75)
legend(x = "topleft", legend = "1-to-1 line", lty = 2, col = adjustcolor(1,0.5), lwd = 2)
abline(0,1, lty = 2, col = adjustcolor(1,0.5), lwd = 2)

plot(y = range(c((as.matrix(n_deg_sigtwas_intersect) / total_number_of_possible_hits))), x = range(c((expected_prob))),
     col = "white", 
     pch = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))),
     ylab = "true sample frequencies", xlab = "expected frequencies")
text(y = c((as.matrix(n_deg_sigtwas_intersect) / total_number_of_possible_hits)), x = c((expected_prob)),
     col = rep(tissue_cols[rownames(expected_prob)], ncol(expected_prob)), 
     labels = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))))
text(mean(par("usr")[1:2]), par("usr")[3], labels = "numbers represent trait indices", pos = 3)
legend(x = "bottomright", legend = rownames(expected_prob), col = tissue_cols[rownames(expected_prob)], 
       pch = "X", ncol = 1, pt.cex = 1, cex = 0.75)
legend(x = "topleft", legend = "1-to-1 line", lty = 2, col = adjustcolor(1,0.5), lwd = 2)
abline(0,1, lty = 2, col = adjustcolor(1,0.5), lwd = 2)

#simulated
plot(y = range(remove_inf(c(logit(as.matrix(sample_count) / total_number_of_possible_hits)))), 
     x = range(remove_inf(c(logit(expected_prob)))),
     col = "white", 
     pch = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))),
     ylab = "simulated sample logodds", xlab = "expected logodds")
text(y = c(logit(as.matrix(sample_count) / total_number_of_possible_hits)), x = c(logit(expected_prob)),
     col = rep(tissue_cols[rownames(expected_prob)], ncol(expected_prob)), 
     labels = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))))
text(partway(par("usr")[1:2]), par("usr")[4], labels = "numbers represent trait indices", pos = 1)
legend(x = "bottomleft", legend = rownames(expected_prob), col = tissue_cols[rownames(expected_prob)], 
       pch = "X", ncol = 1, pt.cex = 1, cex = 0.75)
legend(x = "topleft", legend = "1-to-1 line", lty = 2, col = adjustcolor(1,0.5), lwd = 2)
abline(0,1, lty = 2, col = adjustcolor(1,0.5), lwd = 2)

plot(y = range(c((as.matrix(sample_count) / total_number_of_possible_hits))), x = range(c((expected_prob))),
     col = "white", 
     pch = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))),
     ylab = "simulated sample frequencies", xlab = "expected frequencies")
text(y = c((as.matrix(sample_count) / total_number_of_possible_hits)), x = c((expected_prob)),
     col = rep(tissue_cols[rownames(expected_prob)], ncol(expected_prob)), 
     labels = as.character(unlist(lapply(1:ncol(expected_prob), function(i) rep(i, nrow(expected_prob))))))
text(mean(par("usr")[1:2]), par("usr")[3], labels = "numbers represent trait indices", pos = 3)
legend(x = "bottomright", legend = rownames(expected_prob), col = tissue_cols[rownames(expected_prob)], 
       pch = "X", ncol = 1, pt.cex = 1, cex = 0.75)
legend(x = "topleft", legend = "1-to-1 line", lty = 2, col = adjustcolor(1,0.5), lwd = 2)
abline(0,1, lty = 2, col = adjustcolor(1,0.5), lwd = 2)


#alternate d

d <- list(count_1 = data1$count,
          total_1 = data1$total,
          count_2 = data2$count,
          total_2 = data2$total,
          trait = as.integer(as.factor(data1$trait)),
          tissue = as.integer(as.factor(data1$tissue)),
          n = nrow(data1),
          n_trait = length(unique(data1$trait)),
          n_tissue = length(unique(data1$tissue))
)

trait_subset <- names(sig_twas_by_trait_genes)[order(sig_twas_by_trait_genes, decreasing = T)][1:2]
tissue_subset <- names(n_genes_in_nodes)[order(n_genes_in_nodes, decreasing = T)][1:5]
d <- list(count_1 = data1$count[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
          total_1 = data1$total[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
          count_2 = data2$count[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
          total_2 = data2$total[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
          trait = as.integer(as.factor(data1$trait[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
          tissue = as.integer(as.factor(data1$tissue[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
          n = nrow(data1[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset,]),
          n_trait = length(unique(data1$trait[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
          n_tissue = length(unique(data1$tissue[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset]))
)


# trait_subset <- names(sig_twas_by_trait_genes)[order(sig_twas_by_trait_genes, decreasing = T)][1:20]
# tissue_subset <- names(n_genes_in_nodes)[order(n_genes_in_nodes, decreasing = T)][1:5]
# d <- list(count_1 = data1$count[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
#           total_1 = data1$total[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset],
#           count_2 = data2$count[data2$tissue %in% tissue_subset & data1$trait %in% trait_subset],
#           total_2 = data2$total[data2$tissue %in% tissue_subset & data1$trait %in% trait_subset],
#           trait = as.integer(as.factor(data1$trait[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
#           tissue = as.integer(as.factor(data1$tissue[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
#           n = nrow(data1[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset,]),
#           n_trait = length(unique(data1$trait[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset])),
#           n_tissue = length(unique(data1$tissue[data1$tissue %in% tissue_subset & data1$trait %in% trait_subset]))
# )


#load libraries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

## STAN model
# stan_program <- '
# data {
#     int<lower=1> n;
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=0> count_1[n];
#     int<lower=0> total_1[n];
#     int<lower=0> count_2[n];
#     int<lower=0> total_2[n];
#     int<lower=1> trait[n];
#     int<lower=1> tissue[n];
# }
# transformed data {
# 
# }
# parameters {
#     vector<lower=0,upper=1>[2] mean_beta;
#     vector<lower=0>[2] conc_beta;
#     vector<lower=0,upper=1>[n] prop_1;
#     vector<lower=0,upper=1>[n] prop_2;
# }
# transformed parameters {
#     vector<lower=0>[2] alpha = mean_beta .* conc_beta;
#     vector<lower=0>[2] beta = (1 - mean_beta) .* (conc_beta);
# }
# model {
#     //priors
#     mean_beta ~ beta(1,1);
#     conc_beta ~ exponential(0.01);
#     
#     //likelihood for obs
#     count_1 ~ binomial(total_1, prop_1);
#     count_2 ~ binomial(total_2, prop_2);
#     
#     //hierarchical model for props
#     prop_1 ~ beta(alpha[1], beta[1]);
#     prop_2 ~ beta(alpha[2], beta[2]);
# }
# generated quantities {
#   vector[n] difference_in_props = prop_1 - prop_2;
# }
# '

# # STAN model
# stan_program <- '
# data {
#     int<lower=1> n;
#     int<lower=1> n_cond;
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=0> count[n];
#     int<lower=0> total[n];
#     int<lower=1> trait[n];
#     int<lower=1> tissue[n];
#     int<lower=1> TWAS_hit[n];
# }
# transformed data {
#     vector<lower = 0>[3] concentration_offset = [2,2,2]\';
# }
# parameters {
#     //4th floor
#     real<lower=0,upper=1> mean_beta;
#     real<lower=0> conc_beta;
#     
#     //3rd floor
#     vector<lower=0,upper=1>[n_trait] mean_beta_trait;
#     real<lower=0> conc_beta_trait;
#     
#     //2nd floor
#     matrix<lower=0,upper=1>[n_tissue, n_trait] mean_beta_trait_tissue;
#     
#     real<lower=0> conc_beta_trait_tissue;
#     //real conc_beta_trait_tissue_logmean;
#     //vector[n_tissue] conc_beta_trait_tissue_tissue_eff;
#     //real<lower=0> conc_beta_trait_tissue_tissue_eff_sd;
#     //vector[n_trait] conc_beta_trait_tissue_trait_eff;
#     //real<lower=0> conc_beta_trait_tissue_trait_eff_sd;
#     
#     //1st floor
#     vector<lower=0,upper=1>[n] prop;
# }
# transformed parameters {
#     //4th floor
#     real<lower=0> alpha = mean_beta * (conc_beta + concentration_offset[1]);
#     real<lower=0> beta = (1 - mean_beta) * (conc_beta + concentration_offset[1]);
#     
#     //3rd floor
#     vector[n_trait] alpha_trait = mean_beta_trait * (conc_beta_trait + concentration_offset[2]);
#     vector[n_trait] beta_trait = (1 - mean_beta_trait) * (conc_beta_trait + concentration_offset[2]);
#     
#     //2nd floor
#     matrix<lower=0>[n_tissue, n_trait] alpha_trait_tissue = mean_beta_trait_tissue * (conc_beta_trait_tissue + concentration_offset[2]);
#     matrix<lower=0>[n_tissue, n_trait] beta_trait_tissue = (1 - mean_beta_trait_tissue) * (conc_beta_trait_tissue + concentration_offset[2]);
#     //real<lower=0> conc_beta_trait_tissue;
#     //matrix<lower=0>[n_tissue, n_trait]  conc_beta_trait_tissue = exp(conc_beta_trait_tissue_logmean + 
#     //                                                                 rep_matrix(conc_beta_trait_tissue_tissue_eff, n_trait) + 
#     //                                                                 rep_matrix(conc_beta_trait_tissue_trait_eff\', n_tissue));
#     //matrix<lower=0>[n_tissue, n_trait] alpha_trait_tissue = mean_beta_trait_tissue .* (conc_beta_trait_tissue + concentration_offset[3]);
#     //matrix<lower=0>[n_tissue, n_trait] beta_trait_tissue = (1 - mean_beta_trait_tissue) .* (conc_beta_trait_tissue + concentration_offset[3]);
# }
# model {
#     //priors
#     mean_beta ~ beta(1,1);
#     conc_beta ~ exponential(0.01);
#     
#     mean_beta_trait ~ beta(alpha, beta);
#     conc_beta_trait ~ exponential(0.01);
# 
#     for(tr_i in 1:n_trait){
#       mean_beta_trait_tissue[,tr_i] ~ beta(alpha_trait[tr_i], beta_trait[tr_i]);
#     }
#     conc_beta_trait_tissue ~ exponential(0.01);
#     //conc_beta_trait_tissue_logmean ~ normal(0, 5);
#     //conc_beta_trait_tissue_tissue_eff ~ normal(0, conc_beta_trait_tissue_tissue_eff_sd);
#     //conc_beta_trait_tissue_tissue_eff_sd ~ exponential(0.5);
#     //conc_beta_trait_tissue_trait_eff ~ normal(0, conc_beta_trait_tissue_trait_eff_sd);
#     //conc_beta_trait_tissue_trait_eff_sd ~ exponential(0.5);
#     
#     //for(i in 1:n){
#     //    prop[i] ~ beta(alpha_trait_tissue[tissue[i], trait[i]], beta_trait_tissue[tissue[i], trait[i]]);
#     //}
#     prop[1:(n%/%n_cond)] ~ beta(to_vector(alpha_trait_tissue), to_vector(beta_trait_tissue));
#     prop[(n%/%n_cond+1):n] ~ beta(to_vector(alpha_trait_tissue), to_vector(beta_trait_tissue));
#     //prop ~ beta(alpha_trait_tissue[tissue, trait], beta_trait_tissue[tissue, trait]);
#     
#     //likelihood for obs
#     count ~ binomial(total, prop);
#     
# }
# generated quantities {
#   vector[n%/%n_cond] difference_in_props = prop[1:(n%/%n_cond)] - prop[(n%/%n_cond+1):n];
# }
# '
# 
# 
# # STAN model with logit-normal 
# stan_program <- '
# data {
#     int<lower=1> n;
#     int<lower=1> n_cond;
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=0> count[n];
#     int<lower=0> total[n];
#     int<lower=1> trait[n];
#     int<lower=1> tissue[n];
#     int<lower=1> TWAS_hit[n];
# }
# parameters {
#     //4th floor
#     real mean_all_traits;
#     real<lower=0> sd_all_traits;
#     real mean_log_sd_each_trait;
#     real<lower=0> sd_log_sd_each_trait;
#     
#     //3rd floor
#     vector[n_trait] mean_each_trait;
#     vector[n_trait] log_sd_each_trait;
#     
#     //1st floor
#     vector<lower=0,upper=1>[n] prop;
# }
# transformed parameters {
#     vector[n] logit_prop = logit(prop);
#     vector<lower=0>[n_trait] sd_each_trait = exp(log_sd_each_trait);
# }
# model {
#     //priors
#     mean_all_traits ~ normal(0,5);
#     sd_all_traits ~ normal(0,2);
#     
#     //hiearchical model
#     mean_each_trait ~ normal(mean_all_traits, sd_all_traits);
#     sd_each_trait ~ normal(0,2);
#     
#     logit_prop ~ normal(mean_each_trait[trait], sd_each_trait[trait]);
#     
#     //likelihood for obs
#     count ~ binomial(total, prop);
#     
# }
# generated quantities {
#   vector[n%/%n_cond] difference_in_props = prop[1:(n%/%n_cond)] - prop[(n%/%n_cond+1):n];
# }
# '

# STAN model with logit-normal, alternate + non-centered parameterization
stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=0> count_1[n];
    int<lower=0> total_1[n];
    int<lower=0> count_2[n];
    int<lower=0> total_2[n];
    int<lower=1> trait[n];
    int<lower=1> tissue[n];
}
parameters {
    //3rd floor
    real mean_all_traits;
    real<lower=0> sd_all_traits;
    
    //2nd floor
    vector[n_trait] raw_mean_each_trait;
    real<lower=0> sd_each_trait;
    vector[n_trait] raw_sd_each_trait_log_multiplier;
    real<lower=0> sd_each_trait_log_multipler_sd; 
    
    //1st floor
    vector[n] raw_logit_prop_2;
    vector[n] raw_logit_prop_21;
    real<lower=0> logit_prop_21_sd;
}
transformed parameters {
    //noncentered parameters
    vector[n_trait] sd_each_trait_log_multiplier = raw_sd_each_trait_log_multiplier * sd_each_trait_log_multipler_sd;
    vector<lower=0>[n_trait] sd_each_trait_multiplier = exp(sd_each_trait_log_multiplier);
    vector[n_trait] mean_each_trait = raw_mean_each_trait * sd_all_traits + mean_all_traits;
    vector[n] logit_prop_2 = raw_logit_prop_2 .* sd_each_trait_multiplier[trait] * sd_each_trait + mean_each_trait[trait];
    vector[n] logit_prop_21 = raw_logit_prop_21 * logit_prop_21_sd;
    
    //other transformations
    vector[n] logit_prop_1 = logit_prop_2 + logit_prop_21;
    vector<lower=0,upper=1>[n] prop_2 = inv_logit(logit_prop_2);
    vector<lower=0,upper=1>[n] prop_1 = inv_logit(logit_prop_1);

}
model {
    //priors
    mean_all_traits ~ normal(0,2);
    sd_all_traits ~ normal(0,1);
    
    raw_mean_each_trait ~ std_normal();
    raw_sd_each_trait_log_multiplier ~ std_normal();
    sd_each_trait_log_multipler_sd ~ exponential(1);
    sd_each_trait ~ exponential(1);
    
    logit_prop_21_sd ~ normal(0,1);
    raw_logit_prop_21 ~ std_normal();
    raw_logit_prop_2 ~ std_normal();
    
    //likelihood for obs
    count_1 ~ binomial(total_1, prop_1);
    count_2 ~ binomial(total_2, prop_2);
}
'

# STAN model with logit-normal, non-centered parameterization
stan_program <- '
data {
    int<lower=1> n; //total # cells
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=0> count_1[n];
    int<lower=0> total_1[n];
    int<lower=0> count_2[n];
    int<lower=0> total_2[n];
    int<lower=1> trait[n]; //trait index for counts
    int<lower=1> tissue[n]; //column index for counts
}
parameters {
    //3rd floor
    real mean_all_traits;
    real<lower=0> sd_all_traits;
    
    //2nd floor
    vector[n_trait] raw_mean_each_trait;
    real<lower=0> sd_each_trait; 
    vector[n_trait] raw_sd_each_trait_log_multiplier;
    real<lower=0> sd_each_trait_log_multipler_sd; 
    
    //1st floor
    vector[n] raw_means_tissue_x_trait;
    real<lower=0> sd_tissue_x_trait;
    
    //basement
    vector[n] raw_logit_prop_1;
    vector[n] raw_logit_prop_2;
}
transformed parameters {
    //uncentering 2nd floor
    vector[n_trait] sd_each_trait_log_multiplier = raw_sd_each_trait_log_multiplier * sd_each_trait_log_multipler_sd;
    vector<lower=0>[n_trait] sd_each_trait_multiplier = exp(sd_each_trait_log_multiplier);
    vector[n_trait] mean_each_trait = raw_mean_each_trait * sd_all_traits + mean_all_traits;
    
    //uncentering 1st floor
    vector[n] means_tissue_x_trait = raw_means_tissue_x_trait .* sd_each_trait_multiplier[trait] * sd_each_trait + mean_each_trait[trait];
    
    //uncentering & inv_logit-ing basement
    vector[n] logit_prop_1 = raw_logit_prop_1 * sd_tissue_x_trait + means_tissue_x_trait;
    vector[n] logit_prop_2 = raw_logit_prop_2 * sd_tissue_x_trait + means_tissue_x_trait;
    vector<lower=0,upper=1>[n] prop_1 = inv_logit(logit_prop_1);
    vector<lower=0,upper=1>[n] prop_2 = inv_logit(logit_prop_2);
}
model {
    //3rd floor
    mean_all_traits ~ normal(-1,2);
    sd_all_traits ~ normal(0,0.5);
    
    //2nd floor
    raw_mean_each_trait ~ std_normal();
    raw_sd_each_trait_log_multiplier ~ std_normal();
    sd_each_trait_log_multipler_sd ~ normal(0,0.1);
    sd_each_trait ~ exponential(1);
    
    //1st floor
    raw_means_tissue_x_trait ~ std_normal();
    sd_tissue_x_trait ~ normal(0,0.5);
    
    //basement
    raw_logit_prop_1 ~ std_normal();
    raw_logit_prop_2 ~ std_normal();
    
    //likelihood for obs
    count_1 ~ binomial(total_1, prop_1);
    count_2 ~ binomial(total_2, prop_2);
}
generated quantities {
  vector[n] difference_in_props = prop_1 - prop_2;
}
'

# STAN model with logit-normal, non-centered parameterization, Bob's model
d <- list(count = rbind(data1$count, data2$count),
          total = rbind(data1$total,data2$total),
          trait = as.integer(as.factor(data1$trait)),
          tissue = as.integer(as.factor(data1$tissue)),
          n = nrow(data1),
          n_trait = length(unique(data1$trait)),
          n_tissue = length(unique(data1$tissue)),
          n_cond = 2
)

stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=1> n_cond;
    int<lower=0> count[n_cond,n];
    int<lower=0> total[n_cond,n];
    int<lower=1> trait[n]; //trait index for counts
    int<lower=1> tissue[n]; //column index for counts
}
parameters {
    real grand_mean;

    vector[n_trait] raw_trait_mean;
    real<lower=0> trait_sd;
    
    vector[n_tissue] raw_tissue_mean;
    real<lower=0> tissue_sd;
   
    vector[n] raw_cond_mean;
    real<lower=0> cond_sd;
    
    array[n_cond] vector[n] raw_logodds;
    real<lower=0> logodds_sd;
}
transformed parameters {
    //recenter parameters
    vector[n_trait] trait_mean = raw_trait_mean * trait_sd;
    vector[n_tissue] tissue_mean = raw_tissue_mean * tissue_sd;
    vector[n] cond_mean = raw_cond_mean * cond_sd + grand_mean + trait_mean[trait] + tissue_mean[tissue];
    array[n_cond] vector[n] logodds;
    for (i in 1:n_cond){
      logodds[i] = raw_logodds[i] * logodds_sd + cond_mean;
    }
}
model {
    grand_mean ~ normal(0, 2);

    trait_sd ~ normal(0, 0.5);
    trait_mean ~ std_normal();
    
    tissue_sd ~ normal(0, 0.5);
    tissue_mean ~ std_normal();
    
    cond_sd ~ normal(0, 0.5);
    cond_mean ~ std_normal();
    
    logodds_sd ~ normal(0, 0.5);
    for (i in 1:n_cond){
      raw_logodds[i] ~ std_normal();
      count[i] ~ binomial_logit(total[i], logodds[i]);
    }
    
}
generated quantities {
  vector[n] difference_in_props = inv_logit(logodds[1]) - inv_logit(logodds[2]);
}
'

# STAN model with logit-normal, non-centered parameterization, nest traits in tissues
stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=0> count_1[n];
    int<lower=0> total_1[n];
    int<lower=0> count_2[n];
    int<lower=0> total_2[n];
    int<lower=1,upper=n_trait> trait[n];
    int<lower=1,upper=n_tissue> tissue[n];
}
parameters {
    //3rd floor
    real mean_all_tissues;
    real<lower=0> sd_all_tissues;
    
    //2nd floor
    vector[n_tissue] raw_mean_each_tissue;
    real<lower=0> sd_each_tissue; 
    //multilevel trait effects
    vector[n_tissue] raw_sd_each_tissue_log_multiplier;
    real<lower=0> sd_each_tissue_log_multipler_sd; 
    
    //1st floor
    vector[n] raw_means_trait_x_tissue;
    real<lower=0> sd_trait_x_tissue;
    //multilevel trait effects
    vector[n_tissue] raw_sd_TxT_log_tissue_multiplier;
    real<lower=0> sd_TxT_log_tissue_multipler_sd; 
    vector[n_trait] raw_sd_TxT_log_trait_multiplier;
    real<lower=0> sd_TxT_log_trait_multipler_sd; 
    
    //basement
    vector[n] raw_logit_prop_1;
    vector[n] raw_logit_prop_2;
}
transformed parameters {
    //uncentering 2nd floor
    vector[n_tissue] sd_each_tissue_log_multiplier = raw_sd_each_tissue_log_multiplier * sd_each_tissue_log_multipler_sd;
    vector<lower=0>[n_tissue] sd_each_tissue_multiplier = exp(sd_each_tissue_log_multiplier);
    vector[n_tissue] mean_each_tissue = raw_mean_each_tissue * sd_all_tissues + mean_all_tissues;
    
    //uncentering 1st floor
    vector[n] means_trait_x_tissue = raw_means_trait_x_tissue .* sd_each_tissue_multiplier[tissue] * sd_each_tissue + mean_each_tissue[tissue];
    vector[n_tissue] sd_TxT_log_tissue_multiplier = raw_sd_TxT_log_tissue_multiplier * sd_TxT_log_tissue_multipler_sd;
    vector<lower=0>[n_tissue] sd_TxT_tissue_multiplier = exp(sd_TxT_log_tissue_multiplier);
    vector[n_trait] sd_TxT_log_trait_multiplier = raw_sd_TxT_log_trait_multiplier * sd_TxT_log_trait_multipler_sd;
    vector<lower=0>[n_trait] sd_TxT_trait_multiplier = exp(sd_TxT_log_trait_multiplier);
    
    //uncentering & inv_logit-ing basement
    vector[n] logit_prop_1 = raw_logit_prop_1 .* sd_TxT_tissue_multiplier[tissue] .* sd_TxT_trait_multiplier[trait] * sd_trait_x_tissue + means_trait_x_tissue;
    vector[n] logit_prop_2 = raw_logit_prop_2 .* sd_TxT_tissue_multiplier[tissue] .* sd_TxT_trait_multiplier[trait] * sd_trait_x_tissue + means_trait_x_tissue;
    vector<lower=0,upper=1>[n] prop_1 = inv_logit(logit_prop_1);
    vector<lower=0,upper=1>[n] prop_2 = inv_logit(logit_prop_2);
}
model {
    //3rd floor
    mean_all_tissues ~ normal(-1,2);
    sd_all_tissues ~ normal(0,0.5);
    
    //2nd floor
    raw_mean_each_tissue ~ std_normal();
    raw_sd_each_tissue_log_multiplier ~ std_normal();
    sd_each_tissue_log_multipler_sd ~ normal(0,0.5);
    sd_each_tissue ~ exponential(1);
    
    //1st floor
    raw_means_trait_x_tissue ~ std_normal();
    sd_trait_x_tissue ~ normal(0,0.5);
    raw_sd_TxT_log_tissue_multiplier ~ std_normal();
    raw_sd_TxT_log_trait_multiplier ~ std_normal();
    sd_TxT_log_tissue_multipler_sd ~ normal(0,0.5);
    sd_TxT_log_trait_multipler_sd ~ normal(0,0.5);
    
    //basement
    raw_logit_prop_1 ~ std_normal();
    raw_logit_prop_2 ~ std_normal();
    
    //likelihood for obs
    count_1 ~ binomial(total_1, prop_1);
    count_2 ~ binomial(total_2, prop_2);
}
generated quantities {
  vector[n] difference_in_props = prop_1 - prop_2;
}
'

# stan_program <- '
# data {
#     int<lower=1> n;
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=0> count_1[n];
#     int<lower=0> total_1[n];
#     int<lower=0> count_2[n];
#     int<lower=0> total_2[n];
#     int<lower=1> trait[n];
#     int<lower=1> tissue[n];
# }
# parameters {
#     //3rd floor
#     real mean_all_traits;
#     real<lower=0> sd_all_traits;
#     real mean_log_sd_each_trait;
#     real<lower=0> sd_log_sd_each_trait;
#     
#     //2nd floor
#     vector[n_trait] mean_each_trait;
#     vector[n_trait] log_sd_each_trait;
#     
#     //1st floor
#     vector[n] logit_prop_2;
#     vector[n] logit_prop_21;
#     real<lower=0> logit_prop_21_sd;
# }
# transformed parameters {
#     vector[n] logit_prop_1 = logit_prop_2 + logit_prop_21;
#     vector<lower=0,upper=1>[n] prop_2 = inv_logit(logit_prop_2);
#     vector<lower=0,upper=1>[n] prop_1 = inv_logit(logit_prop_1);
#     vector<lower=0>[n_trait] sd_each_trait = exp(log_sd_each_trait);
# }
# model {
#     //priors
#     mean_all_traits ~ normal(0,2);
#     sd_all_traits ~ normal(0,1);
#     mean_log_sd_each_trait ~ normal(0,0.5);
#     sd_log_sd_each_trait ~ normal(0,0.25);
#     
#     mean_each_trait ~ normal(mean_all_traits, sd_all_traits);
#     log_sd_each_trait ~ normal(mean_log_sd_each_trait,sd_log_sd_each_trait);
#     
#     logit_prop_21_sd ~ normal(0,1);
#     logit_prop_21 ~ normal(0,logit_prop_21_sd);
#     logit_prop_2 ~ normal(mean_each_trait[trait], sd_each_trait[trait]);
#     
#     //likelihood for obs
#     count_1 ~ binomial(total_1, prop_1);
#     count_2 ~ binomial(total_2, prop_2);
# }
# '

#reparameterizing according to joint probability
d <- list(intersect_count = data1$count,
          total = total_number_of_possible_hits,
          trait_count = sig_twas_by_trait_genes[levels(as.factor(data1$trait))],
          tissue_count = n_genes_in_nodes[levels(as.factor(data1$tissue))],
          trait = as.integer(as.factor(data1$trait)),
          tissue = as.integer(as.factor(data1$tissue)),
          n_trait = length(unique(data1$trait)),
          n_tissue = length(unique(data1$tissue))
)

stan_program <- '
data {
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=0> total;
    int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
    int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
    int<lower=0> trait_count[n_trait];
    int<lower=0> tissue_count[n_tissue];
    int<lower=0> intersect_count[n_trait * n_tissue];
}
transformed data {
    int<lower=1> n = n_trait * n_tissue;
}
parameters {
    //main parameters
    real trait_mean;
    real<lower=0> trait_sd;
    vector[n_trait] raw_trait_logodds;

    real tissue_mean;
    real<lower=0> tissue_sd;
    vector[n_tissue] raw_tissue_logodds;
    
    vector[n] raw_intersect_logodds;
    real<lower=0> intersect_sd;
    
    //biases in deviations terms
    vector[n_tissue] tissue_bias_mean;
    vector[n_tissue] tissue_bias_log_sd;
    vector[n_trait] trait_bias_mean;
    vector[n_trait] trait_bias_log_sd;
    
    //multilevel deviation term params
    real<lower=0> tissue_bias_mean_sd;
    real<lower=0> trait_bias_mean_sd;
}
transformed parameters {
    //recenter params
    vector[n_trait] trait_logodds = raw_trait_logodds * trait_sd + trait_mean;
    vector[n_tissue] tissue_logodds = raw_tissue_logodds * tissue_sd + tissue_mean;
    vector[n] intersect_mean = logit(inv_logit(trait_logodds[trait]) .* inv_logit(tissue_logodds[tissue]));
    vector[n] intersect_logodds = raw_intersect_logodds * intersect_sd .* exp(tissue_bias_log_sd[tissue]) .* exp(trait_bias_log_sd[trait]) +
                                  intersect_mean + tissue_bias_mean[tissue] * tissue_bias_mean_sd + 
                                  trait_bias_mean[trait] * trait_bias_mean_sd;
}
model {
    //priors / hyperpriors
    raw_trait_logodds ~ std_normal();
    trait_mean ~ normal(0,2);
    trait_sd ~ normal(0,2);
    
    raw_tissue_logodds ~ std_normal();
    tissue_mean ~ normal(0,2);
    tissue_sd ~ normal(0,2);
    
    raw_intersect_logodds ~ std_normal();
    intersect_sd ~ normal(0,2);
    
    tissue_bias_mean ~ std_normal();
    tissue_bias_log_sd ~ std_normal();
    trait_bias_mean ~ std_normal();
    trait_bias_log_sd ~ std_normal();
    tissue_bias_mean_sd ~ std_normal();
    trait_bias_mean_sd ~ std_normal();
    
    //likelihood
    trait_count ~ binomial_logit(total, trait_logodds);
    tissue_count ~ binomial_logit(total, tissue_logodds);
    intersect_count ~ binomial_logit(total, intersect_logodds);
}
generated quantities {
    vector[n] difference_in_props = inv_logit(intersect_logodds) - inv_logit(intersect_mean);
}
'

#subbing in other joint probability model
d <- list(cell_count = c(cell_count),
          total = total,
          row_count = row_count,
          col_count = col_count,
          row_index = rep(1:row_n, col_n),
          col_index = unlist(lapply(1:col_n, function(i) rep(i, row_n))),
          row_n = row_n,
          col_n = col_n)

stan_program <- '
data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=0> total;
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=0> row_count[row_n];
    int<lower=0> col_count[col_n];
    int<lower=0> cell_count[row_n * col_n];
}
transformed data {
    int<lower=1> n = row_n * col_n;
    int<lower=0,upper=1> smaller_margin[n]; //0 if row, 1 if col
    int<lower=0> marginal_total[n];
    for(i in 1:n){
      //smaller_margin[i] = 1;
      smaller_margin[i] = row_count[row_index[i]] > col_count[col_index[i]];
      marginal_total[i] = smaller_margin[i] * col_count[col_index[i]] + (1-smaller_margin[i]) * row_count[row_index[i]];
    }
    vector[n] smaller_margin_vec = to_vector(smaller_margin);
}
parameters {
    //main parameters
    real col_mean;
    real<lower=0> col_sd;
    vector[col_n] raw_col_logodds;

    real row_mean;
    real<lower=0> row_sd;
    vector[row_n] raw_row_logodds;
    
    vector[n] raw_cell_logodds;
    real<lower=0> cell_sd;
    
    //biases in deviations terms
    vector[row_n] raw_row_bias;
    vector[col_n] raw_col_bias;

    //multilevel deviation term params
    real<lower=0> row_bias_sd;
    real<lower=0> col_bias_sd;
}
transformed parameters {
    //recenter params
    vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
    vector[row_n] row_logodds = raw_row_logodds * row_sd + row_mean;
    vector[n] cell_mean = smaller_margin_vec .* row_logodds[row_index] + (1-smaller_margin_vec) .* col_logodds[col_index];
    vector[n] cell_logodds = raw_cell_logodds * cell_sd + cell_mean + 
                             raw_row_bias[row_index] * row_bias_sd + 
                             raw_col_bias[col_index] * col_bias_sd;
}
model {
    //priors and hyperpriors
    
    //marginal params
    raw_col_logodds ~ std_normal();
    col_mean ~ normal(0,2);
    col_sd ~ std_normal();
    
    raw_row_logodds ~ std_normal();
    row_mean ~ normal(0,2);
    row_sd ~ std_normal();
    
    //bias params
    raw_row_bias ~ std_normal();
    raw_col_bias ~ std_normal();
    row_bias_sd ~ std_normal();
    col_bias_sd ~ std_normal();
    
    //cell params
    raw_cell_logodds ~ std_normal();
    cell_sd ~ std_normal();
    
    //likelihood
    col_count ~ binomial_logit(total, col_logodds);
    row_count ~ binomial_logit(total, row_logodds);
    cell_count ~ binomial_logit(marginal_total, cell_logodds);
}
generated quantities {
    vector[n] cell_bias = cell_logodds - (cell_mean + 
                             raw_row_bias[row_index] * row_bias_sd + 
                             raw_col_bias[col_index] * col_bias_sd);
    vector[n] cell_total_prob_bias = inv_logit(cell_logodds) - inv_logit(cell_mean);
    vector[row_n] row_bias = raw_row_bias * row_bias_sd;
    vector[col_n] col_bias = raw_col_bias * col_bias_sd;
}
'


if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)


#fit model
write_stan_file(stan_program, dir = "~/Desktop/", basename = "deviation_from_expected_logodds")
write_stan_json(d, "~/Desktop/deviation_from_expected_logodds.json")
out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, data = d, parallel_chains = 4, 
                  adapt_delta = 0.9, refresh = 50, init = 0.1, max_treedepth = 20, thin = 5)
# out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.9, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]
samps <- data.frame(as_draws_df(out$draws()))


# #look at pathology
# par(mfrow = c(5,3), mar = c(4,4,2,1))
# for(trait_i in 1:12){
#   plot(c((data1$count / data1$total)[data1$trait == unique(data1$trait)[trait_i]], 
#          (data2$count / data2$total)[data1$trait == unique(data1$trait)[trait_i]]),
#        col = c(rep("blue", length((data1$count / data1$total)[data1$trait == unique(data1$trait)[trait_i]])),
#                rep("red", length((data1$count / data1$total)[data1$trait == unique(data1$trait)[trait_i]]))),
#        pch = 19, cex = 1.5, ylab = "sample proportion", main = trait_categories$new_Phenotype[match(unique(data1$trait)[trait_i], trait_categories$Tag)])
#   legend(x = "bottomright", pt.cex = 1.5, cex = 1, col = c("blue", "red"), legend = c("focal proportion", "complement proportion"), pch = 19)
# }
# 
# par(mfrow = c(5,3), mar = c(4,4,2,1))
# for(tissue_i in 1:15){
#   plot(c((data1$count / data1$total)[data1$tissue == unique(data1$tissue)[tissue_i]], 
#          (data2$count / data2$total)[data1$tissue == unique(data1$tissue)[tissue_i]]),
#        col = c(rep("blue", length((data1$count / data1$total)[data1$tissue == unique(data1$tissue)[tissue_i]])),
#                rep("red", length((data1$count / data1$total)[data1$tissue == unique(data1$tissue)[tissue_i]]))),
#        pch = 19, cex = 1.5, ylab = "sample proportion", main = unique(data1$tissue)[tissue_i])
#   legend(x = "topright", pt.cex = 1.5, cex = 0.75, col = c("blue", "red"), legend = c("focal", "comp."), pch = 19)
# }
# 
# hist(cbind(data1$count / data1$total, data2$count / data2$total)[data1$trait == unique(data1$trait)[trait_i],])

#examine posterior output
# 
# par(mfrow = c(4,2))
# hist(invlogit(samps$mean_all_traits), 
#      main = "invlogit(mean_all_traits)", xlab = "value")
# hist((samps$sd_all_traits), 
#      main = "sd_all_traits", xlab = "value")
# hist(apply(invlogit(samps[,setdiff(grep("mean_each_trait", colnames(samps)), 
#                                    grep("raw_mean_each_trait", colnames(samps)))]), 2, mean), 
#      main = "posterior means of invlogit(mean_each_trait)", xlab = "value")
# hist((samps$sd_each_trait), 
#      main = "sd_each_trait", xlab = "value")
# hist((samps$sd_each_trait_log_multipler_sd), 
#      main = "sd_each_trait_log_multipler_sd", xlab = "value")
# hist((samps$sd_tissue_x_trait), 
#      main = "sd_tissue_x_trait", xlab = "value")
# hist(apply(samps[,grep("difference_in_props", colnames(samps))], 2, mean), 
#      main = "posterior means of difference_in_props", xlab = "value")


prop_greater_than_0 <- function(x) mean(x>0)
logit <- function(p) log(p/(1-p))
invlogit <- function(x) exp(x) / (1+exp(x))

# samps[,grep("logit_prop_21\\.", colnames(samps))]
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) > 0.95)
# sum(apply(samps[,grep("logit_prop_21\\.", colnames(samps))], 2, prop_greater_than_0) < 0.05)

dev.off()
hist(apply(samps[,grep("difference_in_props", colnames(samps))], 2, mean))
prop_greater_than_0 <- function(x) mean(x>0)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) < 0.05)

hist((apply(samps[,grep("difference_in_props", colnames(samps))], 2, mean)) / 
       (apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean)))

hist(apply(samps[,grep("trait_bias_mean", colnames(samps))], 2, mean))

apply(invlogit(samps[,grep("intersect_mean", colnames(samps))]), 2, mean)
hist(apply(invlogit(samps[,grep("trait_logodds", colnames(samps))]), 2, mean))

#generate a "significance" matrix
signif_threshold <- 0.05
signif_df <- data1[1:(nrow(data1)), c("tissue", "trait")]
signif_df$prob_diff_is_positive <- apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)
signif_df$signif <- 0
signif_df$signif[signif_df$prob_diff_is_positive > (1-0.05)] <- 1
signif_df$signif[signif_df$prob_diff_is_positive < 0.05] <- -1

signif_df[signif_df$signif != 0,]

signif_matrix <- reshape(signif_df[,-match("prob_diff_is_positive", colnames(signif_df))], idvar = "tissue", timevar = "trait", direction = "wide")
rownames(signif_matrix) <- signif_matrix$tissue
signif_matrix <- signif_matrix[,-match("tissue", colnames(signif_matrix))]
colnames(signif_matrix) <- gsub(colnames(signif_matrix), pattern = "signif.", replacement = "")
if(!all(colnames(signif_matrix) == colnames(n_deg_sigtwas_intersect)) & all(rownames(signif_matrix) == rownames(n_deg_sigtwas_intersect))){
  stop("something's wrong with the significance matrix")
}

#trait proportion matrix
# prop_twas_are_degs <- prop_twas_are_degs[,colnames(n_deg_sigtwas_intersect)[order(apply(samps[,grep("mean_beta_trait\\.", colnames(samps))], 2, mean), decreasing = T)]]
# par(mfrow = c(2,1))
# hist((data$count / data$total)[data$TWAS_Hit == "YES"], breaks = 0:100/100)
# 
# 
# hist(samps[,grep("mean_all", colnames(samps))])
# hist(samps[,grep("sd_all", colnames(samps))])
# hist(apply(invlogit(samps[,grep("mean_each_trait", colnames(samps))]), 2, mean))
# hist(apply(invlogit(samps[,grep("mean_tiss", colnames(samps))]), 2, mean))
# hist(samps[,grep("sd_tiss", colnames(samps))])
# hist(apply(invlogit(samps[,grep("sd_tiss", colnames(samps))]), 2, mean))




#### now plot the intersect of DEGs & TWAS hits ####

#counts or props?
incl_significance <- T
trait_category_legend_below = T
use_tissue_cols_for_cols <- T
opacity_power_scaler <- 0.25
opacity_white_threshold <- 0.8
use_counts <- T
prop_TWAS <- T
order_by_counts <- F
group_by_tissue_type <- T


if(use_counts){
  table_to_use <- n_deg_sigtwas_intersect  
} else {
    if(prop_TWAS){
      table_to_use <- round(prop_twas_are_degs * 1000)  
    } else{
      table_to_use <- round(prop_degs_are_twas * 1000)  
    }
}

if(order_by_counts){
  table_to_use <- table_to_use[,colnames(n_deg_sigtwas_intersect)]
  signif_matrix <- signif_matrix[,colnames(n_deg_sigtwas_intersect)]
} else {
  table_to_use <- table_to_use[,colnames(prop_twas_are_degs)]
  signif_matrix <- signif_matrix[,colnames(prop_twas_are_degs)]
}


if(group_by_tissue_type){
  tissue_cats <- list(circulation = c("BLOOD", "HEART", "SPLEEN"),
                      skeletal_muscle = c("SKM-GN", "SKM-VL"),
                      adipose = c("WAT−SC"),
                      other = rev(c("ADRNL", "KIDNEY", "LUNG", "LIVER")),
                      brain = c("CORTEX", "HYPOTH", "HIPPOC"),
                      GI = c("SMLINT", "COLON"))
  tissue_cats <- rev(tissue_cats)
  disp_amount <- 0.5
  tissue_disps <- unlist(lapply(1:length(tissue_cats), function(tci) rep(disp_amount * (tci), length(tissue_cats[[tci]]))))
  tissue_cats_bars_ylocs <- cbind(start = (c(0, cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount)) + 2 * disp_amount)[-(length(tissue_cats)+1)], 
                            end = cumsum(unlist(lapply(tissue_cats, function(tc) length(tc))) + disp_amount) + disp_amount)
  
}

category_colors <- RColorBrewer::brewer.pal(length(salient.categories), "Dark2")
names(category_colors) <- salient.categories

cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DEG-TWAS_Intersect", ifelse(use_counts, "_counts", "_permille"),".pdf"), 
          width = 2100 / 72, height = 500 / 72 + ifelse(group_by_tissue_type, disp_amount * 0.75, 0), family="Arial Unicode MS")
par(xpd = T, mar = c(6,0,6 + ifelse(group_by_tissue_type, disp_amount * 4.5, 0),5.5))
plot(1, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", xlim= c(-5,ncol(table_to_use)), ylim = c(-5,nrow(table_to_use)))

if(group_by_tissue_type){
  tissue_cats_bars_xlocs <- sapply(tissue_cats, function(tc) max(strwidth(tc, units = "user"))) + 0.2
  tissue_cats_bars_xlocs <- rep(max(tissue_cats_bars_xlocs), length(tissue_cats_bars_xlocs))
}


if(use_tissue_cols_for_cols){
  heatmap_cols <- sapply((1:max(table_to_use) / max(table_to_use))^opacity_power_scaler, function(opcty) 
    adjustcolor("black", opcty))
} else {
  heatmap_cols <- viridisLite::viridis(n = max(table_to_use, na.rm = T)*100+1)
  heatmap_cols <- heatmap_cols[round(log(1:max(table_to_use, na.rm = T)) / log(max(table_to_use, na.rm = T)) * max(table_to_use, na.rm = T) * 100 + 1)]
}
for(ri in 1:nrow(table_to_use)){
  text(x = 0.5, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 2, labels = rownames(table_to_use)[ri])
  for(ci in 1:ncol(table_to_use)){
    if(ri == 1){
      text(x = ci+0.5, y = -0.9, pos = 2, srt = 45,
                     labels = trait_categories$new_Phenotype[match(colnames(table_to_use)[ci], trait_categories$Tag)])
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 1/2 - 1,
           ytop =  ri + 1/2 - 1,
           col = category_colors[trait_categories$Category[match(colnames(table_to_use)[ci], trait_categories$Tag)]])
    }
    
    #vertical total # options
    if(ri == nrow(table_to_use)){
      text(x = ci-0.5, y = ri+1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 4, srt = 45,
           labels = sig_twas_by_trait_genes[colnames(table_to_use)[ci]])
    }
    
    #horiz total # of options
    if(ci == ncol(table_to_use)){
      text(x = ci+0.45, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), pos = 4,
           labels = n_genes_in_nodes[ri])
    }
    
    #the actual cells
    if(use_tissue_cols_for_cols){
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = adjustcolor(MotrpacBicQC::tissue_cols[rownames(table_to_use)[ri]], (table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler ))  
    } else {
      rect(xleft = ci + 1/2,
           xright = ci - 1/2,
           ybottom = ri - 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           ytop =  ri + 0.5 + ifelse(group_by_tissue_type, tissue_disps[ri], 0),
           col = heatmap_cols[table_to_use[ri, ci]])
    }
    
    
    # rect(xleft = ci + 1/2,
    #      xright = ci - 1/2,
    #      ybottom = ri - 0.475,
    #      ytop =  ri + 0.475,
    #      col = heatmap_cols[table_to_use[ri, ci]],
    #      border = MotrpacBicQC::tissue_cols[rownames(table_to_use)[ri]])
    
    #text inside of cells
    if(table_to_use[ri, ci] != 0){
      text_in_cell <- table_to_use[ri, ci]
      if(incl_significance){
        signif_dir <- c("₋", "","⁺")[match(signif_matrix[ri, ci], -1:1)]
        text_in_cell <- paste0(text_in_cell, signif_dir)
      }
      if(use_tissue_cols_for_cols){
        text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), 
             col = ifelse((table_to_use[ri, ci] / max(table_to_use, na.rm = T))^opacity_power_scaler > opacity_white_threshold, "white", "black"), cex = 0.85)  
      } else {
        text(text_in_cell, x = ci, y = ri + ifelse(group_by_tissue_type, tissue_disps[ri], 0), col = "white", cex = 0.85)
      }
    }
    
  }
}

#tissue category bars
if(group_by_tissue_type){
  for(bi in 1:length(tissue_cats_bars_xlocs)){
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi],
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    bracket_length <- 0.2
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,1],
             y1 = tissue_cats_bars_ylocs[bi,1],
             lwd = 3)
    segments(x0 = -tissue_cats_bars_xlocs[bi],
             x1 = -tissue_cats_bars_xlocs[bi]+bracket_length,
             y0 = tissue_cats_bars_ylocs[bi,2],
             y1 = tissue_cats_bars_ylocs[bi,2],
             lwd = 3)
    
    text(x = -tissue_cats_bars_xlocs[bi], y = mean(tissue_cats_bars_ylocs[bi,]), pos = 2, 
         labels = gsub("Gi", "GI", stringr::str_to_title(gsub("_", " ", names(tissue_cats)[bi]))), cex = 1.25)
  }  
}

#legend for heatmap
x_adj <- 2.25
y_adj <- 0
yb_adjust <- ifelse(incl_significance, 3, 0)
n_legend_rects_to_use <- 30
n_legend_labels_to_use <- 10
legend_yvals <- round(seq(0, max(table_to_use), length.out = n_legend_labels_to_use))
legend_ylocs <- seq(yb_adjust, 1 + nrow(table_to_use) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), length.out = n_legend_rects_to_use)
legend_ycols <- round(seq(1, max(table_to_use), length.out = n_legend_rects_to_use))
for(i in 1:(n_legend_rects_to_use-1)){
  yb = legend_ylocs[i]
  yt = legend_ylocs[i+1]
  print(paste(c(yb,yt)))
  rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
                xright = ncol(table_to_use) + x_adj + 2 - 1/2,
                ybottom = yb,
                ytop =  yt,
                col = heatmap_cols[legend_ycols[i]], border = NA)
}
for(i in 1:n_legend_labels_to_use){
  text(labels = legend_yvals[i], x = ncol(table_to_use) + x_adj + 2.4, pos = 4, cex = 0.75,
              y = -0.25 + yb_adjust + legend_yvals[i] / max(table_to_use) * (1+nrow(table_to_use) - yb_adjust + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0)))
}

#legend for significance
if(incl_significance){
  text(labels = paste0("Pr(Enr.) > ", (1 - signif_threshold), " : X⁺\nPr(Dep.) > ", (1 - signif_threshold), " : X₋"),
       x = ncol(table_to_use) + x_adj,
       y = yb_adjust-2, pos = 4, cex = 0.75)
}

rect(xleft = ncol(table_to_use) + x_adj + 2 + 1/2,
     xright = ncol(table_to_use) + x_adj + 2 - 1/2,
     ybottom = min(legend_ylocs) - diff(range(legend_ylocs)) / max(table_to_use) * 5,
     ytop =  max(legend_ylocs))

#legend for trait categories
if(trait_category_legend_below){
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = -1/2 + i + cumsum(c(0,strwidth(names(category_colors), units = "user") + 1))[i],
         xright = 1/2 + i + cumsum(c(0,strwidth(names(category_colors), units = "user") + 1))[i],
         ybottom = -11,
         ytop = -10,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = -10.5, x = 0.35 + i + cumsum(c(0,strwidth(names(category_colors), units = "user") + 1))[i])
    #draw circle?
    arc(t1 = 3*pi/2, t2 = pi/2, r1 = 10.5 / 2, r2 = 10.5 / 2, center = c(0,-10.5/2), lwd = 2, res = 100, adjx = ifelse(order_by_counts, 1, 1.25))
    points(0, 0, pch = -9658, cex = 2)
  }
} else {
  x_adj2 <- x_adj - 2.5
  y_adj2 <- y_adj - 2.5
  for(i in 1:length(category_colors)){
    rect(xleft = 0 + i,
         xright = 1 + i,
         ybottom = -10,
         ytop = -9,
         col = category_colors[i], border = 1)
    text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
  }
}

#labels
#horiz label for total
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0),
         y1 = nrow(table_to_use) + 1 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = -2, x1 = 0.5, y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 2.5 + ifelse(order_by_counts, 0, -1) + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = -2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 2, labels = "total # of\nTWAS hits")

#vertical label for total
segments(x0 = ncol(table_to_use) + 2, x1 = ncol(table_to_use) + 0.75, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
segments(x0 = ncol(table_to_use) + 2, x1 = ncol(table_to_use) + 1.75, 
         y0 = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), 
         y1 = nrow(table_to_use) + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), lwd = 2)
text(x = ncol(table_to_use) + 2, y = nrow(table_to_use) + 3 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, labels = "total # of DEGs")


#legend and title labels
text(labels = ifelse(use_counts, latex2exp::TeX("n_{intersect}"), latex2exp::TeX("‰_{ intersect}")), pos = 3, font = 2, cex = 1.25,
     x = ncol(table_to_use) + x_adj + 2, y = nrow(table_to_use) + y_adj + 0.75 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0))
if(use_counts){
  text(latex2exp::TeX(paste0("number of genes in 8w - F↓M↓+ 8w - F↑M↑ with IHW significant TWAS at $\\alpha$ = 0.05")), 
     x = 42.5, y = nrow(table_to_use) + 2.5 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, cex = 2.35, font = 2)
} else {
  text(latex2exp::TeX(paste0("proportion (‰) of IHW significant TWAS hits at $\\alpha$ = 0.05 in 8w - F↓M↓ or 8w - F↑M↑")), 
     x = 42.5, y = nrow(table_to_use) + 2.5 + ifelse(group_by_tissue_type, tissue_disps[length(tissue_disps)], 0), pos = 3, cex = 2.35, font = 2)
}
dev.off()

#### plot a few Beta distributions ####

xr <- 0:1000/1000
dens <- dbeta(xr, 5,15)
par(mfrow = c(1,1))
plot(xr, dens, type = "l", lwd = 2, frame = F, yaxt = "n", xlab = "", ylab = "", cex.axis = 1.25)
polygon(xr, dens, lwd = 2, col = adjustcolor(1,0.2))

par(mfrow = c(3,5), mar = c(2,1,0,1))
for(i in 1:15){
  dens <- dbeta(xr, sample(5:15, 1), sample(10:30, 1))
  plot(xr, dens, type = "l", lwd = 2, frame = F, yaxt = "n", xlab = "", ylab = "", cex.axis = 1.25)
  polygon(xr, dens, lwd = 2, col = adjustcolor(tissue_cols[rownames(table_to_use)[i]],0.2))
  
}

par(mfrow = c(10,20), mar = c(0.1,0.1,0.1,0.1))
for(i in 1:200){
  dens <- dbeta(xr, sample(30:40, 1), sample(80:120, 1))
  plot(xr, dens, type = "l", lwd = 0.5, frame = F, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  polygon(xr, dens, lwd = 0.5, col = adjustcolor(tissue_cols[sample(rownames(table_to_use),1)],1))
  
}

#### compute binomial model counts to test for enrichment in + hits ####
twas_with_hits <- colnames(prop_degs_are_twas)
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]

#get "null hypothesis" genes for comparison group
load('~/data/smontgom/transcript_rna_seq_20211008.RData')
rna_dea <- transcript_rna_seq$timewise_dea
rna_dea <- rna_dea[rna_dea$comparison_group == "8w" & rna_dea$selection_fdr > 0.05,]
rna_dea_null <- lapply(setNames(unique(rna_dea$tissue_abbreviation), unique(rna_dea$tissue_abbreviation)), function(tiss) {
  print(tiss)
  subset <- rna_dea[rna_dea$tissue_abbreviation == tiss,]
  genes_in_sub <- table(subset$feature_ID)
  genes_in_sub <- names(genes_in_sub)[genes_in_sub == 2]
  if(length(genes_in_sub) == 0){return(NULL)}
  subset <- subset[subset$feature_ID %in% genes_in_sub,c("feature_ID", "sex", "logFC")]
  subset <- tidyr::spread(subset, key = "sex", value = "logFC")
  subset <- subset[(sign(subset$female) * sign(subset$male)) == 1,]
  x <- sign(subset$female)
  names(x) <- subset$feature_ID
  return(x)
})
do.call(rbind, lapply(rna_dea_null, function(x) table(x)))
rm(transcript_rna_seq)


data <- list()
for(tissue_abbrev in names(twas_intersect)){
  
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){next()}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  nDE_genes_signs <- rna_dea_null[[tissue_abbrev]]
  names(nDE_genes_signs) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(names(nDE_genes_signs), gene_map$RAT_ENSEMBL_ID)]
  nDE_genes_signs <- nDE_genes_signs[!is.na(names(nDE_genes_signs))]
  nDE_genes_signs <- nDE_genes_signs[!is.na(nDE_genes_signs)]
  
  results_1 <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  results_2 <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(nDE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * nDE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  results_1$tissue <- tissue_abbrev
  results_1$trait <- rownames(results_1)
  results_1$group <- "focal"
  
  results_2$tissue <- tissue_abbrev
  results_2$trait <- rownames(results_2)
  results_2$group <- "compl"
  
  data[[tissue_abbrev]] <- rbind(results_1, results_2)
  
}

hist(log2(unlist(lapply(data, function(x) (x$count / x$total)[x$group == "focal"] / (x$count / x$total)[x$group == "compl"]))), breaks = c(-3,0, 3))
x <- log2(unlist(lapply(data, function(x) (x$count / x$total)[x$group == "focal"] / (x$count / x$total)[x$group == "compl"])))
x <- x[x != Inf & x != -Inf & !is.na(x)]
hist(x, breaks = c(-max(abs(x)),0,max(abs(x))))
mean(x > 0)

#### compute binomial model counts to test for enrichment in + hits, alternate version ####
twas_with_hits <- colnames(prop_degs_are_twas)
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
possible_genes <- intersect(all_orthologs_tested, all_twas_genes_tested)
compatible_twas_genes <- some.twas$gene_name %in% possible_genes
twas_by_tissue <- lapply(setNames(unique(some.twas$tissue), tissue_abbr[unique(some.twas$tissue)]), function(tiss) {
  print(tiss)
  some.twas[some.twas$tissue == tiss & compatible_twas_genes,]
})

#get "null hypothesis" genes for overall probabilities
load('~/data/smontgom/transcript_rna_seq_20211008.RData')
rna_dea <- transcript_rna_seq$timewise_dea
rna_dea <- rna_dea[rna_dea$comparison_group == "8w",]
rna_dea_null <- lapply(setNames(unique(rna_dea$tissue_abbreviation), unique(rna_dea$tissue_abbreviation)), function(tiss) {
  print(tiss)
  subset <- rna_dea[rna_dea$tissue_abbreviation == tiss,]
  subset$gene <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  subset <-subset[subset$gene %in% possible_genes,]
  genes_in_sub <- table(subset$feature_ID)
  genes_in_sub <- names(genes_in_sub)[genes_in_sub == 2]
  if(length(genes_in_sub) == 0){return(NULL)}
  subset <- subset[subset$feature_ID %in% genes_in_sub,c("feature_ID", "sex", "logFC")]
  subset <- tidyr::spread(subset, key = "sex", value = "logFC")
  subset <- subset[(sign(subset$female) * sign(subset$male)) == 1,]
  x <- sign(subset$female)
  names(x) <- gene_map$HUMAN_ORTHOLOG_SYMBOL[match(subset$feature_ID, gene_map$RAT_ENSEMBL_ID)]
  return(x)
})
do.call(rbind, lapply(rna_dea_null, function(x) table(x)))
rm(transcript_rna_seq)

tissues <- tissue_code$abbreviation[match(names(motrpac_gtex_map), tissue_code$tissue_name_release)]
tissues <- setdiff(tissues, c("HYPOTH", "TESTES", "OVARY"))
data_null <- lapply(setNames(tissues, tissues), function(tiss){
  print(tiss)
  tws <- twas_by_tissue[[tiss]]
  motr <- rna_dea_null[[tiss]]
  if(is.null(motr) | is.null(tws)){
    return(NULL)
  }
  compatible_genes <- intersect(names(motr), tws$gene_name)
  motr <- motr[compatible_genes]
  tws <- tws[tws$gene_name %in% compatible_genes,]
  tws$motrpac_sign <- motr[tws$gene_name]
  tws$twas_sign <- sign(tws$zscore)
  tws <- tws[tws$twas_sign != 0,]
  output <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    trait_subset <- tws[tws$trait == trait_i, c("motrpac_sign", "twas_sign")]
    c(count_twas = sum(trait_subset$twas_sign > 0), count_motr = sum(trait_subset$motrpac_sign > 0), total = nrow(trait_subset))
  })))
  output$trait <- rownames(output)
  output$tissue <- tiss
  rownames(output) <- NULL
  output
})
data_null <- do.call(rbind, data_null)
hist(data_null$count_twas / data_null$total)
hist(data_null$count_motr / data_null$total)

adj_pvalue_alpha <- 0.05
sig_twas_by_tissue <- lapply(setNames(names(motrpac_gtex_map), names(motrpac_gtex_map)), 
                             function(tissue) some.twas[some.twas$tissue == tissue & some.twas$adj_pvalue < adj_pvalue_alpha,])
data_not_null <- lapply(setNames(tissues, tissues), function(tissue_abbrev){
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){return(NULL)}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  output <- as.data.frame(t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      output <- c(count = sum(shared_genes_signs > 0), total = length(shared_genes_signs))
      return(output)
    } else {return(c(0,0))}
  })))
  
  output$trait <- rownames(output)
  output$tissue <- tissue_abbrev
  rownames(output) <- NULL
  output
})
data_not_null <- do.call(rbind, data_not_null)

#confirm matchup
all(paste0(data_null$trait, "~", data_null$tissue) == paste0(data_not_null$trait, "~", data_not_null$tissue))

#some quick eda
plot(data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
       (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total),
     (data_not_null$count+1) / (data_not_null$total+2), pch = 19, col = adjustcolor(1, 0.5),
     cex = data_not_null$total / max(data_not_null$total) * 5,
     xlab = "expected proportion positives under null",
     ylab = "posterior mean proportion updated from Beta(1,1)")
# points((data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
#        (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total))[c(3, 403, 563, 803)],
#      (data_not_null$count / data_not_null$total)[c(3, 403, 563, 803)], pch = 19, col = adjustcolor(3, 0.5),
#      cex = (data_not_null$total / max(data_not_null$total) * 5)[c(3, 403, 563, 803)])
abline(h = 0.5, lty = 2, col = adjustcolor(2, 0.75), lwd = 2)
abline(v = 0.5, lty = 2, col = adjustcolor(2, 0.75), lwd = 2)

#check iid flat beta model
basic_posteriors_masses <- 1 - pbeta(q = data_null$count_twas / data_null$total * data_null$count_motr / data_null$total + 
        (1-data_null$count_twas / data_null$total) * (1-data_null$count_motr / data_null$total),
      shape1 = 1 + data_not_null$count, shape2 = 1 + data_not_null$total - data_not_null$count)
sum(basic_posteriors_masses > 0.95)
sum(basic_posteriors_masses < 0.05)
head(data_not_null[order(abs(basic_posteriors_masses - 0.5), decreasing = T),], 20)
table(data_not_null[which(basic_posteriors_masses > 0.95),"trait"])
table(data_not_null[which(basic_posteriors_masses < 0.05),"trait"])

hist(sapply(tissues, function(tiss) diff(range(data_null$count_twas[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
hist(sapply(traits, function(trait) diff(range(data_null$count_twas[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))
hist(sapply(tissues, function(tiss) diff(range(data_null$count_motr[data_null$tissue == tiss] / data_null$total[data_null$tissue == tiss]))))
hist(sapply(traits, function(trait) diff(range(data_null$count_motr[data_null$trait == trait] / data_null$total[data_null$trait == trait]))))
trait
data_null$count_motr[data_null$trait == trait] / data_null$total[data_null$trait == trait]
data_null[data_null$trait == trait,]$count_motr / data_null[data_null$trait == trait,]$total
#potential source of heterogeneity here! 


#bayesian model
traits <- unique(data_null$trait)
tissues <- unique(data_null$tissue)
d <- list(intersect_count = data_not_null$count,
          intersect_total = data_not_null$total,
          twas_count = data_null$count_twas,
          motr_count = data_null$count_motr,
          twas_motr_total = data_null$total,
          trait = match(data_null$trait, traits),
          tissue = match(data_null$tissue, tissues),
          n_trait = length(traits),
          n_tissue = length(tissues)
)

# stan_program <- '
# data {
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
#     int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
#     int<lower=0> intersect_count[n_trait * n_tissue];
#     int<lower=0> intersect_total[n_trait * n_tissue];
#     int<lower=0> twas_count[n_trait * n_tissue];
#     int<lower=0> motr_count[n_trait * n_tissue];
#     int<lower=0> twas_motr_total[n_trait * n_tissue];
# }
# transformed data {
#     int<lower=1> n = n_trait * n_tissue;
# }
# parameters {
#     //main parameters
#     real trait_mean;
#     real<lower=0> trait_sd;
#     vector[n_trait] raw_trait_logodds;
# 
#     real tissue_mean;
#     real<lower=0> tissue_sd;
#     vector[n_tissue] raw_tissue_logodds;
#     
#     vector[n] raw_intersect_logodds;
#     real<lower=0> intersect_sd;
#     
#     //biases in deviations terms
#     vector[n_tissue] tissue_bias_mean;
#     vector<lower=0>[n_tissue] tissue_bias_sd;
#     vector[n_trait] trait_bias_mean;
#     vector<lower=0>[n_trait] trait_bias_sd;
#     
#     //multilevel deviation term params
#     real<lower=0> tissue_bias_mean_sd;
#     real<lower=0> trait_bias_mean_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[n_trait] trait_logodds = raw_trait_logodds * trait_sd + trait_mean;
#     vector[n_tissue] tissue_logodds = raw_tissue_logodds * tissue_sd + tissue_mean;
#     
#     vector[n_trait] neg_trait_logodds = logit(1 - inv_logit(trait_logodds));
#     vector[n_tissue] neg_tissue_logodds = logit(1 - inv_logit(tissue_logodds));
#     
#     vector[n] intersect_mean = logit(inv_logit(trait_logodds[trait]) .* inv_logit(tissue_logodds[tissue]) +
#                                      inv_logit(neg_trait_logodds[trait]) .* inv_logit(neg_tissue_logodds[tissue]));
#     vector[n] intersect_logodds = raw_intersect_logodds * intersect_sd .* tissue_bias_sd[tissue] .* trait_bias_sd[trait] +
#                                   intersect_mean + tissue_bias_mean[tissue] * tissue_bias_mean_sd + 
#                                   trait_bias_mean[trait] * trait_bias_mean_sd;
# }
# model {
#     //priors / hyperpriors
#     raw_trait_logodds ~ std_normal();
#     trait_mean ~ normal(0,2);
#     trait_sd ~ normal(0,2);
#     
#     raw_tissue_logodds ~ std_normal();
#     tissue_mean ~ normal(0,2);
#     tissue_sd ~ normal(0,2);
#     
#     raw_intersect_logodds ~ std_normal();
#     intersect_sd ~ normal(0,2);
#     
#     tissue_bias_mean ~ std_normal();
#     tissue_bias_sd ~ std_normal();
#     trait_bias_mean ~ std_normal();
#     trait_bias_sd ~ std_normal();
#     tissue_bias_mean_sd ~ std_normal();
#     trait_bias_mean_sd ~ std_normal();
#     
#     //likelihood
#     twas_count ~ binomial_logit(twas_motr_total, trait_logodds[trait]);
#     motr_count ~ binomial_logit(twas_motr_total, tissue_logodds[tissue]);
#     intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
# }
# generated quantities {
#     vector[n] difference_in_props = inv_logit(intersect_logodds) - inv_logit(intersect_mean);
# }
# '

# stan_program <- '
# data {
#     int<lower=1> n_trait;
#     int<lower=1> n_tissue;
#     int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
#     int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
#     int<lower=0> intersect_count[n_trait * n_tissue];
#     int<lower=0> intersect_total[n_trait * n_tissue];
#     int<lower=0> twas_count[n_trait * n_tissue];
#     int<lower=0> motr_count[n_trait * n_tissue];
#     int<lower=0> twas_motr_total[n_trait * n_tissue];
# }
# transformed data {
#     int<lower=1> n = n_trait * n_tissue;
# }
# parameters {
#     vector[n] raw_intersect_logodds;
#     real intersect_logodds_log_sd;
# }
# transformed parameters {
#     vector[n] intersect_logodds = raw_intersect_logodds * exp(intersect_logodds_log_sd);
# 
# }
# model {
#     raw_intersect_logodds ~ std_normal();
#     intersect_logodds_log_sd ~ normal(0,2);
#     intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
# }
# generated quantities {
#     vector[n] difference_in_props = inv_logit(intersect_logodds) - 0.5;
# }
# '

stan_program <- '
data {
    int<lower=1> n_trait;
    int<lower=1> n_tissue;
    int<lower=1,upper=n_trait> trait[n_trait * n_tissue];
    int<lower=1,upper=n_tissue> tissue[n_trait * n_tissue];
    int<lower=0> intersect_count[n_trait * n_tissue];
    int<lower=0> intersect_total[n_trait * n_tissue];
    int<lower=0> twas_count[n_trait * n_tissue];
    int<lower=0> motr_count[n_trait * n_tissue];
    int<lower=0> twas_motr_total[n_trait * n_tissue];
}
transformed data {
    int<lower=1> n = n_trait * n_tissue;
}
parameters {
    //main parameters
    real trait_mean;
    real<lower=0> trait_sd;
    vector[n_trait] raw_trait_logodds;
    real<lower=0> within_trait_sd;
    vector[n] raw_within_trait_logodds;

    real tissue_mean;
    real<lower=0> tissue_sd;
    vector[n_tissue] raw_tissue_logodds;
    real<lower=0> within_tissue_sd;
    vector[n] raw_within_tissue_logodds;

    vector[n] raw_intersect_logodds;
    real<lower=0> intersect_sd;
    
    //biases in deviations terms
    vector[n_tissue] tissue_bias_mean;
    vector[n_tissue] tissue_bias_log_sd;
    vector[n_trait] trait_bias_mean;
    vector[n_trait] trait_bias_log_sd;
    
    //multilevel deviation term params
    real<lower=0> tissue_bias_log_sd_sd;
    real<lower=0> trait_bias_log_sd_sd;
    real<lower=0> tissue_bias_mean_sd;
    real<lower=0> trait_bias_mean_sd;
}
transformed parameters {
    //recenter params
    vector[n_trait] trait_logodds = raw_trait_logodds * trait_sd + trait_mean;
    vector[n_tissue] tissue_logodds = raw_tissue_logodds * tissue_sd + tissue_mean;
    
    vector[n] within_trait_logodds = raw_within_trait_logodds * within_trait_sd + trait_logodds[trait];
    vector[n] within_tissue_logodds = raw_within_tissue_logodds * within_tissue_sd + tissue_logodds[tissue];
    
    vector<lower=0, upper=1>[n] within_trait_prob = inv_logit(within_trait_logodds);
    vector<lower=0, upper=1>[n] within_tissue_prob = inv_logit(within_tissue_logodds);
    
    vector[n] intersect_mean = logit(within_trait_prob .* within_tissue_prob +
                                     (1-within_trait_prob) .* (1-within_tissue_prob));
    vector[n] intersect_logodds = raw_intersect_logodds * intersect_sd .* exp(tissue_bias_log_sd[tissue] * tissue_bias_log_sd_sd) .* 
                                  exp(trait_bias_log_sd[trait] * trait_bias_log_sd_sd) +
                                  intersect_mean + tissue_bias_mean[tissue] * tissue_bias_mean_sd + 
                                  trait_bias_mean[trait] * trait_bias_mean_sd;
}
model {
    //priors and hyperpriors
    raw_trait_logodds ~ std_normal();
    trait_mean ~ normal(0,2);
    trait_sd ~ normal(0,2);
    raw_within_trait_logodds ~ std_normal();
    within_trait_sd ~ normal(0,2);
    
    raw_tissue_logodds ~ std_normal();
    tissue_mean ~ normal(0,2);
    tissue_sd ~ normal(0,2);
    raw_within_tissue_logodds ~ std_normal();
    within_tissue_sd ~ normal(0,2);
    
    raw_intersect_logodds ~ std_normal();
    intersect_sd ~ normal(0,2);
    
    tissue_bias_mean ~ std_normal();
    tissue_bias_log_sd ~ std_normal();
    tissue_bias_log_sd_sd ~ std_normal();
    trait_bias_mean ~ std_normal();
    trait_bias_log_sd ~ std_normal();
    trait_bias_log_sd_sd ~ std_normal();
    tissue_bias_mean_sd ~ std_normal();
    trait_bias_mean_sd ~ std_normal();
    
    //likelihood
    twas_count ~ binomial_logit(twas_motr_total, within_trait_logodds);
    motr_count ~ binomial_logit(twas_motr_total, within_tissue_logodds);
    intersect_count ~ binomial_logit(intersect_total, intersect_logodds);
}
generated quantities {
    vector[n] difference_in_props = inv_logit(intersect_logodds) - inv_logit(intersect_mean);
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)


#fit model
# write_stan_file(stan_program, dir = "~/Desktop/", basename = "deviation_from_expected_prop_pos")
# write_stan_json(d, "~/Desktop/deviation_from_expected_prop_pos.json")
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                  adapt_delta = 0.9, refresh = 50, init = 0.1, max_treedepth = 20, thin = 5)
# out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.9, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]
samps <- data.frame(as_draws_df(out$draws()))

#take a look at output
hist(apply(samps[,grep("difference_in_props", colnames(samps))], 2, mean))
prop_greater_than_0 <- function(x) mean(x>0)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)
sum((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) < 0.05)
hist(invlogit(samps$intersect_logodds.2.) - invlogit(samps$intersect_mean.2.))
data_not_null[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95),]
hist(invlogit(samps$intersect_mean.1.))
hist(invlogit(samps$intersect_logodds.1.))

#compare estimated mean probs t

sum((apply(samps[,grep("tissue_bias_mean\\.", colnames(samps))], 2, prop_greater_than_0)) > 0.95)
traits[which((apply(samps[,grep("trait_bias_mean\\.", colnames(samps))], 2, prop_greater_than_0)) > 0.95)]

plot(data_null$count_twas / data_null$total,
     apply(invlogit(samps[,setdiff(grep("within_trait_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean),
     pch = 19, col = adjustcolor(match(data_null$trait, traits), 0.75))
abline(0,1,col=2,lty=2)

plot(data_null$count_motr / data_null$total,
     apply(invlogit(samps[,setdiff(grep("within_tissue_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
abline(0,1,col=2,lty=2)

plot((data_not_null$count + 1) / (data_not_null$total + 2), 
     apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean),
     col = adjustcolor(1, 0.25), pch = 19)
points((data_not_null$count / data_not_null$total)[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)], 
     (apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))[which((apply(samps[,grep("difference_in_props", colnames(samps))], 2, prop_greater_than_0)) > 0.95)],
     col = adjustcolor(3, 1), pch = 19)
abline(0,1,col=2,lty=2)

hist(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
hist(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))
data_not_null[which(abs(apply(invlogit(samps[,setdiff(grep("intersect_logodds", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean)-0.5) > 0.1),]

hist(apply(invlogit(samps[,setdiff(grep("intersect_mean", colnames(samps)), grep("raw", colnames(samps)))]), 2, mean))



#### make a summary table for Nicole ####
twas_intersect <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(x) NULL)
twas_with_hits <- colnames(prop_degs_are_twas)
twas_with_hits <- colnames(prop_degs_are_twas)
prop_pos <- lapply(setNames(rownames(n_deg_sigtwas_intersect), rownames(n_deg_sigtwas_intersect)), function(i) NULL)
gene_map <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
tissue_code <- tissue_code[,c("tissue_name_release", "abbreviation")]
tissue_code <- tissue_code[tissue_code$tissue_name_release != "" & !is.na(tissue_code$abbreviation)]
for(tissue_abbrev in names(twas_intersect)){
  
  tissue_i <- tissue_code$tissue_name_release[tissue_code$abbreviation == tissue_abbrev]
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  timepoint = "8w"
  sex = "male"
  
  DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
  
  if(length(DE_genes_in_Nodes) == 0){next()}
  
  DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
    node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
    gene = DE_genes_in_Nodes)
  
  colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
  DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
  DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
  DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
  DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
  
  DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
  names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
  
  results <- lapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
    twas_signs <- sig_twas_by_trait_signs[[trait_i]]
    shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
    if(length(shared_genes) != 0){
      shared_genes_signs <- twas_signs[shared_genes] * DE_genes_signs[shared_genes]
      shared_genes_rat <- gene_map$RAT_SYMBOL[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      shared_genes_rat_ensembl_id <- gene_map$RAT_ENSEMBL_ID[match(shared_genes, gene_map$HUMAN_ORTHOLOG_SYMBOL)]
      
      enrichment_test <- signif_df[which(signif_df$tissue == tissue_abbrev & signif_df$trait == trait_i),]
      
      output <- data.frame(tissue = tissue_abbrev,
                           trait = trait_i,
                           human_ortholog_symbol = shared_genes,
                           rat_symbol = shared_genes_rat,
                           rat_ensembl_ID = shared_genes_rat_ensembl_id,
                           direction_of_DE_on_trait = shared_genes_signs,
                           Pr_Geneset_DiffProp_is_Pos = enrichment_test$prob_diff_is_positive,
                           geneset_enrichment = enrichment_test$signif)
      return(output)
    } else {return(NULL)}
  })
  
  results <- do.call(rbind, results)
  
  twas_intersect[[tissue_abbrev]] <- results
  
}

twas_intersect <- do.call(rbind, twas_intersect)
rownames(twas_intersect) <- NULL
twas_intersect <- as.data.table(twas_intersect)
fwrite(x = twas_intersect, file = "~/data/twas_intersect_table_uncalibrated_Pr.txt")
fread(file = "~/data/twas_intersect_table_uncalibrated_Pr.txt")

#### now plot the proportion of DEGs & TWAS hits in + & - directions ####
twas_with_hits <- colnames(prop_degs_are_twas)
deg_sigtwas_proportion <- array(NA, dim = c(length(motrpac_gtex_map), length(twas_with_hits), 4, 3, 2), 
                                dimnames = list(names(motrpac_gtex_map), twas_with_hits, paste0(2^(0:3), "w"), c("p", "n", "genes"), c("male", "female")))
nuniq <- function(x) length(unique(x))
tissue_code <- MotrpacBicQC::bic_animal_tissue_code
for(tissue_i in names(motrpac_gtex_map)){
  
  print(tissue_i)
  
  sig_twas_by_trait <- lapply(setNames(twas_with_hits, twas_with_hits), 
                              function(trait_i) sig_twas_by_tissue[[tissue_i]][sig_twas_by_tissue[[tissue_i]]$trait == trait_i,])
  sig_twas_by_trait_genes <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i) sig_twas_by_trait[[trait_i]]$gene_name)
  sig_twas_by_trait_signs <- lapply(setNames(twas_with_hits, twas_with_hits), 
                                    function(trait_i){
                                      signs <- sign(sig_twas_by_trait[[trait_i]]$zscore)
                                      names(signs) <- sig_twas_by_trait_genes[[trait_i]]
                                      signs
                                    })
  
  for(timepoint in paste0(2^(0:3), "w")){
    
    for(sex in c("male", "female")){
      
      DE_genes_in_Nodes <- node_metadata_list[[timepoint]]$human_gene_symbol[
        node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]]
      if(length(DE_genes_in_Nodes) == 0){next()}
      DE_genes_in_Nodes_sign <- cbind(as.data.frame(do.call(rbind, strsplit(node_metadata_list[[timepoint]]$node[
        node_metadata_list[[timepoint]]$tissue == tissue_code$abbreviation[tissue_code$tissue_name_release == tissue_i]], split = "_"))),
        gene = DE_genes_in_Nodes)
      colnames(DE_genes_in_Nodes_sign) <- c("time", "female", "male", "gene")
      DE_genes_in_Nodes_sign$female[grep(pattern = "F1", DE_genes_in_Nodes_sign$female)] <- 1
      DE_genes_in_Nodes_sign$female[grep(pattern = "F-1", DE_genes_in_Nodes_sign$female)] <- -1
      DE_genes_in_Nodes_sign$male[grep(pattern = "M1", DE_genes_in_Nodes_sign$male)] <- 1
      DE_genes_in_Nodes_sign$male[grep(pattern = "M-1", DE_genes_in_Nodes_sign$male)] <- -1
      DE_genes_in_Nodes_sign <- DE_genes_in_Nodes_sign[!is.na(DE_genes_in_Nodes_sign$gene),]
      
      DE_genes_signs <- as.integer(DE_genes_in_Nodes_sign[,sex])
      names(DE_genes_signs) <- DE_genes_in_Nodes_sign$gene
      
      results <- t(sapply(setNames(twas_with_hits, twas_with_hits), function(trait_i) {
        twas_signs <- sig_twas_by_trait_signs[[trait_i]]
        shared_genes <- intersect(names(DE_genes_signs), names(twas_signs))
        prop_positive <- mean(twas_signs[shared_genes] * DE_genes_signs[shared_genes] > 0)
        return(list(prop_positive = prop_positive, n = length(shared_genes), genes = paste0(shared_genes, collapse = " ~ ")))
      }))
      
      # DELT <- as.data.frame(deg_eqtl_list_TWAS_cluster_subset[[tissue_i]])
      # 
      # #subset to time, sex, and unique gene entries
      # DELT <- DELT[DELT$comparison_group == time,]
      # DELT <- DELT[DELT$sex == sex,]
      # # DELT <- DELT[-which(is.na(DELT$gene_name)),] #get rid of na genes
      # if(nrow(DELT) == 0){next()}
      # DELT <- DELT[match(unique(DELT$gene_name), DELT$gene_name),] #pull out only first gene entry
      # 
      # 
      # DE_inds <- which(DELT$adj_p_value <= 1.1)
      # # TWAS_inds <- apply(DELT[,grep(colnames(DELT), pattern = "BH_PValue")] < 0.05, 2, which)
      # TWAS_inds <- apply(log(DELT[,grep(colnames(DELT), pattern = ".pvalue")]) <= (log(0.05) - log(n_twas_comparisons_clusters)), 2, which)
      # if(length(TWAS_inds) == 0){next()}
      # 
      # intersect_inds <- lapply(TWAS_inds, function(twi) intersect(DE_inds, twi))
      # # intersect_inds <- lapply(intersect_inds, function(ii) ii[DELT$gene_name[ii] %in% cluster_genes]) #subset to just monotonic sex-homogenous clusters
      # intersect_genes <- lapply(intersect_inds, function(ii) unique(DELT$gene_name[ii]))
      # 
      # 
      # intersect_sign_logFC <- lapply(intersect_inds, function(ii) sign(DELT$logFC[ii]))
      # names(intersect_sign_logFC) <- gsub(names(intersect_sign_logFC), pattern = ".pvalue", replacement = "")
      # intersect_sign_TWAS <- lapply(gsub(names(intersect_inds), pattern = ".pvalue", replacement = ""), function(ii) 
      #                        sign(DELT[intersect_inds[paste0(ii, ".pvalue")][[1]], paste0(ii, ".zscore")]))
      # intersect_sign_match <- lapply(1:length(intersect_sign_logFC), function(i) intersect_sign_logFC[[i]]*intersect_sign_TWAS[[i]])
      # 
      # intersect_genes_strings <- lapply(1:length(intersect_sign_logFC), function(i) 
      #   paste0(intersect_genes[[i]], " (", c("-", "", "+")[intersect_sign_match[[i]] + 2], ")", collapse = " ~ "))
      # intersect_genes_strings <- unlist(intersect_genes_strings)
      # intersect_genes_strings[intersect_genes_strings == " ()"] <- ""
      # 
      # intersect_sign_match <- sapply(intersect_sign_match, function(x) mean(x == 1))
      # intersect_sign_match_n <- sapply(intersect_inds, function(x) length(x))
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "p", sex] <- intersect_sign_match
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "n", sex] <- intersect_sign_match_n
      # deg_sigtwas_proportion[tissue_i,names(n_intersect),time, "genes", sex] <- intersect_genes_strings
      
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "p", sex] <- unlist(results[,"prop_positive"])
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "n", sex] <- unlist(results[,"n"])
      deg_sigtwas_proportion[tissue_i,rownames(results),timepoint, "genes", sex] <- unlist(results[,"genes"])
      
    }
    
  }
  
}


#### the actual plotting ####
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "type_1_diabetes", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "reported_hypertension", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "heart_problems", ignore.case = T)],,"n","male"]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_mass_index_BMI", ignore.case = T)],,"n",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "UKB_20002_1462_self_reported_crohns_disease", ignore.case = T)],,"p",sex]
deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "tag.evrsmk.tbl", ignore.case = T)],,"p",sex]

deg_sigtwas_proportion[,twas_with_hits[grep(twas_with_hits, pattern = "UKB_50_Standing_height", ignore.case = T)],,"n",sex]

#plot lines for proportions
trait <- twas_with_hits[grep(twas_with_hits, pattern = "Body_fat_percentage", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "reported_hypertension", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "heart_problems", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "Body_mass_index_BMI", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "cholesterol", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "ldl", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "CAD", ignore.case = T)][1]
trait <- twas_with_hits[grep(twas_with_hits, pattern = "alzhei", ignore.case = T)][1]

cols = list(Tissue=tissue_cols[names(motrpac_gtex_map)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

#plotting params
plot_gene_names <- F

cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_in_",trait,".pdf"), 
          width = 1400 / 72, height = 500 / 72, family="Arial Unicode MS")
par(mfrow = c(1,2), xpd = NA, mar = c(5,5,5,18.5))
lwd <- 3
tissue_name_cex <- 1.2
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
                                                            c(d[,"8w"] * 10 + rnorm(nrow(d), 0, 0.0001))),
                                             rep.fact = 20, adj.max = 20)
  xylocs_tissue_names$y <- xylocs_tissue_names$y / 100
  xylocs_tissue_names <- as.data.frame(cbind(x = rep(4 + 0.25, nrow(d)), y = c(d[,"8w"])))
  
  #find coordinates to plot gene names
  if(plot_gene_names){
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
  }
  
  
  plot(100,100,xlim = c(1,4), ylim = c(0,1), xpd=NA, 
       ylab = "Proportion Positive Effects on GWAS Trait", xlab = "Timepoint", xaxt = "n", 
       yaxt = "n", bty="n", cex.lab = 1.5, cex.axis = 1.25)
  
  
  #plot faded positive and negative regions
  rect(xl = 1, xr = 4, yb = 0.5, ytop = 1,
       col = grDevices::adjustcolor("red", 0.1), border = NA)
  rect(xl = 1, xr = 4, yb = 0, ytop = 0.5,
       col = grDevices::adjustcolor("blue", 0.1), border = NA)
  
  text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 3, pos = 3,
       x = par("usr")[2] * 0.5 + par("usr")[1] * 0.5, y = par("usr")[3] * 0 + par("usr")[4] * 1)
  
  #horizontal axis
  segments(x0 = 1:4, x1 = 1:4, y0 = - 0.02, y1 = - 0.04, lwd = 2)
  segments(x0 = 1, x1 = 4, y0 = - 0.02, y1 = - 0.02, lwd = 2)
  text(x = 1:4, y = - 0.07, labels = paste0(2^(0:3), "w"), pos = , cex = 1.251)
  
  #vertical axis
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
  text(x = 1 - 3 * 0.035, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.25)
  
  
  for(tissue in rownames(d)){
    
    if(all(is.na(d[tissue,]))){
      next()
    }
    
    lines(1:4, d[tissue,], lwd = lwd, col = cols$Tissue[tissue])
    
    #plot tissue names
    text(x = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue], cex = tissue_name_cex,
         y = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue], 
         labels = paste0(tissue, " (", dn[tissue,"8w"], ")"), col = cols$Tissue[tissue], pos = 4)
    segments(x0 = 4, y0 = d[tissue,"8w"], 
             x1 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 0.075, 
             y1 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
             col = cols$Tissue[tissue], lty = 3)
    
    #plot gene names
    if(plot_gene_names){
      tissue_genes <- xylocs_tissues_genes[xylocs_tissues_genes$tissue == tissue,]
      text(labels = tissue_genes$gene,
           x = tissue_genes$x,
           y = tissue_genes$y,
           cex = 0.875, col = cols$Tissue[tissue], pos = 4)
      for(gene_i in 1:nrow(tissue_genes)){
        segments(x0 = xylocs_tissue_names$x[rownames(xylocs_tissue_names) == tissue] + 
                   strwidth(paste0(tissue, " (", dn[tissue,"8w"], ")   "), cex = tissue_name_cex, units = "user"), 
                 y0 = xylocs_tissue_names$y[rownames(xylocs_tissue_names) == tissue],
                 x1 = tissue_genes$x[gene_i] + 0.75*strwidth("  ", cex = tissue_name_cex, units = "user"),
                 y1 = tissue_genes$y[gene_i],
                 col = cols$Tissue[tissue], lty = 3)
      }
    }
    
  }
  # legend(x = 1, y = 1.5, legend = tissue_names, 
  #        col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
  # segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")
  
}
text(trait_categories$new_Phenotype[trait_categories$Tag == trait], x =0.3, y = 1.175, cex = 2, font = 2)
dev.off()


#### plot scatterplot ####

cols = list(Tissue=tissue_cols[names(motrpac_gtex_map)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue<- cols$Tissue[order(match(MotrpacBicQC::tissue_abbr[names(cols$Tissue)], MotrpacBicQC::tissue_order))]

tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

#plotting params
trait_category_legend_below <- T
deg_sigtwas_proportion[,order(),,,]

cairo_pdf(paste0("~/Documents/Documents - nikolai/pass1b_fig8_DE_protective_effect_scatterplot.pdf"), 
          width = 2000 / 72, height = 900 / 72, family="Arial Unicode MS")
par(mfrow = c(1,1), xpd = NA, mar = c(14,9,5,7))
lwd <- 3
tissue_name_cex <- 1.2

plot(100,100,xlim = c(3,length(twas_with_hits)), ylim = c(0,1), xpd=NA, 
     ylab = "Proportion Positive Effects on GWAS Trait", xlab = "", xaxt = "n", 
     yaxt = "n", bty="n", cex.lab = 2, cex.axis = 1.25)

for(sex in c("male", "female")[1]){
  
  # mean_freq <- apply(sapply(twas_with_hits, function(trait) as.numeric(deg_sigtwas_proportion[, trait,"8w","p",sex])), 2, weighted.mean, na.rm = T)
  
  mean_freq <- sapply(setNames(1:length(twas_with_hits), twas_with_hits), function(trait_i){
     weighted.mean(as.numeric(deg_sigtwas_proportion[, twas_with_hits[trait_i],"8w","p",sex]), 
                   w = as.numeric(deg_sigtwas_proportion[, twas_with_hits[trait_i],"8w","n",sex]),
                  na.rm = T)})
  
  order_twas_with_hits <- order(mean_freq, decreasing = T)
  
  #plot faded positive and negative regions
  rect(xl = 1, xr = length(twas_with_hits) + 2, yb = 0.5, ytop = 1,
       col = grDevices::adjustcolor("red", 0.1), border = NA)
  rect(xl = 1, xr = length(twas_with_hits) + 2, yb = 0, ytop = 0.5,
       col = grDevices::adjustcolor("blue", 0.1), border = NA)
  
  #vertical axis
  segments(x0 = 1 - length(twas_with_hits) * 0.005, x1 = 1 - 3 * 0.04, y0 = 0:5/5, y1 = 0:5/5, lwd = 2)
  segments(x0 = 1 - 3 * 0.02, x1 = 1 - 3 * 0.02, y0 = 0, y1 = 1, lwd = 2)
  text(x = 1 - 3 * 0.05, y = 0:5/5, labels = 0:5/5, pos = 2, cex = 1.5)
  
  #horizontal axis
  text(1:length(twas_with_hits) + 1.5, -0.07, srt = 45, pos = 2,
       trait_categories$new_Phenotype[match(twas_with_hits[order_twas_with_hits], trait_categories$Tag)])
  segments(1:length(twas_with_hits) + 1, 0, 1:length(twas_with_hits) + 1, 1, col = adjustcolor(1, 0.2), lty = 2)
  
  #plot points
  points(1:length(twas_with_hits)+1, y = sort(mean_freq, decreasing = T), pch = "*", cex = 3)
  
  for(trait_i in (1:length(twas_with_hits))){
    
    #plot category blocks
    rect(xleft = which(order_twas_with_hits == trait_i) + 1/2, xright = which(order_twas_with_hits == trait_i) + 3/2,
         ybottom = -0.06, ytop = -0.02,
         col = category_colors[trait_categories$Category[match(twas_with_hits[trait_i], trait_categories$Tag)]])
    

    d <- as.numeric(deg_sigtwas_proportion[, twas_with_hits[trait_i],"8w","p",sex])
    dn <- as.numeric(deg_sigtwas_proportion[, twas_with_hits[trait_i],"8w","n",sex])
    names(d) <- names(dn) <- rownames(deg_sigtwas_proportion)
    dn <- dn[!is.na(d)]
    d <- d[!is.na(d)]
    #plot white points first
    points(x = rep(which(order_twas_with_hits == trait_i) + 1, length(d)), y = d, pch = 19, col = "white", cex = dn^0.25)
    #and then the actual points
    points(x = rep(which(order_twas_with_hits == trait_i) + 1, length(d)), y = d, pch = 19, col = adjustcolor(cols$Tissue[names(d)], 0.5), cex = dn^0.25)
    
  }
  
  
  #legend for tissues
  points(x = rep(length(twas_with_hits) + 3, length(cols$Tissue)), 
         y = seq(0.5, 0.99, length.out = length(cols$Tissue)), 
         col = adjustcolor(cols$Tissue, 0.5),
         pch = 19, cex = 2)
  text(x = rep(length(twas_with_hits) + 3.25, length(cols$Tissue)), 
       y = seq(0.5, 0.99, length.out = length(cols$Tissue)), 
       pos = 4, pch = 19, cex = 1,
       labels = MotrpacBicQC::tissue_abbr[names(cols$Tissue)])
  
  points(x = length(twas_with_hits) + 3, 
         y = min(seq(0.5, 0.99, length.out = length(cols$Tissue))) - 0.0675, 
         col = "black", pch = "*", cex = 3)
  text(x = length(twas_with_hits) + 3.25, cex = 1.1,
       y = min(seq(0.5, 0.99, length.out = length(cols$Tissue))) - 0.075, 
       pos = 4, labels = "Weighted\nMean")
  
  #legend for tissue size
  n_pts <- 5
  pt_size_logs <- seq(1, log(max(as.numeric(deg_sigtwas_proportion[,,"8w","n",sex]), na.rm = T)) / log(2), length.out = n_pts)
  pt_size_legend <- round(2^pt_size_logs)
  text(x = length(twas_with_hits) + 2.25, 
       y = 0.1875 + 0.01 * n_pts + sum(pt_size_legend^0.25/100), 
       pos = 4, pch = 19, cex = 1.1,
       labels = "Sample Size")
  points(x = rep(length(twas_with_hits) + 3.25, n_pts), 
         y = 0.15 + cumsum(pt_size_legend^0.25/100) + cumsum(rep(0.01, n_pts)), 
         col = adjustcolor(1, 0.5),
         pch = 19, cex = pt_size_legend^0.25)
  text(x = rep(length(twas_with_hits) + 3.25, n_pts) + pt_size_legend^0.25/10, 
       y = 0.15 + cumsum(pt_size_legend^0.25/100) + cumsum(rep(0.01, n_pts)), 
       pos = 4, pch = 19, cex = 1,
       labels = pt_size_legend)
  
  #legend for categories
  if(trait_category_legend_below){
    for(i in 1:length(category_colors)){
      rect(xleft = 0 + i + cumsum(c(0,strwidth(names(category_colors), cex = 1.25, units = "user") + 1))[i],
           xright = 1 + i + cumsum(c(0,strwidth(names(category_colors), cex = 1.25, units = "user") + 1))[i],
           ybottom = -0.375,
           ytop = -0.335,
           col = category_colors[i], border = 1)
      text(labels = names(category_colors)[i], pos = 4, y = -0.355, cex = 1.25, x = 0.85 + i + cumsum(c(0,strwidth(names(category_colors), cex = 1.25, units = "user") + 1))[i])
      #draw circle?
      arc(t1 = 3*pi/2, t2 = pi/2, r1 = (0.355-0.04) / 2, r2 = (0.355-0.04) / 2, center = c(0.5,(-0.355-0.04)/2), lwd = 3, res = 50, adjx = 30)
      points(0.5, -0.04, pch = -9658, cex = 3)
    }
  } else {
    x_adj2 <- x_adj - 2.5
    y_adj2 <- y_adj - 2.5
    for(i in 1:length(category_colors)){
      rect(xleft = 0 + i,
           xright = 1 + i,
           ybottom = -10,
           ytop = -9,
           col = category_colors[i], border = 1)
      text(labels = names(category_colors)[i], pos = 4, y = i + y_adj2 - 7, x = ncol(table_to_use) - 1/2 + x_adj2 + 1/3)
    }
  }

}

text("Proportion of Positive Effects Across Tissues and Traits at 8W Timepoint", x =40, y = 1.075, cex = 3, font = 2)

dev.off()

#### stitch pdfs togetehr ####

pdftools::pdf_combine(paste0("~/Documents/Documents - nikolai/", 
                             c("pass1b_fig8_DEG-TWAS_Intersect_counts",
                             "pass1b_fig8_DEG-TWAS_Intersect_permille", 
                             "pass1b_fig8_DE_protective_effect_scatterplot"),".pdf"),
                      output = paste0("~/Documents/Documents - nikolai/", ifelse(use_tissue_cols_for_cols, "option_1", "option_2"), ".pdf"))


#####

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


relative_expression_data <- 
  list(phenotypic_expression = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
    sapply(tissues, function(tissue) quantile(
      x = deg_eqtl_list[[tissue]]$phenotypic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
      probs = qs2use, na.rm = T))),
      
      genetic_expression = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
    sapply(tissues, function(tissue) quantile(
      x = deg_eqtl_list[[tissue]]$genetic_expression_Z[deg_eqtl_list[[tissue]]$comparison_group == ti], 
      probs = qs2use, na.rm = T))),
      
      genetic_expression_p2SE = lapply(setNames(paste0(2^(0:3), "w"),paste0(2^(0:3), "w")), function(ti) 
    sapply(tissues, function(tissue) quantile(
      x = (deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z[deg_eqtl_list[[tissue]]$comparison_group == ti]), 
      probs = qs2use, na.rm = T)))
  )
save(file = "~/data/smontgom/relative_expression_motrpac_gtex", relative_expression_data)

tissue <- tissues[3]
plot(deg_eqtl_list[[tissue]]$logFC, deg_eqtl_list[[tissue]]$genetic_expression_plus2SE_Z)

#### just the plotting ####
#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
# squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
# unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
squish_middle_x <- function(x,f) asinh(x*f)
unsquish_middle_x <- function(x,f) sinh(x)/f
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}

load("~/data/smontgom/relative_expression_motrpac_gtex")

EZ_PZ <- relative_expression_data[[1]]

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
