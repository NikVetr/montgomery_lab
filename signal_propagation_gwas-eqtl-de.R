library(ks)
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
library(data.table)

if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
} else {
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}

getDoParWorkers()

library(data.table)


####load deg eqtl data ####
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
# system("gsutil cp gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt ~/data/smontgom/")
GTEx_logo <- readJPEG("~/Documents/GTEx_logo.jpg")
eqtl = '~/data/smontgom/GTEx_Analysis_v8_eQTL'
#deg = 'gs://mawg-data/pass1b-06/transcript-rna-seq/dea/'
deg = '~/data/smontgom/dea/'
#map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t')
map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt', sep='\t', header=T)
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



if(!exists("rna_dea")){
  load('~/data/smontgom/rna_dea_20210114.RData')
  load('~/data/smontgom/dea/transcript_rna_seq_20210126.RData')
  old_rna_dea <- rna_dea
  rna_dea$training_dea <- as.data.table(transcript_rna_seq$training_dea)
  rna_dea$timewise_dea <- as.data.table(transcript_rna_seq$timewise_dea)
  new_rna_dea <- rna_dea
}

if(!exists("deg_eqtl_list")){
  deg_eqtl_list = list()
  for(motrpac_tissue in unique(rna_dea$timewise_dea$tissue)){
    
    if(!motrpac_tissue %in% names(motrpac_gtex_map)){next}
    
    cat(paste0(motrpac_tissue, "\n"))
    
    # read in eQTLs
    gtex_tissue = motrpac_gtex_map[[motrpac_tissue]]
    gtex_egene = fread(sprintf('%s/%s.v8.egenes.txt.gz',eqtl,gtex_tissue), sep='\t', header=T)
    gtex_egene[, human_ensembl_gene := gsub('\\..*','',gene_id)]
    
    # match human genes with rat ensembl genes 
    gtex_egene = merge(gtex_egene, map, by='human_ensembl_gene')
    gtex_motrpac = merge(gtex_egene, rna_dea$timewise_dea[tissue == motrpac_tissue], by='feature_ID')
    gtex_motrpac$abs_slope <- abs(gtex_motrpac$slope)
    
    deg_eqtl_list[[motrpac_tissue]] = gtex_motrpac

  }
}

cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

#DE in rat
deg_eqtl_list[[1]]$logFC

#lead eQTL effect 
deg_eqtl_list[[1]]$slope

#lead eQTL effect variant position
deg_eqtl_list[[1]]$chr
deg_eqtl_list[[1]]$variant_pos

#colocalization
pp4_threshold <- 1e-4
pp4_threshold <- 0
load(file=paste0('~/data/smontgom/coloc_list_pp4threshold_', pp4_threshold,'.RData'))
# coloc_list

#GWAS effect size
ldsc_directory <- "~/repos/ldsc/"
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
# gi = 1
# fread(paste0(ldsc_directory, "gwas_sumstats/proper_format/signed/", c(gwas_summary_files[gi]), ".sumstats.gz"))
# gwas_data <- fread(paste0(gwas_dir, gwas_summary_files[10]))
# sum(!is.na(gwas_data$effect_size)) / length((gwas_data$effect_size))

#now merge everything together

for(gi in 1:length(gwas_summary_files)){
# foreach(gi=1:length(gwas_summary_files), .packages = c("data.table")) %dopar% {
  
  trait <- gsub(gwas_summary_files[gi], pattern = ".txt.gz", replacement = "")
  cat(paste0(trait, "\n"))
  
  gwas_data_full <- fread(paste0(gwas_dir, gwas_summary_files[gi]))
  gwas_data <- gwas_data_full[,c("chromosome", "position", "effect_allele", "non_effect_allele", "zscore", "pvalue", "effect_size", "sample_size", "standard_error", "variant_id", "frequency")]
  colnames(gwas_data) <- paste0("gwas_", colnames(gwas_data))
  gwas_data$chrloc <- paste0(gwas_data$gwas_chromosome, "-", gwas_data$gwas_position)
  
  for(motrpac_tissue in names(deg_eqtl_list)){
    cat(paste0(" ", match(motrpac_tissue, names(deg_eqtl_list))))
  
    motrpac_data <- deg_eqtl_list[[motrpac_tissue]]
    motrpac_data$chrloc <-  paste0(motrpac_data$chr, "-", motrpac_data$variant_pos)
  
    #  match lead eQTL locus with GWAS locus
    merged_data <- merge(motrpac_data, gwas_data, by = "chrloc")
    merged_data$coloc_p4 <- coloc_list_trait_tissue$p4[match(merged_data$gene_id, coloc_list_trait_tissue$gene_id)]
    
    if(!dir.exists("~/data/smontgom/signal_propagation/")){dir.create("~/data/smontgom/signal_propagation")}
    if(!dir.exists(paste0("~/data/smontgom/signal_propagation/", motrpac_tissue))){
      dir.create(paste0("~/data/smontgom/signal_propagation/", motrpac_tissue))}
    fwrite(merged_data, paste0("~/data/smontgom/signal_propagation/", motrpac_tissue, "/", trait, ".txt"))
    
  }
}

#get GWAS sample sizes

#now read it back in and add the colocalization data whoops
for(gi in 1:length(gwas_summary_files)){
# foreach(gi=1:length(gwas_summary_files), .packages = c("data.table")) %dopar% {
    
  trait <- gsub(gwas_summary_files[gi], pattern = ".txt.gz", replacement = "")
  coloc_list_trait <- coloc_list[coloc_list$gwas_trait == gsub("imputed_", "", trait), ]
  
  cat(paste0(trait, "\n"))
  
  for(motrpac_tissue in names(deg_eqtl_list)){
    coloc_list_trait_tissue <- coloc_list_trait[coloc_list_trait$motrpac_tissue == motrpac_tissue, ]  
    data <- fread(paste0("~/data/smontgom/signal_propagation/", motrpac_tissue, "/", trait, ".txt"))
    # if(!("coloc_p4" %in% colnames(data))){
    if(T){
      data$coloc_p4 <- coloc_list_trait_tissue$p4[match(data$gene_id, coloc_list_trait_tissue$gene_id)]
      fwrite(data, paste0("~/data/smontgom/signal_propagation/", motrpac_tissue, "/", trait, ".txt"))
    }
    
  }
}


#### now read it back in and do some basic number crunching ####
selection_alpha <- 0.1
proportion_positive_effects <- array(data = 0, dim = c(length(deg_eqtl_list),length(gwas_summary_files), 4, 3, 2), 
           dimnames = list(tissue = names(deg_eqtl_list), trait = gsub(".txt.gz", "", gsub(pattern = "imputed_", "", gwas_summary_files)), 
                           time = paste0(2^(0:3), "w"), type = c("prop", "est_effect", "n"),sex = c("male", "female")))


for(gi in 1:length(gwas_summary_files)){
  
  trait <- gsub(gwas_summary_files[gi], pattern = ".txt.gz", replacement = "")
  cat(paste0(trait, "\n"))
  
  for(motrpac_tissue in names(deg_eqtl_list)){
    
    data <- fread(paste0("~/data/smontgom/signal_propagation/", motrpac_tissue, "/", trait, ".txt"))
    
    data$alleles_match <- NA
    data$alleles_match[data$alt == data$gwas_effect_allele & data$ref == data$gwas_non_effect_allele] <- 1
    data$alleles_match[data$alt == data$gwas_non_effect_allele & data$ref == data$gwas_effect_allele] <- -1
    # sum(is.na(data$alleles_match))
    
    data$eQTL_DE_same_direction <- sign(data$logFC) * sign(data$slope)
    
    #get standardized effect per allele
    data$gwas_standardized_effect <- tanh(qnorm(1-data$gwas_pvalue/2) / sqrt(data$gwas_sample_size-3)) / sqrt(data$gwas_frequency * (1-data$gwas_frequency) * 2) * sign(data$gwas_zscore)
    # plot(data$gwas_standardized_effect, data$gwas_effect_size)
    
    # data_filtered <- data[!is.na(data$coloc_p4) & data$selection_fdr < selection_alpha,]
    data_filtered <- data[data$selection_fdr < selection_alpha,]
    data_filtered <- data_filtered[!is.na(data_filtered$alleles_match),]
    #var(B_hat) = var(Y) / var(X) * (1-p^2) / (N-3), https://stats.stackexchange.com/questions/102786/sampling-distribution-of-regression-coefficients-for-normally-distributed-random/104821#104821
    
    

    
    if(nrow(data_filtered) > 0){
      
      #compute simple direction of effect
      data_filtered$direction_of_effect <- (data_filtered$eQTL_DE_same_direction * data_filtered$alleles_match * sign(data_filtered$gwas_zscore))
      
      #invert slope in case of mismatched alleles
      data_filtered$slope_matched_alleles <- data_filtered$slope
      data_filtered$slope_matched_alleles[data_filtered$alleles_match < 0] <- -data_filtered$slope_matched_alleles[data_filtered$alleles_match < 0]
      
      #flip gwas standardized effect to get everything on a positive effect size scale
      data_filtered$logFC_strictlyPos <- abs(data_filtered$logFC)
      
      data_filtered$slope_matched_alleles_strictlyPos_intermediate <- data_filtered$slope_matched_alleles
      data_filtered$slope_matched_alleles_strictlyPos_intermediate[data_filtered$logFC < 0] <- -data_filtered$slope_matched_alleles[data_filtered$logFC < 0] #flip sign for when we flipped logFC sign
      data_filtered$slope_matched_alleles_strictlyPos <- abs(data_filtered$slope_matched_alleles_strictlyPos_intermediate)
      
      data_filtered$gwas_standardized_effect_strictlyPos <- data_filtered$gwas_standardized_effect
      data_filtered$gwas_standardized_effect_strictlyPos[data_filtered$slope_matched_alleles_strictlyPos_intermediate < 0] <- 
        -data_filtered$gwas_standardized_effect[data_filtered$slope_matched_alleles_strictlyPos_intermediate < 0] #flip effect sign for when we flipped eQTL sign
      
      data_filtered$slope_matched_alleles_strictlyPos <- abs(data_filtered$slope_matched_alleles_strictlyPos)
      
      data_filtered$est_effect_of_DE <- data_filtered$logFC_strictlyPos / data_filtered$slope_matched_alleles_strictlyPos * data_filtered$gwas_standardized_effect_strictlyPos
      # data_filtered$est_effect_of_DE <- log(2^data_filtered$logFC_strictlyPos) / log(2^data_filtered$slope_matched_alleles_strictlyPos) * data_filtered$gwas_standardized_effect_strictlyPos #equivalently
      # done to fine log_(2^slope){2^logFC}, ie allele equivalents
      
      for(sex in c("male", "female")){
      
        prop_positive_and_expected_effect <- sapply(paste0(2^(0:3), "w"), function(time) 
          c(sum(data_filtered$direction_of_effect[data_filtered$comparison_group == time & data_filtered$sex == sex & data_filtered$gwas_pvalue < 0.05] > 0, na.rm = T) / 
              sum(!is.na(data_filtered$direction_of_effect[data_filtered$comparison_group == time & data_filtered$sex == sex & data_filtered$gwas_pvalue < 0.05])),
          sum(data_filtered$est_effect_of_DE[data_filtered$comparison_group == time & data_filtered$sex == sex], na.rm = T),
            sum(!is.na(data_filtered$direction_of_effect[data_filtered$comparison_group == time & data_filtered$sex == sex]))))
        
        proportion_positive_effects[motrpac_tissue, gsub("imputed_", "", trait),,,sex] <- t(prop_positive_and_expected_effect)
        
      }
    }

  }
}

#### do some plotting ####

traits <- gsub(".txt.gz", "", gsub(pattern = "imputed_", "", gwas_summary_files))
trait <- "CARDIoGRAM_C4D_CAD_ADDITIVE"
trait <-  traits[grep("hypertension", traits)]

tissue_names <- sapply(strsplit(names(cols$Tissue), "-"), function(x) paste0(x[ifelse(length(x) == 2, c(2), list(2:3))[[1]]], collapse = " "))

grDevices::cairo_pdf(filename = paste0("~/Documents/signal_prop/", trait,".pdf"), width = 1000 / 72, height = 500 / 72, family="Arial Unicode MS")


par(xpd = NA, mfrow = c(2,3), mar = c(4,5,2,1))
for(sex in c("male", "female")){

trait_esteffect_range <- quantile(proportion_positive_effects[, trait,,"est_effect",sex], c(0.0, 0.99))
trait_esteffect_range <- range(proportion_positive_effects[, trait,,"est_effect",sex])

  
plot(100,100,xlim = c(1,4), ylim = c(0,1), xpd=NA, 
     ylab = "Proportion Positive Effects on GWAS Scale", xlab = "Timepoint", xaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25
)

text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 2, pos = 3,
     x = par("usr")[2] * 0.95 + par("usr")[1] * 0.05, y = par("usr")[3] * 0.95 + par("usr")[4] * 0.05)
  
segments(x0 = 1:4, x1 = 1:4, y0 = 0, y1 = - 0.02)
segments(x0 = 1, x1 = 4, y0 = 0, y1 = 0)
text(x = 1:4, y = -0.02, labels = paste0(2^(0:3), "w"), pos = 1)
for(i in 1:nrow(proportion_positive_effects[, trait,,"prop",sex])){
  
  if(all(is.na(proportion_positive_effects[i, trait,,"prop",sex]))){
    next()
  }
  
  lines(1:4, proportion_positive_effects[i, trait,,"prop",sex], 
        lwd = 2, col = cols$Tissue[rownames(proportion_positive_effects)[i]])
  
}
legend(x = 1, y = 1.1, legend = tissue_names, 
       col = cols$Tissue, lwd = 3, ncol = 4, cex = 1, border = NA, seg.len = 1, bg = NA, bty = "n", x.intersp = 0.25, text.width = 0.65)
segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")

plot(100,100, xlim = c(1,4), ylim = trait_esteffect_range, xpd=NA, main = ifelse(sex == "male", trait, ""),
     ylab = "Expected Linear Std. Effect on GWAS Trait", xlab = "Timepoint", xaxt = "n", bty="n", cex.lab = 1.25, cex.axis = 1.25,
)

text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 2, pos = 3,
     x = par("usr")[2] * 0.95 + par("usr")[1] * 0.05, y = par("usr")[3] * 0.95 + par("usr")[4] * 0.05)

segments(x0 = 1:4, x1 = 1:4, y0 = trait_esteffect_range[1], y1 = trait_esteffect_range[1]-0.2)
segments(x0 = 1, x1 = 4, y0 = trait_esteffect_range[1], y1 = trait_esteffect_range[1])
text(x = 1:4, y = trait_esteffect_range[1]-0.2, labels = paste0(2^(0:3), "w"), pos = 1)
for(i in 1:nrow(proportion_positive_effects[, trait,,"est_effect",sex])){
  
  if(all(is.na(proportion_positive_effects[i, trait,,"est_effect",sex]))){
    next()
  }
  
  lines(1:4, proportion_positive_effects[i, trait,,"est_effect",sex], 
        lwd = 2, col = cols$Tissue[rownames(proportion_positive_effects)[i]])
  
}
# legend(x = "topright", legend = names(cols$Tissue), col = cols$Tissue, lwd = 3, ncol = 2, cex = 0.95)
segments(x0 = 0.9, y0 = 0, x1 = 4.1, y1 = 0, lty = 3, lwd = 3, col = "lightgrey")

plot(c(proportion_positive_effects[, trait,,"prop",sex]), c(proportion_positive_effects[, trait,,"est_effect",sex]),
     xlab = "proportion of positive effects", ylab = "estimated effect (in units SD)", pch = 19, cex = 1.5, cex.lab = 1.25, cex.axis = 1.25,
     col = adjustcolor(cols$Tissue[rep(rownames(proportion_positive_effects[, trait,,"prop",sex]), 4)], 0.85))
text(c("\u2642", "\u2640")[c("male", "female") == sex], col = cols$Sex[sex], cex = 3.5, font = 2, pos = 3,
     x = par("usr")[2] * 0.95 + par("usr")[1] * 0.05, y = par("usr")[3] * 0.98 + par("usr")[4] * 0.02)

}

dev.off()




# now look at it by tissue
tissue <- "t58-heart"
par(xpd = NA, mfrow = c(1,1))
plot(100,100,xlim = c(1,4), ylim = c(0,1), xpd=NA, 
     ylab = "Proportion Positive Effects on GWAS Scale", xlab = "Timepoint", xaxt = "n", bty="n",
)
segments(x0 = 1:4, x1 = 1:4, y0 = 0, y1 = - 0.02)
segments(x0 = 1, x1 = 4, y0 = 0, y1 = 0)
text(x = 1:4, y = -0.02, labels = paste0(2^(0:3), "w"), pos = 1)
for(i in 1:nrow(proportion_positive_effects[tissue,,,"prop",sex])){
  
  if(all(is.na(proportion_positive_effects[tissue,i,,"prop",sex]))){
    next()
  }
  
  lines(1:4, proportion_positive_effects[tissue,i,,"prop",sex], 
        lwd = 2, col = adjustcolor(cols$Tissue[tissue], 0.5))
  
}
segments(x0 = 0.9, y0 = 0.5, x1 = 4.1, y1 = 0.5, lty = 3, lwd = 3, col = "lightgrey")


#### graph of sign(al) propagation workflow ####

library(DiagrammeR)
library(DiagrammeRsvg)  # for conversion to svg
library(magrittr)
library(rsvg)  # for saving svg
std_equation <- latex2exp::TeX('\\frac{\\frac{z}{|z|}(\\frac{tanh(Q_{N}(1-\\frac{p}{2}))}{\\sqrt{n-3}})}{\\sqrt{2f(1-f)}}')[[1]]
graph_spec <- "digraph {graph [layout = dot, rankdir = TB, compound=true, splines=ortho]
newrank='true'
# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = White]

DEG [label=<<B><FONT COLOR='#8B0000'>ASDF-123</FONT></B><br/><B>DE</B>  log₂FC = 3>, height = 0.65, shape = egg, penwidth = 2]

eQTL [shape = record, label = <lead <B>cis-eQTL</B><br/>log₂FC = 2 | 
  ref <B><FONT COLOR='red'>A</FONT></B><br/>alt <B><FONT COLOR='blue'>G</FONT></B>>,
  height = 0.65]

GWAS [shape = record, label = <<B>GWAS</B>  Effect Size<br/><i>p-value</i>  = 0.017<br/><i>z-score</i>  = -2.12 | 
  ref <B><FONT COLOR='blue'>G</FONT></B><br/>alt <B><FONT COLOR='red'>A</FONT></B>>,
  height = 0.65]

GWAS_flipped [shape = record, label = <<B>GWAS</B>  Effect Size<br/><i>p-value</i>  = 0.017<br/><i>z-score</i>  = +2.12 | 
  ref <B><FONT COLOR='red'>A</FONT></B><br/>alt <B><FONT COLOR='blue'>G</FONT></B>>,
  height = 0.65]


GWAS_info [label = <<B>GWAS</B><br/><i>n</i>  = 20,000<br/><i>freq</i><SUB><B><FONT COLOR='blue'>G</FONT></B></SUB>  = 0.13>,
  height = 0.85]

gene_equiv [label=<<B>DE</B>  corresponds to<br/>effect of 1.5 <B><FONT COLOR='blue'>G</FONT></B>  alleles>]

GWAS_std [label=<<B>GWAS</B>  effect per<B><FONT COLOR='blue'>G</FONT></B><br/> allele in units SD<br/> = 0.0355>]

exp_effect [label=<Expected effect<br/>of <B>DE</B>  at<B><FONT COLOR='#8B0000'>ASDF-123</FONT></B>:<br/> 1.5 x 0.0355<br/> = 0.0532>]

summation [label = <sum over all <B>DE <FONT COLOR='#8B0000'>genes</FONT></B><br/>for particular <B>GWAS</B>>]

subgraph cluster_singlegene {
  label = <<B>Single Gene</B>>
  fontsize = 20
  graph[style=dashed]
  width = 200
  
  DEG -> eQTL [label = ' match to cis-eQTL', decorate=true]
  eQTL -> GWAS [label = ' Same Locus (rs12345)', decorate=true]
  DEGeQTL [shape=point,width=0.01,height=0.01];
  DEG -> DEGeQTL [dir=none]
  eQTL -> DEGeQTL [dir=none] 
  DEGeQTL -> gene_equiv [label = ' 3 ÷ 2', decorate = true]
  
  flip_alleles [label = <flip sign to<br/>match alleles>, shape = none]
  GWAS -> flip_alleles [dir=none]
  flip_alleles -> GWAS_flipped 
  GWAS_flipped -> GWAS_info [dir=none]
  GWAS_comb [shape=point,width=0.01,height=0.01];
  GWAS_flipped -> GWAS_comb [dir=none]
  GWAS_info -> GWAS_comb [dir=none]
  GWAS_comb -> GWAS_std
  
  all_comb [shape=point,width=0.01,height=0.01];
  GWAS_std -> all_comb [dir=none]
  gene_equiv -> all_comb [dir=none]
  all_comb -> exp_effect 

  {rank = same; DEG; eQTL;}    
  {rank = same; GWAS_flipped; GWAS_info;}    
  {rank = same; gene_equiv; GWAS_std;}    

}

exp_effect -> summation


}
"

library(magick)


grViz(graph_spec)
export_graph(graph_spec,
             file_name = "pic.png",
             file_type = "png")
grViz(graph_spec) %>% export_svg %>% charToRaw %>% rsvg_png("~/Documents/DExEQTL/composite_signal_workflow.png")
# image_read_pdf("~/Documents/DExEQTL/composite_signal_workflow.pdf", density = 80)


