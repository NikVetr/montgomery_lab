library(foreach)
library(doParallel)
library(parallel)
library(MotrpacBicQC)
library(data.table)

# paste0("python2 ./run_mesc.py ",
#        "--compute-expscore-indiv ",
#        "--plink-path /Users/nikgvetr/repos/plink/plink ",
#        "--expression-matrix /Users/nikgvetr/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz  ",
#        "--columns 4,1,2,5 ",
#        "--exp-bfile /Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8 ",
#        "--geno-bfile /Users/nikgvetr/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.1 ",
#        "--chr 1 ",
#        "--covariates /Users/nikgvetr/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.transpose.txt ",
#        "--out testrun"
# )
# 
#   
# 
# #get gtex files into appropriate format
# #GTEx sumstats files?
# # d <- fread("/Volumes/SSD500GB/GTEx_Analysis_v8_eQTL_all_associations/Adipose_Subcutaneous.allpairs.txt")
# # str(d)
# 
# #get covariates files in the proper format -- do for all tissues
# covf <- fread("~/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.txt")
# covf <- transpose(covf, keep.names = "ID")
# colnames(covf) <- as.character(covf[1,])
# covf <- as.data.frame(covf[-1,])
# rownames(covf) <- as.character(covf[,1])
# covf <- covf[,-1]
# famf <- fread("~/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr.fam")
# covf <- cbind(famf[match(rownames(covf), unlist(famf[,2])),1:2], covf)
# colnames(covf)[1:2] <- c("FID", "IID")
# fwrite(covf, "~/repos/mesc/GTEx_covariates/Muscle_Skeletal.v8.covariates.txt", sep = "\t", row.names = F, col.names = T)


#get genesets for all tissues
load("~/data/smontgom/node_metadata_list.RData") #from "~/scripts/montgomery_lab/motrpac-DEG_barbeira-twas_intersection.R"
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

tissues <- names(motrpac_gtex_map)
mesc_compatible_genesets <- lapply(setNames(tissues, tissues), function(tiss){
  if(!file.exists(paste0("~/data/smontgom/GTEx_v8_ExpressionScores/tissues/", 
                         motrpac_gtex_map[tiss], "/", motrpac_gtex_map[tiss], ".", 1, ".gannot.gz"))){
    return(NULL)
  }
  genes <- node_metadata_list[["8w"]]$human_gene_symbol[node_metadata_list[["8w"]]$tissue == tissue_abbr[tiss]]
  genes <- genes[!is.na(genes)]
  mesc_genes <- unlist(lapply(1:22, function(cri) fread(paste0("~/data/smontgom/GTEx_v8_ExpressionScores/tissues/", 
                               motrpac_gtex_map[tiss], "/", motrpac_gtex_map[tiss], ".", cri, ".gannot.gz"))$Gene))
  print(paste0(tiss, ", ", round(mean(genes %in% mesc_genes), 3)))
  intersect(genes, mesc_genes)
})

geneset_dir <- "~/repos/mesc/genesets/"
for(tiss in names(mesc_compatible_genesets)){
  print(tiss)
  if(!dir.exists(geneset_dir)){dir.create(geneset_dir)}
  sink(file = paste0(geneset_dir, tiss, ".txt"))
  cat(paste0(tiss, "\t", paste(mesc_compatible_genesets[[tiss]], collapse = "\t")))
  sink()
}

#get geneset expression scores with mesc
output_dir <- "~/repos/mesc/output/Expression_Scores/"
if(!dir.exists(output_dir)){dir.create(output_dir)}
command_setup <- paste0("source ~/.bash_profile; cd ~/repos/mesc/; source activate mesc; ")

if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
}

tissues_to_analyze <- names(mesc_compatible_genesets)[sapply(mesc_compatible_genesets, function(x) length(x) > 0)]
foreach(tiss = tissues_to_analyze) %dopar% {
# for(tiss in tissues_to_analyze){
  
  tissue_dir <- paste0("~/repos/mesc/output/Expression_Scores/", tiss)
  if(!dir.exists(tissue_dir)){dir.create(tissue_dir)}
  
  for(cri in 1:22){
    cat(paste0(tiss, ": ", cri, "\n"))
    if(readLines("~/foreach_continue.txt", warn = F) != "T"){
      next()
    }
    
    command_mesc <- paste0("python2 ./gene_set_analysis.py ",
                           "--input-prefix /Users/nikgvetr/data/smontgom/GTEx_v8_ExpressionScores/tissues/", motrpac_gtex_map[tiss], "/", motrpac_gtex_map[tiss], " ",
                           "--gene-sets genesets/", tiss, ".txt ",
                           "--bfile /Users/nikgvetr/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, " ",
                           "--chr ", cri, " ",
                           "--out output/Expression_Scores/", tiss, "/", tiss, "; ")
    command <- paste0(command_setup, command_mesc)
    system(command)
    
  }
}

#compute heritability enrichment for all of the above expression scores
gwas_dir <- "~/repos/ldsc/gwas_sumstats/proper_format/"
gwas_files <- list.files(gwas_dir)
gwas_files <- gwas_files[grep("sumstats.gz", gwas_files)]
tissues_to_analyze <- names(mesc_compatible_genesets)[sapply(mesc_compatible_genesets, function(x) length(x) > 0)]
foreach(tiss = tissues_to_analyze) %dopar% {
  
  #create output directory
  tissue_dir <- paste0("~/repos/mesc/output/Enrichment/", tiss)
  if(!dir.exists(tissue_dir)){dir.create(tissue_dir)}
  
  for(gwas in gwas_files){
    
    #bookkeeping
    cat(paste0(tiss, ": ", gwas, "\n"))
    if(readLines("~/foreach_continue.txt", warn = F) != "T"){
      next()
    }
    
    #run mesc
    command_mesc <- paste0("python2 ./run_mesc.py ",
                           "--h2med /Users/nikgvetr/repos/ldsc/gwas_sumstats/proper_format/", gwas, " ",
                           "--exp-chr output/Expression_Scores/", tiss, "/", tiss, " ",
                           "--out output/Enrichment/", tiss, "/", gsub("imputed_", "", gsub(pattern = ".txt.gz.sumstats.gz", "", gwas)), "; ")
    command <- paste0(command_setup, command_mesc)
    system(command)
    
  }
}


#### read in results ####
mesc_output_categories <- do.call(rbind, lapply(tissues_to_analyze, function(tiss){
  print(tiss)
  tissue_dir <- paste0("~/repos/mesc/output/Enrichment/", tiss, "/")
  mesc_out <- do.call(rbind, lapply(gwas_files, function(gwas){
    enrich <- fread(paste0(tissue_dir, gsub("imputed_", "", gsub(pattern = ".txt.gz.sumstats.gz", "", gwas)), ".categories.h2med"))
    enrich$trait <- gwas
    enrich <- enrich[!grepl("h2cis_bin", enrich$Gene_category),]
    enrich
  }))
  mesc_out$tissue <- tiss
  mesc_out
}))

mesc_output_basic <- do.call(rbind, lapply(tissues_to_analyze, function(tiss){
  print(tiss)
  tissue_dir <- paste0("~/repos/mesc/output/Enrichment/", tiss, "/")
  mesc_out <- do.call(rbind, lapply(gwas_files, function(gwas){
    enrich <- fread(paste0(tissue_dir, gsub("imputed_", "", gsub(pattern = ".txt.gz.sumstats.gz", "", gwas)), ".all.h2med"))
    enrich$trait <- gwas
    enrich
  }))
  mesc_out$tissue <- tiss
  mesc_out
}))
hist(mesc_output_basic$Estimate[mesc_output_basic$Quantity == "h2med"])

#subset to sufficiently well-powered genesets
mesc_output_categories <- mesc_output_categories[mesc_output_categories$Num_genes > 200,]
hist(mesc_output_categories$h2med_enrichment, breaks = 1000)

hist(mesc_output_categories$h2med_enrichment_pvalue)
mesc_output_categories$h2med_enrichment_pvalue_BHadjust <- p.adjust(mesc_output_categories$h2med_enrichment_pvalue, method = "BH")
hist(mesc_output_categories$h2med_enrichment_pvalue_BHadjust, breaks = 20)
ihw_data <- data.frame(trait = as.factor(mesc_output_categories$trait), 
                                tissue = as.factor(mesc_output_categories$tissue), 
                                pvalue = mesc_output_categories$h2med_enrichment_pvalue)
mesc_output_categories$enrich_pval_ihw_trait <- IHW::ihw(pvalue ~ trait, data = ihw_data, alpha = 0.05)@df$adj_pvalue
mesc_output_categories$enrich_pval_ihw_tissue <- IHW::ihw(pvalue ~ tissue, data = ihw_data, alpha = 0.05)@df$adj_pvalue
pairs(cbind(raw_pval = mesc_output_categories$h2med_enrichment_pvalue,
            ihw_trait = mesc_output_categories$enrich_pval_ihw_trait, 
            ihw_tissue = mesc_output_categories$enrich_pval_ihw_tissue,
            bh = mesc_output_categories$h2med_enrichment_pvalue_BHadjust))
sum(mesc_output_categories$enrich_pval_ihw_trait < 0.05)
sum(mesc_output_categories$enrich_pval_ihw_tissue < 0.05)
sum(mesc_output_categories$h2med_enrichment_pvalue_BHadjust < 0.05)

mesc_output_categories[apply(cbind(ihw_trait = mesc_output_categories$enrich_pval_ihw_trait, 
      ihw_tissue = mesc_output_categories$enrich_pval_ihw_tissue,
      bh = mesc_output_categories$h2med_enrichment_pvalue_BHadjust) < 0.05, 1, sum) > 0,]


#### computing expressions cores from scratch ####

# bim_file <- fread("/Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr.bim")
# fwrite(bim_file, "/Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr_noRSID.bim", col.names = F)
bim_file <- fread("/Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr_noRSID.bim")
bim_file_rsids <- do.call(rbind, lapply(unique(bim_file$V1), function(cri){
  print(cri)
  sub <- bim_file[bim_file$V1 == cri,]
  vids <- sub$V2
  vids <- substr(vids, start = 4, nchar(vids))
  vids <- strsplit(vids, "_", T)
  vids <- as.data.table(data.frame(data.table::transpose(vids)))
  colnames(vids) <- c("CHROM", "POS", "REF", "ALT", "BUILD")
  
  RSID_POS_MAP <- fread(paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
  vids$RSID <- RSID_POS_MAP$ID[match(vids$POS, RSID_POS_MAP$POS)]
  sub$V2 <- vids$RSID
  print(1-mean(is.na(sub$V2)))
  sub
}))
1-mean(is.na(bim_file_rsids$V2))
bim_file_rsids$V2[is.na(bim_file_rsids$V2)] <- "NA"
bim_file_rsids$V2[bim_file_rsids$V2 == "rs123456789"] <- "NA"
bim_file_rsids[bim_file_rsids$V2 != "NA",]
fwrite(bim_file_rsids, "/Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr.bim", col.names = F, sep = "\t")

bim_file_1000G <- fread("/Users/nikgvetr/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bim")
head(bim_file_1000G)
head(bim_file_rsids)
mean(bim_file_1000G$V2 %in% bim_file_rsids$V2)
mean(bim_file_rsids$V2[bim_file_rsids$V1 == "1"] %in% bim_file_1000G$V2)

#TODO: change all 'chr1's to '1's, 'chr2's to '2's, etc.
#TODO: transpose covariates matrix