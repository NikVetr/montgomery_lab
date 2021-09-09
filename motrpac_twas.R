#### TWAS ####

library(foreach)
library(doParallel)
library(parallel)

#initialize parallelization
if(!exists("cl")){
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

library(data.table)

#get sumstats & twas directories
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
twas_directory <- "~/repos/fusion_twas-master/"
command_changedir <- paste0("source ~/.bash_profile; cd ", twas_directory, "; ")

#undo sparsity in the the sumstats
for(gi in 1:length(gwas_summary_files)){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  sumstats <- fread(paste0("~/repos/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
  sumstats <- sumstats[!(sumstats$A1 == "" | is.na(sumstats$Z)),]
  fwrite(x = sumstats, sep = "\t",
         paste0(twas_directory, "sumstats/", gsub(x = gwas_summary_files[gi], pattern = ".txt.gz", replacement = ""), ".sumstats.gz"))
  fwrite(x = sumstats, sep = "\t",
         paste0(twas_directory, "sumstats/", gsub(x = gwas_summary_files[gi], pattern = ".txt.gz", replacement = ""), ".sumstats"))
}

#unzip (untar?) all tissue-specific expression data
tissue_tars <- list.files("~/repos/fusion_twas-master/WEIGHTS/")
tissue_tars <- tissue_tars[grep(tissue_tars, pattern = "tar")]
tissue_tars <- tissue_tars[-grep(tissue_tars, pattern = "ALL.tar")]
for(tt in tissue_tars){print(tt); system(paste0(command_changedir, "cd WEIGHTS; tar xjf ", tt))}

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

sapply(motrpac_gtex_map, function(tt) tissue_tars[grep(tt, tissue_tars)])

tissue_position_files <- list.files("~/repos/fusion_twas-master/WEIGHTS/")
tissue_position_files <- tissue_position_files[grep("[.]pos", tissue_position_files)]
motrpac_gtex_map <- motrpac_gtex_map[sapply(motrpac_gtex_map, function(tt) length(tissue_position_files[grep(tt, tissue_position_files)])) == 1]

#now let's try to do TWAS
foreach(gi=1:length(gwas_summary_files)) %dopar% {
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  # for(gi in 1:1){
  for(cri in 1:22){
    for(tissue in motrpac_gtex_map){
      if(!dir.exists(paste0("~/repos/fusion_twas-master/output/", tissue))){dir.create(paste0("~/repos/fusion_twas-master/output/", tissue))}
      command_twas <- paste0("Rscript FUSION.assoc_test.R --sumstats ", twas_directory, "sumstats/", 
                             gsub(x = gwas_summary_files[gi], pattern = ".txt.gz", replacement = ""), 
                             ".sumstats.gz --weights ./WEIGHTS/", tissue_position_files[grep(tissue, tissue_position_files)], " --weights_dir ./WEIGHTS/", 
                             " --ref_ld_chr ./LDREF/1000G.EUR. --chr ", cri, 
                             " --out output/", tissue, "/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "-chr_", cri, ".dat")
      command <- paste0(command_changedir, command_twas)
      system(command)      
    }
  }
}                         

#read in all of the results 
if(!file.exists("~/repos/fusion_twas-master/output/all_results.txt")){
  twas_tissues <- data.table()
  for(gi in 1:length(gwas_summary_files)){
    cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
    for(cri in 1:22){
      for(tissue in motrpac_gtex_map){
        temp <- fread(paste0(twas_directory, "output/", tissue, "/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "-chr_", cri, ".dat"))
        temp$tissue <- tissue
        temp$trait <- strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1]
        twas_tissues <- rbind(twas_tissues, temp)
      }
    }
  }
  fwrite(x = twas_tissues, file = "~/repos/fusion_twas-master/output/all_results.txt")
} else {
  twas_tissues <- fread(file = "~/repos/fusion_twas-master/output/all_results.txt")
}
hist(log(twas_tissues$TWAS.P), breaks = -(0:800))
mean(log(twas_tissues$TWAS.P) < (log(0.05) - log(nrow(twas_tissues))), na.rm = T)  * 100
sort(sapply(unique(twas_tissues$trait), function(trt) mean(log(twas_tissues$TWAS.P[twas_tissues$trait == trt]) < (log(0.05) - log(nrow(twas_tissues))), na.rm = T)) * 100, T)

#could also do metaxscan
#first gotta separate out loci on the basis of chromosome
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
metaxscan_dir <- "~/repos/MetaXcan/"
for(gi in 1:length(gwas_summary_files)){
  
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  
  if(!dir.exists(paste0(metaxscan_dir, "sumstats/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1]))){
    dir.create(paste0(metaxscan_dir, "sumstats/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1]))}
  
  sumstats <- fread(paste0(gwas_dir, c(gwas_summary_files[gi])))
  sumstats_sub <- sumstats[,c("chromosome","variant_id", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue")]
  colnames(sumstats_sub) <- c("CHR", "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P")
  sumstats_sub <- lapply(1:22, function(cri) sumstats_sub[sumstats_sub$CHR == paste0("chr", cri),])
  
  for(cri in 1:22){
    sumstats_sub_cri <- sumstats_sub[[cri]][,-c("CHR")]
    sumstats_sub_cri <- sumstats_sub_cri[!is.na(sumstats_sub_cri$SNP),]
    sumstats_sub_cri <- as.data.frame(sumstats_sub_cri)
    for(ci in 1:ncol(sumstats_sub_cri)){
      sumstats_sub_cri[is.na(sumstats_sub_cri[,ci]),ci] <- "NA"
    }
    
    fwrite(x = sumstats_sub_cri, sep = "\t",
           paste0(metaxscan_dir, "software/sumstats/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "/chr", cri, ".assoc.dosage.gz"))
    
  }
  
}

#run basic metaxscan
command_changedir <- paste0("source ~/.bash_profile; cd ", metaxscan_dir, "/software; ")
foreach(gi=1:length(gwas_summary_files)) %dopar% {
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  command_metaxscan <- paste0("python3 SPrediXcan.py ",
                              "--model_db_path data/DGN-WB_0.5.db ", 
                              "--covariance data/covariance.DGN-WB_0.5.txt.gz ", 
                              "--gwas_folder sumstats/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], " ",
                              "--gwas_file_pattern \".*gz\" ", 
                              "--snp_column SNP ", 
                              "--effect_allele_column A1 ", 
                              "--non_effect_allele_column A2 ", 
                              "--beta_column BETA ", 
                              "--pvalue_column P ", 
                              "--output_file results/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1],".csv")
  command <- paste0(command_changedir, command_metaxscan)
  system(command)      
}


#read in all of the results 
if(!file.exists("~/repos/MetaXcan/software/results/all_results.txt")){
  metaxscan_results <- data.table()
  for(gi in 1:length(gwas_summary_files)){
    
    cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
    
    temp <- fread(paste0(metaxscan_dir, "software/results/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], ".csv"))
    temp$trait <- strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1]
    metaxscan_results <- rbind(metaxscan_results, temp)
      
  }
  fwrite(x = metaxscan_results, file = "~/repos/MetaXcan/software/results/all_results.txt")
} else {
  metaxscan_results <- fread(file = "~/repos/MetaXcan/software/results/all_results.txt")
}
hist(log(metaxscan_results$pvalue), breaks = -(0:800))
mean(log(metaxscan_results$pvalue) < (log(0.05) - log(nrow(metaxscan_results))), na.rm = T)  * 100
sort(sapply(unique(metaxscan_results$trait), function(trt) mean(log(metaxscan_results$pvalue[metaxscan_results$trait == trt]) < (log(0.05) - log(nrow(metaxscan_results))), na.rm = T)) * 100, T)
