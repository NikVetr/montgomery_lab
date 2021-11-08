#### revTWMR ####
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

#specify basic commands
revTWMR_directory <- "~/repos/revTWMR/"
gtex_pipeline_directory <- "~/repos/gtex-pipeline/"
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
command_changedir <- paste0("source ~/.bash_profile; cd ", revTWMR_directory, "; ")

#take a look at what these are doing
example_betaGWAS <- fread(paste0(revTWMR_directory, "bmi.matrix.betaGWAS"))
example_genes.N <- fread(paste0(revTWMR_directory, "genes.N"), header = T)
dim(example_betaGWAS)
"R < revTWMR.R --no-save bmi.matrix.betaGWAS bmi genes.N"

#retrieve GTEx eQTLs
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

GTEx_eQTLs = c("Adipose_Subcutaneous.allpairs.txt.gz", "
                Adipose_Visceral_Omentum.allpairs.txt.gz", "
                Adrenal_Gland.allpairs.txt.gz", "
                Artery_Aorta.allpairs.txt.gz", "
                Artery_Coronary.allpairs.txt.gz", "
                Artery_Tibial.allpairs.txt.gz", "
                Brain_Amygdala.allpairs.txt.gz", "
                Brain_Anterior_cingulate_cortex_BA24.allpairs.txt.gz", "
                Brain_Caudate_basal_ganglia.allpairs.txt.gz", "
                Brain_Cerebellar_Hemisphere.allpairs.txt.gz", "
                Brain_Cerebellum.allpairs.txt.gz", "
                Brain_Cortex.allpairs.txt.gz", "
                Brain_Frontal_Cortex_BA9.allpairs.txt.gz", "
                Brain_Hippocampus.allpairs.txt.gz", "
                Brain_Hypothalamus.allpairs.txt.gz", "
                Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt.gz", "
                Brain_Putamen_basal_ganglia.allpairs.txt.gz", "
                Brain_Spinal_cord_cervical_c-1.allpairs.txt.gz", "
                Brain_Substantia_nigra.allpairs.txt.gz", "
                Breast_Mammary_Tissue.allpairs.txt.gz", "
                Cells_Cultured_fibroblasts.allpairs.txt.gz", "
                Cells_EBV-transformed_lymphocytes.allpairs.txt.gz", "
                Colon_Sigmoid.allpairs.txt.gz", "
                Colon_Transverse.allpairs.txt.gz", "
                Esophagus_Gastroesophageal_Junction.allpairs.txt.gz", "
                Esophagus_Mucosa.allpairs.txt.gz", "
                Esophagus_Muscularis.allpairs.txt.gz", "
                Heart_Atrial_Appendage.allpairs.txt.gz", "
                Heart_Left_Ventricle.allpairs.txt.gz", "
                Kidney_Cortex.allpairs.txt.gz", "
                Liver.allpairs.txt.gz", "
                Lung.allpairs.txt.gz", "
                Minor_Salivary_Gland.allpairs.txt.gz", "
                Muscle_Skeletal.allpairs.txt.gz", "
                Nerve_Tibial.allpairs.txt.gz", "
                Ovary.allpairs.txt.gz", "
                Pancreas.allpairs.txt.gz", "
                Pituitary.allpairs.txt.gz", "
                Prostate.allpairs.txt.gz", "
                Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz", "
                Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz", "
                Small_Intestine_Terminal_Ileum.allpairs.txt.gz", "
                Spleen.allpairs.txt.gz", "
                Stomach.allpairs.txt.gz", "
                Testis.allpairs.txt.gz", "
                Thyroid.allpairs.txt.gz", "
                Uterus.allpairs.txt.gz", "
                Vagina.allpairs.txt.gz", "
                Whole_Blood.allpairs.txt.gz")

GTEx_eQTLs <- GTEx_eQTLs[sapply(motrpac_gtex_map, function(x) grep(x, GTEx_eQTLs))]
GTEx_eQTLs <- gsub(x = GTEx_eQTLs, pattern = "\\n                ", replacement = "")
for(tissue in GTEx_eQTLs){
  print(tissue)
  if(!file.exists(paste0("~/data/smontgom/GTEx_Analysis_v8_eQTL_all_associations/", tissue))){
      system(paste0("scp -r nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/", tissue, " ~/data/smontgom/GTEx_Analysis_v8_eQTL_all_associations"))  
  }
}


if(!dir.exists(paste0(revTWMR_directory, "betaGWAS"))){dir.create(paste0(revTWMR_directory, "betaGWAS"))}
if(!dir.exists(paste0(revTWMR_directory, "genesN"))){dir.create(paste0(revTWMR_directory, "genesN"))}

#parse SNP IDs
RSID_POS_MAP <- fread("~/data/smontgom/RSID_POS_MAP.txt")
colnames(RSID_POS_MAP)[1] <- "CHROM"
RSID_POS_MAP <- lapply(setNames(unique(RSID_POS_MAP$CHROM), unique(RSID_POS_MAP$CHROM)), function(cri) RSID_POS_MAP[CHROM == cri])

#slice up the positional map for easier processing
for(cri in names(RSID_POS_MAP)){
  print(cri)
  fwrite(RSID_POS_MAP[[cri]], paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
}

#slice up the eQTL files by chromosome for easier processsing
for(tissue in names(motrpac_gtex_map)){
  print(tissue)
  eQTL_sumstats <- fread(paste0("~/data/smontgom/GTEx_Analysis_v8_eQTL_all_associations/", 
                                GTEx_eQTLs[grep(motrpac_gtex_map[match(tissue, names(motrpac_gtex_map))], GTEx_eQTLs)])[1],
                         select = c("gene_id", "variant_id", "slope"))
  CHROMs <- substr(eQTL_sumstats$variant_id, 4, 5)
  CHROMs <- gsub("_", "", CHROMs)
  CHROMs <- lapply(setNames(c(1:22, "X"), c(1:22, "X")), function(cri) which(CHROMs == cri))
  for(cri in c(1:22, "X")){
    fwrite(eQTL_sumstats[CHROMs[[cri]]], 
           file = paste0("~/data/smontgom/GTEx_Analysis_v8_eQTL_all_associations/", tissue, "_", cri, ".txt.gz"))
  }
  
  rm(eQTL_sumstats)
  gc()
}

#process all the input files needed for revTWMR
for(tissue in names(motrpac_gtex_map)){
  
  print(tissue)
  
  for(cri in c(1:22, "X")){
   
    cat(paste0(" ", cri))
    
    RSID_POS_MAP <- fread(paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
    eQTL_sumstats <- fread(paste0("~/data/smontgom/GTEx_Analysis_v8_eQTL_all_associations/", tissue, "_", cri, ".txt.gz"))
    vids <- eQTL_sumstats$variant_id
    vids <- substr(vids, start = 4, nchar(vids))
    vids <- strsplit(vids, "_", T)
    vids <- as.data.table(data.frame(data.table::transpose(vids)))
    colnames(vids) <- c("CHROM", "POS", "REF", "ALT", "BUILD")
    vids$RSID <- RSID_POS_MAP$ID[match(vids$POS, RSID_POS_MAP$POS)]
    eQTL_sumstats <- cbind(eQTL_sumstats, vids)
    eQTL_sumstats <- eQTL_sumstats[!is.na(eQTL_sumstats$RSID)]
    eQTL_sumstats$ENSG <- gsub('\\..*','',eQTL_sumstats$gene_id)
    
    for(gi in 1:length(gwas_summary_files)){
      cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
      gwas_sumstats <- fread(paste0("~/repos/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
      gwas_sumstats <- gwas_sumstats[!(gwas_sumstats$A1 == "" | is.na(gwas_sumstats$Z)),]
      #A1 = REF
      
      
      
    }
    
    #clean up memory
    rm(vids)
    rm(RSID_POS_MAP)
    rm(eQTL_sumstats)
    gc()
    
  }
  
}

for(gi in 1:length(gwas_summary_files)){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  gwas_sumstats <- fread(paste0("~/repos/ldsc/gwas_sumstats/proper_format/signed/", gwas_summary_files[gi], ".sumstats.gz"))
  gwas_sumstats <- gwas_sumstats[!(gwas_sumstats$A1 == "" | is.na(gwas_sumstats$Z)),]
  #A1 = REF
}
