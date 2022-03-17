library(data.table)

paste0("python2 ./run_mesc.py ",
       "--compute-expscore-indiv ",
       "--plink-path /Users/nikgvetr/repos/plink/plink ",
       "--expression-matrix /Users/nikgvetr/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_expression_matrices/Muscle_Skeletal.v8.normalized_expression.bed.gz  ",
       "--columns 4,1,2,5 ",
       "--exp-bfile /Users/nikgvetr/repos/gtex-pipeline/expression_genotypes/GTEx_v8 ",
       "--geno-bfile /Users/nikgvetr/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.1 ",
       "--chr 1 ",
       "--covariates /Users/nikgvetr/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.transpose.txt ",
       "--out testrun"
)

  

#get gtex files into appropriate format
#GTEx sumstats files?
# d <- fread("/Volumes/SSD500GB/GTEx_Analysis_v8_eQTL_all_associations/Adipose_Subcutaneous.allpairs.txt")
# str(d)

#get covariates files in the proper format -- do for all tissues
covf <- fread("~/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.txt")
covf <- transpose(covf, keep.names = "ID")
colnames(covf) <- as.character(covf[1,])
covf <- as.data.frame(covf[-1,])
rownames(covf) <- as.character(covf[,1])
covf <- covf[,-1]
famf <- fread("~/repos/gtex-pipeline/expression_genotypes/GTEx_v8_noChr.fam")
covf <- cbind(famf[match(rownames(covf), unlist(famf[,2])),1:2], covf)
colnames(covf)[1:2] <- c("FID", "IID")
fwrite(covf, "~/repos/mesc/GTEx_covariates/Muscle_Skeletal.v8.covariates.txt", sep = "\t", row.names = F, col.names = T)

#get genesets for all tissues
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

#

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