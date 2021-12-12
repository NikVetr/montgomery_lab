#load packages
library(arrow)
library(data.table)
library(cmdstanr)
library(posterior)

#functions
removeNA <- function(x) x[!is.na(x)]

## essential metadata

#eqtls
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
GTEx_SampleSize <- read.csv("~/data/smontgom/GTEx_SampleSizes.csv", header = T)
gene_files <- list.files(paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri, "/"))
gene_files <- gene_files[grep("full", gene_files)]

#general
RSID_POS_MAP <- fread(paste0("~/data/smontgom/RSID_POS_MAP_",cri,".txt"))
ENSG_coord <- read.table(file = "~/repos/ldsc/ENSG_coord.txt", header = T)

#gwas
ldsc_directory <- "~/repos/ldsc/"
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[grep(gwas_summary_files, pattern = ".txt.gz")]
traitnames <- gwas_summary_files
traitnames <- sapply(traitnames, function(trait) strsplit(trait, ".txt.gz")[[1]][1])
traitnames <- sapply(traitnames, function(trait) strsplit(trait, "imputed_")[[1]][2])
traitnames <- as.character(traitnames)

#reprocess eqtl sumstats
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
      eQTL_sumstats$ID <- apply(vids[,c("CHROM", "POS", "REF", "ALT")], 1, paste0, collapse = "_")
      eQTL_sumstats <- eQTL_sumstats[!is.na(eQTL_sumstats$RSID),]
      eQTL_sumstats$Z <- round(eQTL_sumstats$slope / eQTL_sumstats$slope_se, 5)
      eQTL_sumstats$N <-  paste0(GTEx_SampleSize$sample_size[GTEx_SampleSize$tissue == tissue], ".00")
      
      genes <- unique(eQTL_sumstats$ENSG)
      eQTL_sumstats <- split(eQTL_sumstats , f = eQTL_sumstats$ENSG)
      eQTL_sumstats <- lapply(eQTL_sumstats, function(gid){
        gene_info <- gid[,c("RSID", "REF", "ALT", "Z", "N", "af", "slope", "slope_se", "POS", "ID")] #A2 / ALT is the effect allele
        colnames(gene_info) <- c("SNP", "A1", "A2", "Z", "N", "ALLELE_FREQ", "SLOPE", "SLOPE_SE", "POS", "ID")
        return(gene_info)
      })
      
      for(gene_name in names(eQTL_sumstats)){
        
        fwrite(eQTL_sumstats[[gene_name]], file = paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri, "/", gene_name, "_full.sumstats.gz"), sep = "\t", append = F, col.names = T)
      }
      
    }
  }
}

#specify targets
tissue = tissues[1]
cri = 1
gwas_trait <- gwas_summary_files[1]
gwas_traitname <- gsub(".txt.gz", "", gsub("imputed_", "", gwas_trait))
gene_file <- gene_files[1]
gene_name <- gsub("_full.sumstats.gz", "", gene_file)
window_size <- 1E5

#load data
gwas_sumstats <- fread(paste0(gwas_dir, gwas_trait))
gwas_sumstats$chromosome <- as.integer(substr(gwas_sumstats$chromosome, 4, nchar(gwas_sumstats$chromosome)))
gwas_sumstats <- split(x = gwas_sumstats, by = "chromosome")
gene_sumstats <- fread(file = paste0(ldsc_directory, "/GTEx_v8_log2norm_sumstats/", tissue, "/", cri, "/", gene_file))

#match loci
loc <- ENSG_coord[ENSG_coord$GENE == strsplit(gene_file, split = "[_]")[[1]][1],c("CHR", "START", "END")]
RSIDs_in_window <- RSID_POS_MAP$ID[(RSID_POS_MAP$POS > (loc$START - window_size) & RSID_POS_MAP$POS < loc$START) | (RSID_POS_MAP$POS > (loc$END) & RSID_POS_MAP$POS < loc$END + window_size)]
in_gene_sumstats <- removeNA(gene_sumstats$SNP[match(RSIDs_in_window, gene_sumstats$SNP)])
in_gwas_sumstats <- removeNA(gwas_sumstats[[cri]]$variant_id[match(RSIDs_in_window, gwas_sumstats[[cri]]$variant_id)])
in_both <- intersect(in_gene_sumstats, in_gwas_sumstats)

#merge sumstats
gwas_sub <- gwas_sumstats[[cri]][match(in_both, gwas_sumstats[[cri]]$variant_id),
                                 c("variant_id", "effect_allele", "non_effect_allele", "frequency", "sample_size", "effect_size", "standard_error")]
gwas_sub <- gwas_sub[complete.cases(gwas_sub),]
gene_sub <- gene_sumstats[match(gwas_sub$variant_id, gene_sumstats$SNP),]

sign_flip <- rep(NA, nrow(gwas_sub))
sign_flip[gwas_sub$effect_allele == gene_sub$A2 & gwas_sub$non_effect_allele == gene_sub$A1] <- 1
sign_flip[gwas_sub$effect_allele == gene_sub$A1 & gwas_sub$non_effect_allele == gene_sub$A2] <- -1
gene_sub$SLOPE <- gene_sub$SLOPE * sign_flip
gene_sub <- gene_sub[!is.na(sign_flip),]
gwas_sub <- gwas_sub[!is.na(sign_flip),]

plot(gene_sub$SLOPE, gwas_sub$effect_size)
summary(lm(gene_sub$SLOPE~ gwas_sub$effect_size))



#get LD matrix
if(!dir.exists(paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname))){dir.create(paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname))}
if(!dir.exists(paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname, "/chr", cri, "/"))){dir.create(paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname, "/chr", cri, "/"))}
sink(file = paste0("/Volumes/SSD500GB/all.EUR.hg38a/", gwas_traitname,"/chr", cri, "/", gene_name,"_", as.integer(window_size),"BP-Window_SNPs.txt"))
cat(paste0(gene_sub$ID, "\n"), sep = "")
sink()
long <- fread("/Volumes/SSD500GB/all.EUR.hg38a/test_SNPs_LD.ld.gz")
pos <- as.character( sort(unique(c(long$BP_A, long$BP_B))) )

library(reshape2)
irregular <- acast(long, BP_A ~ BP_B, value.var="R")  # wide format but irregular shape

mat.r2 <- matrix(NA, length(pos), length(pos), dimnames=list(pos, pos))
mat.r2[ rownames(irregular), colnames(irregular) ] <- irregular
mat.r2[ lower.tri(mat.r2) ] <- t(mat.r2)[ lower.tri( t(mat.r2) ) ]
diag(mat.r2) <- 1
rm(tmp)
#index genome vcf file
# read_lines_gzip = function(infile, lines, outfile = NULL){
#   lines = as.integer(lines)
#   sed_arg = paste(lines, "p", sep = "", collapse = ";")
#   sed_arg = paste(sed_arg, ";", max(lines) + 1, "q", sep = "")
#   sed_command = sprintf("sed -n '%s'", sed_arg)
#   out_command = if(is.null(outfile)) "" else sprintf("> %s", as.character(outfile))
#   command = sprintf("gunzip -c %s | %s %s", infile, sed_command, out_command)
#   system(command, intern = is.null(outfile))
# }
# 
# gtex_dir <- "~/repos/gtex-pipeline/"
# infile = paste0(gtex_dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz")

sapply(1:10, function(x) )
lines = 1

cat(read_lines_gzip(infile, 30000))

#fit stan model
file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)

# names correspond to the data block in the Stan program
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 5E3,
  iter_sampling = 5E3
)

fit$summary()

fit_vb <- mod$variational(data = data_list, seed = 123, output_samples = 4000) 
