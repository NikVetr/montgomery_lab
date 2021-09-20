library(data.table)
gtex_pipeline_directory <- "~/repos/gtex-pipeline/"

sample_ids <- readLines("~/repos/gtex-pipeline/sample_ids.txt", )
sample_ids <- strsplit(sample_ids, "\t")[[1]]
sample_ids <- sample_ids[grep("GTEX", sample_ids)]
participant_ids <- readLines("~/repos/gtex-pipeline/participant_ids.txt")

participant_to_samples <- sapply(participant_ids, function(pid) grep(pid, sample_ids))

length(sample_ids)
length(unique(unlist(participant_to_samples)))

unaccounted_for_participant_IDs <- unique(sapply(sample_ids[setdiff(1:length(sample_ids), as.integer(unlist(participant_to_samples)))], 
                                                 function(x) paste0(strsplit(x, "-")[[1]][1:2], collapse = "-")))

unaccounted_for_participant_IDs
intersect(unaccounted_for_participant_IDs, participant_ids)

#make file
participant_to_samples <- lapply(participant_to_samples, function(pid) sample_ids[pid])
sample_participant_lookup <- do.call(rbind, lapply(seq_along(participant_to_samples), function(i) cbind(sample_id = participant_to_samples[[i]], participant_id = names(participant_to_samples)[i])))
data.table::fwrite(x = sample_participant_lookup, file = "~/repos/gtex-pipeline/sample_participant_lookup.txt", sep = "\t")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## now let us do the E(count) gct files ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#trim out extra samples
gene_expected_count <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count.gct.gz")
cols_to_discard <- unlist(sapply(unaccounted_for_participant_IDs, function(pid) grep(pid, colnames(gene_expected_count))))
gene_expected_count_subset <- gene_expected_count[,-..cols_to_discard]
rm(gene_expected_count)

#write partial data to disk for troubleshooting
gene_expected_count_subset_subset <- gene_expected_count_subset[1:10,]
writeLines(paste0("#1.2\n", nrow(gene_expected_count_subset_subset), "\t", ncol(gene_expected_count_subset_subset)-2),
           con = "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset_subset.gct")
fwrite(gene_expected_count_subset_subset[,-"transcript_id(s)"], "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset_subset.gct", sep = "\t", 
       append = T, col.names = T)

#write full data to disc
# fwrite(gene_expected_count_subset, "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset.gct.gz", sep = "\t")
writeLines(paste0("#1.2\n", nrow(gene_expected_count_subset), "\t", ncol(gene_expected_count_subset)-2),
           con = "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset.gct")
fwrite(gene_expected_count_subset[,-"transcript_id(s)"], "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset.gct", sep = "\t", 
       append = T, col.names = T)
rm(gene_expected_count_subset)
system("cd ~/repos/gtex-pipeline/rna_seq/; gzip -f GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset.gct")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## now let us do the tpm gct files ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gene_tpm <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz")
cols_to_discard <- unlist(sapply(unaccounted_for_participant_IDs, function(pid) grep(pid, colnames(gene_tpm))))
gene_tpm_subset <- gene_tpm[,-..cols_to_discard]
rm(gene_tpm)

#troubleshooting
gene_tpm_subset_subset <- gene_tpm_subset[1:10,]
writeLines(paste0("#1.2\n", nrow(gene_tpm_subset_subset), "\t", ncol(gene_tpm_subset_subset)-2),
           con = "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset_subset.gct")
fwrite(gene_tpm_subset_subset[,-"transcript_id(s)"], "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset_subset.gct", sep = "\t",
       append = T, col.names = T, nThread = 8)

#full data
# fwrite(gene_tpm_subset, "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset.gct.gz", sep = "\t")
writeLines(paste0("#1.2\n", nrow(gene_tpm_subset), "\t", ncol(gene_tpm_subset)-2),
           con = "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset.gct")
fwrite(gene_tpm_subset[,-"transcript_id(s)"], "~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset.gct", sep = "\t",
       append = T, col.names = T, nThread = 8)
system("cd ~/repos/gtex-pipeline/rna_seq/; gzip -f GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset.gct")
rm(gene_tpm_subset)

# #sample ids must match?
# gene_tpm <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm_subset.gct.gz")
# gene_expected_count <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count_subset.gct.gz")
# lookup_table <- fread("~/repos/gtex-pipeline/sample_participant_lookup.txt")
# all.equal(colnames(gene_tpm)[-(1:2)], lookup_table$sample_id)
# all.equal(colnames(gene_expected_count)[-(1:2)], lookup_table$sample_id)

#ok, now let's partition each of these by tissue lol

sample_lookup_table <- fread("~/data/smontgom/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
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

#first let's do expected counts
gene_expected_count <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count.gct.gz")
sample_ids <- colnames(gene_expected_count)
sample_tissue_match <- sapply(sample_ids, function(sample) sample_lookup_table$SMTSD[match(sample, sample_lookup_table$SAMPID)])
sample_tissue_match <- names(motrpac_gtex_map)[match(sample_tissue_match, motrpac_gtex_map)]
sample_tissue_match[is.na(sample_tissue_match)] <- "irrelevant"
cols_to_discard <- unlist(sapply(unaccounted_for_participant_IDs, function(pid) grep(pid, colnames(gene_expected_count))))
if(!dir.exists(paste0(gtex_pipeline_directory, "expression_data"))){dir.create(paste0(gtex_pipeline_directory, "expression_data"))}

for(tissue in names(motrpac_gtex_map)){
        
        print(tissue)
        
        #subset the data
        samps_to_use <- which(sample_tissue_match == tissue)
        samps_to_use <- setdiff(samps_to_use, cols_to_discard)
        subdata <- gene_expected_count[,..samps_to_use]
        subdata <- cbind(gene_id = gene_expected_count$gene_id, subdata)
        
        #write to file
        writeLines(paste0("#1.2\n", nrow(subdata), "\t", ncol(subdata)-1),
                   con = paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count.gct"))
        fwrite(subdata, paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count.gct"), sep = "\t", 
               append = T, col.names = T)
        
        #gzip the file
        system(paste0("cd ",gtex_pipeline_directory, "expression_data/; gzip -f ", "GTEx_Analysis_v8_", tissue, "_expected_count.gct"))
        
}
rm(gene_expected_count)

#now let's do TPM
gene_tpm <- fread("~/repos/gtex-pipeline/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz")
sample_ids <- colnames(gene_tpm)
sample_tissue_match <- sapply(sample_ids, function(sample) sample_lookup_table$SMTSD[match(sample, sample_lookup_table$SAMPID)])
sample_tissue_match <- names(motrpac_gtex_map)[match(sample_tissue_match, motrpac_gtex_map)]
sample_tissue_match[is.na(sample_tissue_match)] <- "irrelevant"
cols_to_discard <- unlist(sapply(unaccounted_for_participant_IDs, function(pid) grep(pid, colnames(gene_tpm))))
if(!dir.exists(paste0(gtex_pipeline_directory, "expression_data"))){dir.create(paste0(gtex_pipeline_directory, "expression_data"))}

for(tissue in names(motrpac_gtex_map)){
        
        print(tissue)
        
        #subset the data
        samps_to_use <- which(sample_tissue_match == tissue)
        samps_to_use <- setdiff(samps_to_use, cols_to_discard)
        subdata <- gene_tpm[,..samps_to_use]
        subdata <- cbind(gene_id = gene_tpm$gene_id, subdata)
        
        #write to file
        writeLines(paste0("#1.2\n", nrow(subdata), "\t", ncol(subdata)-1),
                   con = paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm.gct"))
        fwrite(subdata, paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm.gct"), sep = "\t", 
               append = T, col.names = T)
        
        #gzip the file
        system(paste0("cd ",gtex_pipeline_directory, "expression_data/; gzip -f ", "GTEx_Analysis_v8_", tissue, "_tpm.gct"))
        
}
rm(gene_tpm)

#now to run these expression normalization scripts
library(data.table)
gtex_pipeline_directory <- "~/repos/gtex-pipeline/"
command_changedir <- paste0("source ~/.bash_profile; cd ", gtex_pipeline_directory, "; ")
if(!dir.exists(paste0(gtex_pipeline_directory, "log2-normalized-expression"))){dir.create(paste0(gtex_pipeline_directory, "log2-normalized-expression"))}
for(tissue in names(motrpac_gtex_map)){
        # foreach(tissue = names(motrpac_gtex_map)) %dopar% {
        
        print(tissue)
        
        command_gtex <- paste0("python3 qtl/src/eqtl_prepare_expression.py ", 
                               paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm.gct.gz "), 
                               paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count.gct.gz "), 
                               "gencode.v26.GRCh38.genes.gtf ", 
                               "sample_participant_lookup.txt ", 
                               "GTEx_Analysis_vcf_chr_list.txt ", 
                               paste0("log2-normalized-expression_", tissue, " "), 
                               "--tpm_threshold 0.1 ", 
                               "--count_threshold 6 ", 
                               "--sample_frac_threshold 0.2 ", 
                               "--normalization_method tmm; ")
        command_movefile_1 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".expression.bed.gz.tbi "), 
                                     paste0("log2-normalized-expression/log2-normalized-expression_", tissue, ".expression.bed.gz.tbi; "))
        command_movefile_2 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".expression.bed.gz "), 
                                     paste0("log2-normalized-expression/log2-normalized-expression_", tissue, ".expression.bed.gz;"))
        command <- paste0(command_changedir, command_gtex, command_movefile_1, command_movefile_2)
        system(command)      
}

#now process the next step of the pipeline, including covariates
motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map, pattern = "-", replacement = "_")
motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map_covariates, pattern = " ", replacement = "_")
motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map_covariates, pattern = "___", replacement = "_")
all(paste0(motrpac_gtex_map_covariates, ".v8.covariates.txt") %in% list.files("~/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates"))
if(!dir.exists(paste0(gtex_pipeline_directory, "fastQTL_output"))){dir.create(paste0(gtex_pipeline_directory, "fastQTL_output"))}
for(tissue in names(motrpac_gtex_map)){
        command_gtex_1 <- paste0("python3 fastqtl/python/run_FastQTL_threaded.py ", 
                                 "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz ", 
                                 paste0("log2-normalized-expression/log2-normalized-expression_", tissue, ".expression.bed.gz "), 
                                 paste0("log2-normalized-expression_", tissue, " "), 
                                 paste0("--covariates GTEx_Analysis_v8_eQTL_covariates/", motrpac_gtex_map_covariates[match(tissue, names(motrpac_gtex_map_covariates))], ".v8.covariates.txt", " "), 
                                 "--window 1e6 --chunks 100 --threads 16; ")
        command_gtex_2 <- paste0("python3 fastqtl/python/run_FastQTL_threaded.py ", 
                                 "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz ", 
                                 paste0("log2-normalized-expression/log2-normalized-expression_", tissue, ".expression.bed.gz "), 
                                 paste0("log2-normalized-expression_", tissue, " "), 
                                 paste0("--covariates GTEx_Analysis_v8_eQTL_covariates/", motrpac_gtex_map_covariates[match(tissue, names(motrpac_gtex_map_covariates))], ".v8.covariates.txt", " "), 
                                 "--window 1e6 --chunks 100 --threads 16 --permute 1000 10000; ")
        command_movefile_1 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".egenes.txt.gz "),
                                     paste0("fastQTL_output/log2-normalized-expression_", tissue, ".egenes.txt.gz; "))
        command_movefile_2 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".allpairs.txt.gz "),
                                     paste0("fastQTL_output/log2-normalized-expression_", tissue, ".allpairs.txt.gz;"))
        command <- paste0(command_changedir, command_gtex_1, command_gtex_2, command_movefile_1, command_movefile_2)
        command <- paste0(command_changedir, command_gtex_1, command_movefile_1, command_movefile_2)
        system(command)      
}
        
#offset all expected count and tpm by 1?
for(tissue in names(motrpac_gtex_map)){
        print(tissue)
        
        #read files
        ec <- as.data.frame(fread(paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_", tissue, "_expected_count.gct.gz")))
        tpm <- as.data.frame(fread(paste0("~/repos/gtex-pipeline/expression_data/GTEx_Analysis_v8_", tissue, "_tpm.gct.gz")))
        
        #modify files
        ec[,-1] <- ec[,-1] + 1
        tpm[,-1] <- tpm[,-1] + 1
        ec <- as.data.table(ec)
        tpm <- as.data.table(tpm)
        
        #write files
        writeLines(paste0("#1.2\n", nrow(ec), "\t", ncol(ec)-1),
                   con = paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count_OFFSET1.gct"))
        fwrite(ec, paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count_OFFSET1.gct"), sep = "\t", 
               append = T, col.names = T)
        system(paste0("cd ",gtex_pipeline_directory, "expression_data/; gzip -f ", "GTEx_Analysis_v8_", tissue, "_expected_count_OFFSET1.gct"))
        
        
        writeLines(paste0("#1.2\n", nrow(tpm), "\t", ncol(tpm)-1),
                   con = paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm_OFFSET1.gct"))
        fwrite(tpm, paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm_OFFSET1.gct"), sep = "\t", 
               append = T, col.names = T)
        system(paste0("cd ",gtex_pipeline_directory, "expression_data/; gzip -f ", "GTEx_Analysis_v8_", tissue, "_tpm_OFFSET1.gct"))

}

#now rerun normalization script with offset
command_changedir <- paste0("source ~/.bash_profile; cd ", gtex_pipeline_directory, "; ")
if(!dir.exists(paste0(gtex_pipeline_directory, "log2-normalized-expression"))){dir.create(paste0(gtex_pipeline_directory, "log2-normalized-expression"))}
for(tissue in names(motrpac_gtex_map)){
        # foreach(tissue = names(motrpac_gtex_map)) %dopar% {
        
        print(tissue)
        
        command_gtex <- paste0("python3 qtl/src/eqtl_prepare_expression.py ", 
                               paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_tpm_OFFSET1.gct.gz "), 
                               paste0(gtex_pipeline_directory, "expression_data/GTEx_Analysis_v8_", tissue, "_expected_count_OFFSET1.gct.gz "), 
                               "gencode.v26.GRCh38.genes.gtf ", 
                               "sample_participant_lookup.txt ", 
                               "GTEx_Analysis_vcf_chr_list.txt ", 
                               paste0("log2-normalized-expression_", tissue, "_OFFSET1 "), 
                               "--tpm_threshold 0.1 ", 
                               "--count_threshold 6 ", 
                               "--sample_frac_threshold 0.2 ", 
                               "--normalization_method tmm; ")
        command_movefile_1 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, "_OFFSET1.expression.bed.gz.tbi "), 
                                     paste0("log2-normalized-expression/log2-normalized-expression_", tissue, "_OFFSET1.expression.bed.gz.tbi; "))
        command_movefile_2 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, "_OFFSET1.expression.bed.gz "), 
                                     paste0("log2-normalized-expression/log2-normalized-expression_", tissue, "_OFFSET1.expression.bed.gz;"))
        command <- paste0(command_changedir, command_gtex, command_movefile_1, command_movefile_2)
        system(command)      
}

#now run tensorqtl on these
library(foreach)
library(doParallel)
library(parallel)

#initialize parallelization
if(!exists("cl")){
        cl <- makeCluster(4, outfile="")
        registerDoParallel(cl)
}
getDoParWorkers()

motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map, pattern = "-", replacement = "_")
motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map_covariates, pattern = " ", replacement = "_")
motrpac_gtex_map_covariates <- gsub(x = motrpac_gtex_map_covariates, pattern = "___", replacement = "_")
all(paste0(motrpac_gtex_map_covariates, ".v8.covariates.txt") %in% list.files("~/repos/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates"))
if(!dir.exists(paste0(gtex_pipeline_directory, "tensorQTL_output"))){dir.create(paste0(gtex_pipeline_directory, "tensorQTL_output"))}
# for(tissue in names(motrpac_gtex_map)){
command_changedir <- paste0("source ~/.bash_profile; cd ", gtex_pipeline_directory, "; ")
foreach(tissue = names(motrpac_gtex_map)) %dopar% {
        print(tissue)
        command_gtex_1 <- paste0("python3 -m tensorqtl ",
                                 "GTEx_v8 ",
                                 paste0("log2-normalized-expression/log2-normalized-expression_", tissue,".expression.bed.gz "), #whoops! just redoing hypothalamus expression lol
                                 paste0("log2-normalized-expression_", tissue," "),
                                 paste0("--covariates GTEx_Analysis_v8_eQTL_covariates/", motrpac_gtex_map_covariates[names(motrpac_gtex_map_covariates) == tissue], ".v8.covariates.txt"," "),
                                 "--mode cis ",
                                 "--window 1000000 --maf_threshold 0.01 &> ", paste0("tensorQTL_output/", tissue, "_permutation_screen-out.txt"),"; ")
        command_gtex_2 <- paste0("python3 -m tensorqtl ",
                                 "GTEx_v8 ",
                                 paste0("log2-normalized-expression/log2-normalized-expression_", tissue,".expression.bed.gz "),
                                 paste0("log2-normalized-expression_", tissue," "),
                                 paste0("--covariates GTEx_Analysis_v8_eQTL_covariates/", motrpac_gtex_map_covariates[names(motrpac_gtex_map_covariates) == tissue], ".v8.covariates.txt"," "),
                                 "--mode cis_nominal ",
                                 "--window 1000000 --maf_threshold 0.01 &> ", paste0("tensorQTL_output/", tissue, "_nominal_screen-out.txt"),"; ")
        command_movefiles_1 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".tensorQTL.cis_nominal.log "), 
                                     paste0("tensorQTL_output/log2-normalized-expression_", tissue, ".tensorQTL.cis_nominal.log; "))
        command_movefiles_2 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".tensorQTL.cis.log "), 
                                    paste0("tensorQTL_output/log2-normalized-expression_", tissue, ".tensorQTL.cis.log; "))
        command_movefiles_3 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".cis_qtl.txt.gz "), 
                                      paste0("tensorQTL_output/log2-normalized-expression_", tissue, ".cis_qtl.txt.gz; "))
        command_movefiles_4 <- paste0("mv ", paste0("log2-normalized-expression_", tissue, ".cis_qtl_pairs.chr",c(1:22, "X"),".parquet "), 
                                      paste0("tensorQTL_output/log2-normalized-expression_", tissue, ".cis_qtl_pairs.chr",c(1:22, "X"),".parquet; "))
        command_movefiles <- c(command_movefiles_1, command_movefiles_2, command_movefiles_3, command_movefiles_4)
        command_movefiles <- paste0(command_movefiles, collapse = "")
        command <- paste0(command_changedir, command_gtex_1, command_gtex_2, command_movefiles)
        system(command)      
        
}
