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

#custom ldsc annotations
if(!exists("cluster_membership")){
  load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
}

#summon mapping of rat -> human genes
geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")
DE_genes <- cluster_membership[cluster_membership$ome == "TRNSCRPT",]
cluster_ids <- sort(unique(DE_genes$cluster))
DE_genes <- lapply(cluster_ids, function(cli) DE_genes[DE_genes$cluster == cli,])
for(cli in cluster_ids){
  DE_genes[[cli]]$human_ensembl_gene <- geneID_map$human_ensembl_gene[match(DE_genes[[cli]]$feature_ID, geneID_map$feature_ID)]
  DE_genes[[cli]] <- DE_genes[[cli]][-which(is.na(DE_genes[[cli]]$human_ensembl_gene)),]
}
names(DE_genes) <- cluster_ids
str(DE_genes)

#confirm ENSG coordinate mapping
ENSG_coord <- read.table(file = "~/repos/ldsc/ENSG_coord.txt", header = T)
100-sapply(cluster_ids, function(cli) length(setdiff(DE_genes[[cli]]$human_ensembl_gene, ENSG_coord$GENE)) / length(DE_genes[[cli]]$human_ensembl_gene) * 100)

sum(sapply(cluster_ids, function(cli) length(unique(DE_genes[[cli]]$human_ensembl_gene))))
length(unique(unlist(sapply(cluster_ids, function(cli) ((DE_genes[[cli]]$human_ensembl_gene))))))

#write gene sets to file
ldsc_directory <- "~/repos/ldsc/"
geneset_directory <- "~/repos/ldsc/custom_genesets/"

if(!dir.exists(geneset_directory)){dir.create(geneset_directory)}
for(cli in cluster_ids){
  sink(paste0(geneset_directory, "/cluster_", cli, ".GeneSet"))
  cat(paste0(unique(DE_genes[[cli]]$human_ensembl_gene), "\n"), sep = "")
  sink()
}

for(cli in cluster_ids){
  print(length(unique(DE_genes[[cli]]$human_ensembl_gene)))
}

ncl <- length(cluster_ids)
sapply(1:(ncl-1), function(cl1) sapply((cl1+1):ncl, function(cl2) 
  length(intersect(unique(DE_genes[[cl1]]$human_ensembl_gene), unique(DE_genes[[cl2]]$human_ensembl_gene))) / 
  length(union(unique(DE_genes[[cl1]]$human_ensembl_gene), unique(DE_genes[[cl2]]$human_ensembl_gene)))))

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){
    cat(paste0("(chr: ", cri, ", clus: ", cli, ") "))
    command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets/cluster_", cli, ".GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "cluster_", cli, ".chr_", cri,".annot")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command, ignore.stdout = T)
  }
}

#thicken up those .annots
subset_baseline_to_BIM <- T
# foreach(cri=1:n_chromosomes) %dopar% {
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
# for(cri in 1:n_chromosomes){
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("~/repos/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("~/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  if(subset_baseline_to_BIM){
    SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
    baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
    baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))
  } else {
    SNP_inds <- match(bim_file$SNP, baseline_annot$SNP)
  }
  
  missing_SNP_inds <- which(is.na(SNP_inds))
  
  # intersect(baseline_annot$SNP, bim_file$SNP)
  for(cli in cluster_ids){
    cat(paste0(cli, " "))
    
    cluster_code <- fread(paste0("~/repos/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
    cluster_annot <- cbind(bim_file, cluster_code)
    # cluster_annot <- merge(baseline_annot, cluster_annot, by = "SNP")
    # sum(cluster_annot$BP.x != cluster_annot$BP.y)
    # cor(cluster_annot$CM.x, cluster_annot$CM.y)
    
    if(subset_baseline_to_BIM){
      
      merged_annot <- cluster_annot
      merged_annot <- cbind(merged_annot, base = 1)
      merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
      merged_annot[baseline_present_SNPs,7:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
      
      colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_", cli)
      colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]
      # tail(merged_annot, n = 100)
      
    } else {
      
      merged_annot <- baseline_annot
      merged_annot$cluster_TEMP <- 0
      merged_annot$cluster_TEMP[SNP_inds[!is.na(SNP_inds)]] <- cluster_annot$ANNOT[!is.na(SNP_inds)]
      
      missing_SNPS <- cluster_annot[missing_SNP_inds,]
      missing_SNPS <- cbind(missing_SNPS[,-"ANNOT"], base = rep(1, nrow(missing_SNPS)), 
                            matrix(0, nrow = nrow(missing_SNPS), ncol = ncol(baseline_annot) - 5), missing_SNPS[,"ANNOT"])
      merged_annot <- rbind(merged_annot, missing_SNPS, use.names = F)
      colnames(merged_annot)[colnames(merged_annot) == "cluster_TEMP"] <- paste0("Cluster_", cli)
    }
    
    merged_annot$CM <- as.character(merged_annot$CM)
    merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
    
    #write to disk
    fwrite(merged_annot, paste0("~/repos/ldsc/custom_genesets/thickCluster_", cli, ".chr_", cri, ".annot"), sep = "\t")
    
  }
}

#process ld files
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()
# detach(“package:data.table”, unload=TRUE) 
foreach(cri=1:n_chromosomes) %dopar% {
# for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){ #redundant across clusters
    print(paste0("now processing cluster ", cli))
    # command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --thin-annot --out ", geneset_directory, "cluster_", cli, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --out ", geneset_directory, "cluster_", cli, ".chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}  

for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){ 
    file.rename(from = paste0("~/repos/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"), 
                to = paste0("~/repos/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot"))
    file.rename(from = paste0("~/repos/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot.gz"), 
                to = paste0("~/repos/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot.gz"))
    file.rename(from = paste0("~/repos/ldsc/custom_genesets/thickCluster_", cli, ".chr_", cri, ".annot"), 
                to = paste0("~/repos/ldsc/custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
  }
}

#python2 ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr baseline/baseline. --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr 1000G_frq/1000G.mac5eur. --out BMI_baseline

#now partition heritability for 2010 BMI GWAS
for(cli in cluster_ids){ #redundant across clusters
  print(paste0("now processing cluster ", cli))
  command_ldsc <- paste0("python2 ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr custom_genesets/cluster_", cli, ".chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out output/BMI_Cluster_", cli)
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#let's get the gwas summary stats into a format that munge_sumstats.py can process
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]
for(gi in 1:length(gwas_summary_files)){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")"))
  sumstats_file <- fread(paste0(gwas_dir, gwas_summary_files[gi]))
  
  # sumstats_file_subset <- sumstats_file[,c("variant_id", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "n_cases")]
  # colnames(sumstats_file_subset) <- c("MarkerName", "Allele1", "Allele2", "Freq.Allele1.HapMapCEU",	"b", "se", "p", "N")
  # if(all(is.na(sumstats_file_subset[,N]))){
  #   sumstats_file_subset <- sumstats_file_subset[,-"N"]
  # }
  # if(all(is.na(sumstats_file_subset[,se]))){
  #   sumstats_file_subset <- sumstats_file_subset[,-"se"]
  # }
  # if(all(is.na(sumstats_file_subset[,b]))){
  #   sumstats_file_subset <- sumstats_file_subset[,-"b"]
  # }
  
  sumstats_file_subset <- sumstats_file[,c("variant_id", "effect_allele", "non_effect_allele", "frequency", "pvalue", "sample_size")]
  colnames(sumstats_file_subset) <- c("MarkerName", "Allele1", "Allele2", "Freq.Allele1.HapMapCEU", "p", "N")
  fwrite(sumstats_file_subset, file = paste0(ldsc_directory, "gwas_sumstats/", gwas_summary_files[gi]), sep = " ")
}

#now let's munge those sumstats
for(gi in 1:length(gwas_summary_files)){
# for(gi in 1:1){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  command_ldsc <- paste0("python2 munge_sumstats.py --sumstats ~/repos/ldsc/gwas_sumstats/", gwas_summary_files[gi], " --merge-alleles w_hm3.snplist --out gwas_sumstats/proper_format/", gwas_summary_files[gi]," --a1-inc")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#now compute partitioned heritabilities
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
foreach(gi=1:length(gwas_summary_files)) %dopar% {
# for(gi in 1:length(gwas_summary_files)){
  for(cli in cluster_ids){ #redundant across clusters
    cat(paste0(" (", cli, ": ", gi, " / ", length(gwas_summary_files), ")")) 
    command_ldsc <- paste0("python2 ldsc.py --h2 ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_", cli, ".chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out output/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_", cli)
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}

#### joint / conditional analysis ####

#create master .annot
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("~/repos/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("~/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
  baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
  baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))

  missing_SNP_inds <- which(is.na(SNP_inds))
  
  #initialize at 1st cluster
  cli = cluster_ids[1]
  cluster_code <- fread(paste0("~/repos/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot"))
  cluster_annot <- cbind(bim_file, cluster_code)
  merged_annot <- cluster_annot
  colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_", cli)
  
  #add in remaining clusters
  for(cli in cluster_ids[-1]){
    cat(paste0(cli, " "))
    cluster_code <- fread(paste0("~/repos/ldsc/custom_genesets/thincluster_", cli, ".chr_", cri, ".annot"))
    merged_annot <- cbind(merged_annot, cluster_code)
    colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_", cli)
  }
   
  merged_annot <- cbind(merged_annot, base = 1)
  merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
  merged_annot[baseline_present_SNPs,21:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
  colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]

  merged_annot$CM <- as.character(merged_annot$CM)
  merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
  
  #write to disk
  fwrite(merged_annot, paste0("~/repos/ldsc/custom_genesets/cluster_1-15.chr_", cri, ".annot"), sep = "\t")
  
}

#process ld files
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()
# detach(“package:data.table”, unload=TRUE) 
foreach(cri=1:n_chromosomes) %dopar% {
  print(paste0("now processing cr: ", cri))
  # command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --thin-annot --out ", geneset_directory, "cluster_", cli, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "cluster_1-15.chr_", cri,".annot --out ", geneset_directory, "cluster_1-15.chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#now compute partitioned heritabilities
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
foreach(gi=1:length(gwas_summary_files)) %dopar% {
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  command_ldsc <- paste0("python2 ldsc.py --h2 ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_1-15.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out output/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#### pairwise genetic correlations ####

ldsc_directory <- "~/repos/ldsc/"

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#gwas sumstats locations
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]

#let's get the gwas summary stats appropriately signed
naughty_focal_alleles_traits <- numeric(0)
for(gi in 1:length(gwas_summary_files)){
  cat(paste0(gi, " "))
  orig_sumstats <- fread(paste0(gwas_dir, c(gwas_summary_files[gi])))
  unsigned_sumstats <- fread(paste0(ldsc_directory, "gwas_sumstats/proper_format/", c(gwas_summary_files[gi]), ".sumstats.gz"))
  a1s <- orig_sumstats$effect_allele[match(unsigned_sumstats$SNP, orig_sumstats$variant_id)]
  matching_focal_alleles <- unsigned_sumstats$A1 == a1s | is.na(unsigned_sumstats$A1 == a1s)  | unsigned_sumstats$A1 == ""
  if(!all(matching_focal_alleles)){
    print(paste0("error, error, focal alleles do not match for: ", gi))
    naughty_focal_alleles_traits <- c(naughty_focal_alleles_traits, gi)
    next()
  }
  zscore_signs <- sign(orig_sumstats$zscore[match(unsigned_sumstats$SNP, orig_sumstats$variant_id)])
  signed_zscores <- zscore_signs * unsigned_sumstats$Z
  signed_sumstats <- as.data.frame(unsigned_sumstats)
  signed_sumstats$Z <- signed_zscores
  signed_sumstats[signed_sumstats == ""] <- NA
  signed_sumstats$N[!is.na(signed_sumstats$N)] <- paste0(signed_sumstats$N[!is.na(signed_sumstats$N)], ".000")
  fwrite(x = signed_sumstats, sep = "\t", 
         file = paste0(ldsc_directory, "gwas_sumstats/proper_format/signed/", c(gwas_summary_files[gi]), ".sumstats.gz"))
}

#run the analyses across gwas
for(gi in 1:length(gwas_summary_files)){
  # for(gi in 1:1){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  command_ldsc <- paste0("python2 ldsc.py --rg ", paste0("gwas_sumstats/proper_format/signed/", c(gwas_summary_files[gi], gwas_summary_files[-gi]), ".sumstats.gz", collapse = ","),
                         " --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out output/signed_gcor/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_pairwise-Gcorrs")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### partitioned heritabilities with cts flag ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#get list of all genes in transcriptome DEA
if(!exists("genes_tested_in_transcriptome_DEA")){
  load("~/data/smontgom/genes_tested_in_transcriptome_DEA.RData")
}

#summon mapping of rat -> human genes
geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")

human_ensembl_genes_in_DEA <- geneID_map$human_ensembl_gene[match(genes_tested_in_transcriptome_DEA, geneID_map$feature_ID)]
human_ensembl_genes_in_DEA <- human_ensembl_genes_in_DEA[-which(is.na(human_ensembl_genes_in_DEA))]

#confirm ENSG coordinate mapping
ENSG_coord <- read.table(file = "~/repos/ldsc/ENSG_coord.txt", header = T)
100-length(setdiff(human_ensembl_genes_in_DEA, ENSG_coord$GENE)) / length(human_ensembl_genes_in_DEA) * 100

#write gene sets to file
ldsc_directory <- "~/repos/ldsc/"
geneset_directory <- "~/repos/ldsc/custom_genesets/"
sink(paste0(geneset_directory, "/cluster_control.GeneSet"))
cat(paste0(human_ensembl_genes_in_DEA, "\n"), sep = "")
sink()

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  cat(paste0("(chr: ", cri, ") "))
  command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets/cluster_control.GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "cluster_control.chr_", cri,".annot")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command, ignore.stdout = T)
}

#thicken up those .annots
subset_baseline_to_BIM <- T
# foreach(cri=1:n_chromosomes) %dopar% {
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("~/repos/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("~/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  if(subset_baseline_to_BIM){
    SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
    baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
    baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))
  } else {
    SNP_inds <- match(bim_file$SNP, baseline_annot$SNP)
  }
  
  missing_SNP_inds <- which(is.na(SNP_inds))
    
  cluster_code <- fread(paste0("~/repos/ldsc/custom_genesets/cluster_control.chr_", cri,".annot"))
  cluster_annot <- cbind(bim_file, cluster_code)
  # cluster_annot <- merge(baseline_annot, cluster_annot, by = "SNP")
  # sum(cluster_annot$BP.x != cluster_annot$BP.y)
  # cor(cluster_annot$CM.x, cluster_annot$CM.y)
  
  if(subset_baseline_to_BIM){
    
    merged_annot <- cluster_annot
    merged_annot <- cbind(merged_annot, base = 1)
    merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
    merged_annot[baseline_present_SNPs,7:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
    
    colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_Control")
    colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]
    # tail(merged_annot, n = 100)
    
  } else {
    
    merged_annot <- baseline_annot
    merged_annot$cluster_TEMP <- 0
    merged_annot$cluster_TEMP[SNP_inds[!is.na(SNP_inds)]] <- cluster_annot$ANNOT[!is.na(SNP_inds)]
    
    missing_SNPS <- cluster_annot[missing_SNP_inds,]
    missing_SNPS <- cbind(missing_SNPS[,-"ANNOT"], base = rep(1, nrow(missing_SNPS)), 
                          matrix(0, nrow = nrow(missing_SNPS), ncol = ncol(baseline_annot) - 5), missing_SNPS[,"ANNOT"])
    merged_annot <- rbind(merged_annot, missing_SNPS, use.names = F)
    colnames(merged_annot)[colnames(merged_annot) == "cluster_TEMP"] <- paste0("Cluster_Control")
  }
  
  merged_annot$CM <- as.character(merged_annot$CM)
  merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
  
  #write to disk
  fwrite(merged_annot, paste0("~/repos/ldsc/custom_genesets/cluster_control.chr_", cri,".annot"), sep = "\t")

}


#process ld files
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()
# detach(“package:data.table”, unload=TRUE) 
foreach(cri=1:n_chromosomes) %dopar% {
  print(paste0("(", cri, ")"))
  # command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "thickCluster_", cli, ".chr_", cri,".annot --thin-annot --out ", geneset_directory, "cluster_", cli, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "cluster_control.chr_", cri,".annot --out ", geneset_directory, "cluster_control.chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}  


#create .ldcts files
ldcts <- cbind(paste0("Cluster_", 1:15, "  "), paste0("custom_genesets/cluster_", 1:15, ".chr_"), ",", paste0("custom_genesets/cluster_control.chr_"))
ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))

sink(paste0(ldsc_directory, "cluster_cts.ldcts"))
cat(paste0(ldcts, "\n"), sep = "")
sink()

#let's convert all of my .annot's to .annot.gz's???
annot_files <- list.files(paste0(ldsc_directory, "custom_genesets"))
annot_files <- annot_files[grep("annot", annot_files)]
annot_files_GZ <- annot_files[grep("gz", annot_files)]
annot_files_GZ <- stringr::str_remove_all(annot_files_GZ, ".gz")
annot_files_noGZ <- annot_files[-grep("gz", annot_files)]
annot_files_to_zip <- setdiff(annot_files_noGZ, annot_files_GZ)
for(i in annot_files_to_zip){
  print(which(annot_files_to_zip == i) / length(annot_files_to_zip) * 100)
  command_ldsc <- paste0("gzip -k ", paste0("custom_genesets/", i))
  command <- paste0(command_changedir, command_ldsc)
  system(command)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#hmmm still not working... what if we try to get everything into the exact format of the cts tutorial??
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){
    cat(paste0("(", cri, ", ", cli, ")"))
    #make cluster annotation thin
    thick_annot <- fread(paste0(ldsc_directory, "custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
    target_col <- which(colnames(thick_annot) == paste0("Cluster_", cli))
    thin_annot <- thick_annot[,..target_col]
    colnames(thin_annot) <- "ANNOT"
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/cluster_", cli, ".chr_", cri, ".annot"))
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/cluster_", cli, ".chr_", cri, ".annot.gz"))
  }
}
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids[1]){
    cat(paste0("(", cri, ", ", cli, ")"))
    #make baseline annotation its own thing
    thick_annot <- fread(paste0(ldsc_directory, "custom_genesets/cluster_", cli, ".chr_", cri, ".annot"))
    target_col <- which(colnames(thick_annot) != paste0("Cluster_", cli))
    baseline_annot <- thick_annot[,..target_col]
    fwrite(baseline_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/baseline.chr_", cri, ".annot"), sep = "\t")
    fwrite(baseline_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/baseline.chr_", cri, ".annot.gz"), sep = "\t")
  }
}
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
# for(cri in 1:n_chromosomes){
  cat(paste0("(", cri, ", ", cli, ")"))
  #make cluster annotation thin
  thick_annot <- fread(paste0(ldsc_directory, "custom_genesets/cluster_control.chr_", cri, ".annot"))
  target_col <- which(colnames(thick_annot) == paste0("Cluster_Control"))
  thin_annot <- thick_annot[,..target_col]
  colnames(thin_annot) <- "ANNOT"
  fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/cluster_control.chr_", cri, ".annot"))
  fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets/cts_thin_annot/cluster_control.chr_", cri, ".annot.gz"))
}


if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
#process ld files?? I guess ugh
foreach(cri=1:n_chromosomes) %dopar% {
  
  cat(paste0("(", cri, ",", "b", ")"))
  command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --annot ", geneset_directory, "cts_thin_annot/baseline.chr_", cri,".annot --out ", geneset_directory, "cts_thin_annot/baseline.chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
  
  for(cli in cluster_ids){
    cat(paste0("(", cri, ",", cli, ")"))
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --thin-annot --annot ", geneset_directory, "cts_thin_annot/cluster_", cli, ".chr_", cri,".annot --out ", geneset_directory, "cts_thin_annot/cluster_", cli, ".chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
    }
}

#get control annotation processed for ld scores
foreach(cri=1:n_chromosomes) %dopar% {
  cat(paste0("(", cri, ",", cli, ")"))
  command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --thin-annot --annot ", geneset_directory, "cts_thin_annot/cluster_control.chr_", cri,".annot --out ", geneset_directory, "cts_thin_annot/cluster_control.chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}

#and then remake that cluster_cts.ldcts file
ldcts <- cbind(paste0("Cluster_", 1:15, "  "), paste0("custom_genesets/cts_thin_annot/cluster_", 1:15, ".chr_"), ",", paste0("custom_genesets/cts_thin_annot/cluster_control.chr_"))
ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))

sink(paste0(ldsc_directory, "cluster_cts_2.ldcts"))
cat(paste0(ldcts, "\n"), sep = "")
sink()

#now compute partitioned heritabilities??????
if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

foreach(gi=1:length(gwas_summary_files)) %dopar% {
# for(gi in 1:1){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_1-15.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 200 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cts_thin_annot/baseline.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts_2.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 200 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 1000 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_1.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}


#### hmm what if I did the DE genes per sex at timepoint 4? or across the other timepoints? or tissues at all timepoints? ####
if(!exists("rna_dea")){
  load('~/data/smontgom/dea/transcript_rna_seq_20210126.RData')
  rna_dea <- as.data.frame(transcript_rna_seq$timewise_dea)
}

sapply(unique(rna_dea$tissue), function(tiss) length(unique(rna_dea$feature_ID[rna_dea$tissue == tiss & rna_dea$selection_fdr < 0.05])))
sapply(unique(rna_dea$tissue), function(tiss1) sapply(unique(rna_dea$tissue), function(tiss2){
  t1g <- unique(rna_dea$feature_ID[rna_dea$tissue == tiss1 & rna_dea$selection_fdr < 0.05])
  t2g <- unique(rna_dea$feature_ID[rna_dea$tissue == tiss2 & rna_dea$selection_fdr < 0.05])
  length(intersect(t1g, t2g)) / length(union(t1g, t2g))
}))
 
sapply(unique(rna_dea$sex), function(sx) length(unique(rna_dea$feature_ID[rna_dea$sex == sx & rna_dea$selection_fdr < 0.05])))
length(intersect(unique(rna_dea$feature_ID[rna_dea$sex == "male" & rna_dea$selection_fdr < 0.05]), 
       unique(rna_dea$feature_ID[rna_dea$sex == "female" & rna_dea$selection_fdr < 0.05]))) / 
  length(union(unique(rna_dea$feature_ID[rna_dea$sex == "male" & rna_dea$selection_fdr < 0.05]), 
                   unique(rna_dea$feature_ID[rna_dea$sex == "female" & rna_dea$selection_fdr < 0.05])))

sapply(unique(rna_dea$comparison_group), function(time) length(unique(rna_dea$feature_ID[rna_dea$comparison_group == time & rna_dea$adj_p_value < 0.05])))
sapply(unique(rna_dea$comparison_group), function(t1) sapply(unique(rna_dea$comparison_group), function(t2){
  t1g <- unique(rna_dea$feature_ID[rna_dea$comparison_group == t1 & rna_dea$adj_p_value < 0.05])
  t2g <- unique(rna_dea$feature_ID[rna_dea$comparison_group == t2 & rna_dea$adj_p_value < 0.05])
  length(intersect(t1g, t2g)) / length(union(t1g, t2g))
}))

sapply(unique(rna_dea$sex), function(sx) sapply(unique(rna_dea$comparison_group), function(time) 
  length(unique(rna_dea$feature_ID[rna_dea$sex == sx & rna_dea$comparison_group == time & rna_dea$adj_p_value < 0.05]))))
sapply(unique(rna_dea$sex), function(sx1) sapply(unique(rna_dea$comparison_group), function(time1) 
sapply(unique(rna_dea$sex), function(sx2) sapply(unique(rna_dea$comparison_group), function(time2){
  t1g <- unique(rna_dea$feature_ID[rna_dea$sex == sx1 & rna_dea$comparison_group == time1 & rna_dea$adj_p_value < 0.05])
  t2g <- unique(rna_dea$feature_ID[rna_dea$sex == sx2 & rna_dea$comparison_group == time2 & rna_dea$adj_p_value < 0.05])
  length(intersect(t1g, t2g)) / length(union(t1g, t2g))
}))))

# OK so new gene sets! let's do tissue-specific all-time, time-specific, time+sex-specific, & all DE genes
genesets <- lapply(setNames(unique(rna_dea$tissue),unique(rna_dea$tissue)), function(tiss) data.frame(feature_ID = unique(rna_dea$feature_ID[rna_dea$tissue == tiss & rna_dea$selection_fdr < 0.05])))
genesets <- c(genesets, lapply(setNames(unique(rna_dea$comparison_group), unique(rna_dea$comparison_group)), 
                               function(time) data.frame(feature_ID = unique(rna_dea$feature_ID[rna_dea$comparison_group == time & rna_dea$adj_p_value < 0.05]))))
timesex <- paste0(unique(rna_dea$sex), "-", unique(rna_dea$comparison_group))
genesets <- c(genesets, lapply(setNames(timesex, timesex), function(ts) 
                        data.frame(feature_ID = unique(rna_dea$feature_ID[rna_dea$sex == strsplit(ts, "-")[[1]][1] & rna_dea$comparison_group == strsplit(ts, "-")[[1]][2] & rna_dea$adj_p_value < 0.05]))))
genesets <- c(genesets, list(all_DE = data.frame(feature_ID = unique(rna_dea$feature_ID[rna_dea$selection_fdr < 0.05]))))


#alternatively, what if we do the things that nicole and dan proposed?
# ie cluster 7, 15 + gastroc, heart, liver
if(!exists("cluster_membership")){
  load("~/data/dea_clustering_0.1-FDR-ftest_kmeans-15.RData")
}

#summon mapping of rat -> human genes
geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")
DE_genes <- cluster_membership[cluster_membership$ome == "TRNSCRPT",]
cluster_ids <- c(7, 15)
DE_genes <- lapply(cluster_ids, function(cli) DE_genes[DE_genes$cluster == cli,])
for(cli in cluster_ids){
  DE_genes[[cli]]$human_ensembl_gene <- geneID_map$human_ensembl_gene[match(DE_genes[[cli]]$feature_ID, geneID_map$feature_ID)]
  DE_genes[[cli]] <- DE_genes[[cli]][-which(is.na(DE_genes[[cli]]$human_ensembl_gene)),]
}
names(DE_genes) <- cluster_ids
str(DE_genes)
geneset_info <- expand.grid(cluster_ids, c("t55-gastrocnemius", "t68-liver", "t58-heart"))
geneset_names <- paste0("Cluster_", trimws(apply(geneset_inf, 1, paste0, collapse = "-")), "")
rownames(geneset_info) <- geneset_names
colnames(geneset_info) <- c("cluster", "tissue")
geneset_info$tissue_common_name <- MotrpacBicQC::bic_animal_tissue_code$abbreviation[match(geneset_info$tissue, MotrpacBicQC::bic_animal_tissue_code$tissue_name_release)]
genesets <- lapply(setNames(geneset_names, geneset_names), function(group){ 
    all_genes <- DE_genes[[as.character(geneset_info[group, "cluster"])]]
    data.frame(feature_ID = all_genes$feature_ID[all_genes$tissue == geneset_info[group, "tissue_common_name"]])
  }
)
sapply(genesets, function(x1) sapply(genesets, function(x2) length(intersect(x1[,1], x2[,1])) / length(union(x1[,1], x2[,1])) * 100))

#summon mapping of rat -> human genes
geneID_map <- read.table("~/data/smontgom/motrpac_geneID_map.txt")
cluster_ids <- names(genesets)
DE_genes <- genesets
for(cli in cluster_ids){
  DE_genes[[cli]]$human_ensembl_gene <- geneID_map$human_ensembl_gene[match(DE_genes[[cli]]$feature_ID, geneID_map$feature_ID)]
  DE_genes[[cli]] <- DE_genes[[cli]][-which(is.na(DE_genes[[cli]]$human_ensembl_gene)),]
}
names(DE_genes) == cluster_ids
str(DE_genes)

#confirm ENSG coordinate mapping
ENSG_coord <- read.table(file = "~/repos/ldsc/ENSG_coord.txt", header = T)
100-sapply(cluster_ids, function(cli) length(setdiff(DE_genes[[cli]]$human_ensembl_gene, ENSG_coord$GENE)) / length(DE_genes[[cli]]$human_ensembl_gene) * 100)

sum(sapply(cluster_ids, function(cli) length(unique(DE_genes[[cli]]$human_ensembl_gene))))
length(unique(unlist(sapply(cluster_ids, function(cli) ((DE_genes[[cli]]$human_ensembl_gene))))))

#write gene sets to file
ldsc_directory <- "~/repos/ldsc/"
geneset_directory <- "~/repos/ldsc/custom_genesets_2/"

if(!dir.exists(geneset_directory)){dir.create(geneset_directory)}
for(cli in cluster_ids){
  sink(paste0(geneset_directory, "/cluster_", cli, ".GeneSet"))
  cat(paste0(unique(DE_genes[[cli]]$human_ensembl_gene), "\n"), sep = "")
  sink()
}

for(cli in cluster_ids){
  print(length(unique(DE_genes[[cli]]$human_ensembl_gene)))
}

ncl <- length(cluster_ids)
sapply(1:(ncl-1), function(cl1) sapply((cl1+1):ncl, function(cl2) 
  length(intersect(unique(DE_genes[[cl1]]$human_ensembl_gene), unique(DE_genes[[cl2]]$human_ensembl_gene))) / 
    length(union(unique(DE_genes[[cl1]]$human_ensembl_gene), unique(DE_genes[[cl2]]$human_ensembl_gene)))))

#generic commands
command_changedir <- paste0("source ~/.bash_profile; cd ", ldsc_directory, "; ")
command_initiateLDSC <- "source activate ldsc; "

#construct annotation files
windowsize <- "100000"
n_chromosomes <- 22

#make .annot files
foreach(cri=1:n_chromosomes) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){
    cat(paste0("(chr: ", cri, ", clus: ", cli, ") "))
    command_ldsc <- paste0("python2 make_annot.py --gene-set-file custom_genesets_2/cluster_", cli, ".GeneSet --gene-coord-file ENSG_coord.txt --windowsize ", windowsize, " --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri,".bim --annot-file ", geneset_directory, "cluster_", cli, ".chr_", cri,".annot")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command, ignore.stdout = T)
  }
}

#thicken up those .annots
subset_baseline_to_BIM <- T
# foreach(cri=1:n_chromosomes) %dopar% {
foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  cat(paste0("\n(chr: ", cri, "): "))
  baseline_annot <- fread(paste0("~/repos/ldsc/baseline/baseline.", cri, ".annot"))
  bim_file <- fread(paste0("~/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri, ".bim"), header = F)
  colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "minor_allele", "major_allele")
  bim_file <- bim_file[,-c("minor_allele", "major_allele")]
  
  if(subset_baseline_to_BIM){
    SNP_inds <- match(baseline_annot$SNP, bim_file$SNP)
    baseline_present_SNPs <- SNP_inds[!is.na(SNP_inds)]
    baseline_absent_SNPs <- which(is.na(match(bim_file$SNP, baseline_annot$SNP)))
  } else {
    SNP_inds <- match(bim_file$SNP, baseline_annot$SNP)
  }
  
  missing_SNP_inds <- which(is.na(SNP_inds))
  
  # intersect(baseline_annot$SNP, bim_file$SNP)
  for(cli in cluster_ids){
    cat(paste0(cli, " "))
    
    cluster_code <- fread(paste0("~/repos/ldsc/custom_genesets_2/cluster_", cli, ".chr_", cri, ".annot"))
    cluster_annot <- cbind(bim_file, cluster_code)
    # cluster_annot <- merge(baseline_annot, cluster_annot, by = "SNP")
    # sum(cluster_annot$BP.x != cluster_annot$BP.y)
    # cor(cluster_annot$CM.x, cluster_annot$CM.y)
    
    if(subset_baseline_to_BIM){
      
      merged_annot <- cluster_annot
      merged_annot <- cbind(merged_annot, base = 1)
      merged_annot <- cbind(merged_annot, matrix(0, nrow = nrow(merged_annot), ncol = ncol(baseline_annot) - 5))
      merged_annot[baseline_present_SNPs,7:ncol(merged_annot)] <- baseline_annot[!is.na(SNP_inds), 6:ncol(baseline_annot)]
      
      colnames(merged_annot)[colnames(merged_annot) == "ANNOT"] <- paste0("Cluster_", cli)
      colnames(merged_annot)[startsWith(colnames(merged_annot), "V")] <- colnames(baseline_annot)[6:ncol(baseline_annot)]
      # tail(merged_annot, n = 100)
      
    } else {
      
      merged_annot <- baseline_annot
      merged_annot$cluster_TEMP <- 0
      merged_annot$cluster_TEMP[SNP_inds[!is.na(SNP_inds)]] <- cluster_annot$ANNOT[!is.na(SNP_inds)]
      
      missing_SNPS <- cluster_annot[missing_SNP_inds,]
      missing_SNPS <- cbind(missing_SNPS[,-"ANNOT"], base = rep(1, nrow(missing_SNPS)), 
                            matrix(0, nrow = nrow(missing_SNPS), ncol = ncol(baseline_annot) - 5), missing_SNPS[,"ANNOT"])
      merged_annot <- rbind(merged_annot, missing_SNPS, use.names = F)
      colnames(merged_annot)[colnames(merged_annot) == "cluster_TEMP"] <- paste0("Cluster_", cli)
    }
    
    merged_annot$CM <- as.character(merged_annot$CM)
    merged_annot$CM[merged_annot$CM == "0"] <- "0.0"
    
    #write to disk
    fwrite(merged_annot, paste0("~/repos/ldsc/custom_genesets_2/thickCluster_", cli, ".chr_", cri, ".annot"), sep = "\t")
    
  }
}

for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){ 
    file.rename(from = paste0("~/repos/ldsc/custom_genesets_2/cluster_", cli, ".chr_", cri, ".annot"), 
                to = paste0("~/repos/ldsc/custom_genesets_2/thincluster_", cli, ".chr_", cri, ".annot"))
    file.rename(from = paste0("~/repos/ldsc/custom_genesets_2/cluster_", cli, ".chr_", cri, ".annot.gz"), 
                to = paste0("~/repos/ldsc/custom_genesets_2/thincluster_", cli, ".chr_", cri, ".annot.gz"))
    file.rename(from = paste0("~/repos/ldsc/custom_genesets_2/thickCluster_", cli, ".chr_", cri, ".annot"), 
                to = paste0("~/repos/ldsc/custom_genesets_2/cluster_", cli, ".chr_", cri, ".annot"))
  }
}

#get list of all genes in transcriptome DEA
if(!exists("genes_tested_in_transcriptome_DEA")){
  load("~/data/smontgom/genes_tested_in_transcriptome_DEA.RData")
}

foreach(cri=1:n_chromosomes, .packages = c("data.table")) %dopar% {
  # for(cri in 1:n_chromosomes){
  for(cli in cluster_ids){
    cat(paste0("(", cri, ", ", cli, ")"))
    #make cluster annotation thin
    thick_annot <- fread(paste0(ldsc_directory, "custom_genesets_2/cluster_", cli, ".chr_", cri, ".annot"))
    target_col <- which(colnames(thick_annot) == paste0("Cluster_", cli))
    thin_annot <- thick_annot[,..target_col]
    colnames(thin_annot) <- "ANNOT"
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets_2/cts_thin_annot/cluster_", cli, ".chr_", cri, ".annot"))
    fwrite(thin_annot, paste0(ldsc_directory, "custom_genesets_2/cts_thin_annot/cluster_", cli, ".chr_", cri, ".annot.gz"))
  }
}

if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(2, outfile="")
  registerDoParallel(cl)
}
#process ld files?? I guess ugh
foreach(cri=1:n_chromosomes) %dopar% {
  for(cli in cluster_ids){
    cat(paste0("(", cri, ",", cli, ")"))
    command_ldsc <- paste0("python2 ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.", cri," --ld-wind-cm 1 --thin-annot --annot ", geneset_directory, "cts_thin_annot/cluster_", cli, ".chr_", cri,".annot --out ", geneset_directory, "cts_thin_annot/cluster_", cli, ".chr_", cri, " --print-snps hapmap3_snps/hm.", cri,".snp")
    command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
    system(command)
  }
}


#and then remake that cluster_cts.ldcts file
ldcts <- cbind(paste0("cluster_", cluster_ids, "  "), paste0("custom_genesets_2/cts_thin_annot/cluster_", cluster_ids, ".chr_"), ",", paste0("custom_genesets_2/cts_thin_annot/cluster_control.chr_"))
ldcts <- sapply(1:nrow(ldcts), function(ri) paste0(ldcts[ri,], collapse = ""))
# ldcts_name <- list.files(ldsc_directory, pattern = "ldcts")

sink(paste0(ldsc_directory, "cluster_cts_4.ldcts"))
cat(paste0(ldcts, "\n"), sep = "")
sink()

#now compute partitioned heritabilities??????
#gwas sumstats locations
gwas_dir <- "~/data/smontgom/imputed_gwas_hg38_1.1/"
gwas_summary_files <- list.files(gwas_dir)
gwas_summary_files <- gwas_summary_files[-grep(gwas_summary_files, pattern = "README")]

if(exists("cl")){
  stopCluster(cl)
  remove(cl)
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

foreach(gi=1:length(gwas_summary_files)) %dopar% {
  # for(gi in 1:1){
  cat(paste0(" (", gi, " / ", length(gwas_summary_files), ")")) 
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_1-15.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 200 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets_2/cts_thin_annot/baseline.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts_4.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 200 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_Addenda")
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --n-blocks 1000 --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1-15")
  # command_ldsc <- paste0("python2 ldsc.py --h2-cts ", "gwas_sumstats/proper_format/", gwas_summary_files[gi], ".sumstats.gz --ref-ld-chr custom_genesets/cluster_1.chr_ --w-ld-chr weights_hm3_no_hla/weights. --overlap-annot --ref-ld-chr-cts cluster_cts.ldcts --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out output/cts/", strsplit(gwas_summary_files[gi], ".txt.gz")[[1]][1], "_Cluster_1")
  command <- paste0(command_changedir, command_initiateLDSC, command_ldsc)
  system(command)
}
