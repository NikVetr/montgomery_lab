
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
  cat(paste0(gi, " "))
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
save(gcor_mat, file = "~/data/smontgom/est_gcor_mat.RData")
