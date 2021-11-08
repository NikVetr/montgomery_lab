#nicole's other code
library(data.table)
library(ggplot2)
library(viridis)
library(MotrpacBicQC)
library(edgeR)
library(limma)
library(testit)
library(grid)
library(gridExtra)
# knitr::opts_chunk$set(echo = TRUE)
# knitr::opts_knit$set(root.dir = "/projects/motrpac/PASS1A-1C-06/RNA")
# source("/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/tools/get_fx.R") # for dl_read_gcp
# source("/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/integrative/outliers_covariates/pi1_cook_fx.R") # for run_deseq

gsutil = "~/google-cloud-sdk/bin/gsutil"
tmp = "/oak/stanford/groups/smontgom/nicolerg/tmp"
repeated_groups = c('acute_0.5h','acute_1h','control_0h')

group_cols = c(acute_0.5h='orange',
               acute_1h='yellow3',
               acute_24h='purple',
               acute_48h='brown',
               acute_4h='green3',
               acute_7h='blue',
               acute_IPEh='red',
               control_0.5h='gray40',
               control_0h='gray60',
               control_4h='gray20',
               control_7h='gray0',
               control_IPEh='gray80')

dl_format_pass1a_pheno = function(scratch, gsutil_path, parallel=F){
  dmaqc_metadata = 'gs://motrpac-internal-release2-results/pass1a-06/phenotype/motrpac_pass1a-06_pheno_viallabel-data.txt'
  dmaqc_dict = 'gs://motrpac-internal-release2-results/pass1a-06/phenotype/motrpac_pass1a_pheno_merged-dictionary.txt'
  
  # download and format phenotypic data
  dmaqc_metadata = dl_read_gcp(dmaqc_metadata, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
  cols = dl_read_gcp(dmaqc_dict, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
  old_cols = colnames(dmaqc_metadata)
  new_cols = tolower(cols[match(old_cols, BICUniqueID), FullName])
  colnames(dmaqc_metadata) = new_cols # this isn't perfect, but we don't care about the columns it doesn't work for for now
  
  # make some variables human-readable
  # create new variables "protocol", "agegroup", "intervention", "sacrificetime", "sex" with readable strings
  for (var in c('key.protocol','key.agegroup','key.intervention','key.sacrificetime','registration.sex')){
    d = cols[Field.Name == gsub('.*\.','',var)]
    keys=unname(unlist(strsplit(d[,Categorical.Values],'\|')))
    values=tolower(unname(unlist(strsplit(d[,Categorical.Definitions],'\|'))))
    names(values) = keys
    # match keys to values; create new column
    new_var = gsub(".*\.","",var)
    dmaqc_metadata[,(new_var) := unname(values)[match(get(var), names(values))]]
  }
  dmaqc_metadata[,time_to_freeze := calculated.variables.frozetime_after_acute - calculated.variables.deathtime_after_acute]
  
  # clean up "sacrificetime"
  dmaqc_metadata[,sacrificetime := sapply(sacrificetime, function(x) gsub(' hour.*','h',x))]
  dmaqc_metadata[grepl('immediate', sacrificetime), sacrificetime := 'IPEh']
  
  # clean up 'intervention'
  dmaqc_metadata[grepl('exercise',intervention), intervention := 'acute']
  
  # make "group" "acute_0h"
  dmaqc_metadata[,group := paste0(intervention, '_', sacrificetime)]
  dmaqc_metadata[group=="acute_0h", group := "acute_IPEh"]
  dmaqc_metadata[group=="control_IPEh", group := "control_0h"]
  
  # make a few things strings
  dmaqc_metadata[,specimen.processing.techid := paste0('tech',specimen.processing.techid)]
  dmaqc_metadata[,pid := paste0("pid_", pid)]
  dmaqc_metadata[,bid := paste0("bid_", bid)]
  
  # make viallabel char
  dmaqc_metadata[,viallabel := as.character(viallabel)]
  
  return(dmaqc_metadata)
}

tmm_norm = function(filt_counts){
  raw_dge = edgeR::DGEList(counts=filt_counts)
  dge = edgeR::calcNormFactors(raw_dge, method='TMM')
  tmm = edgeR::cpm(dge,log=TRUE)
  return(as.data.frame(tmm, check.names=F))
}

pass1a_1c_pca = function(counts, meta, tissue, tmm=NULL, repeated_groups = c('acute_0.5h','acute_1h','control_0h'), return_prcomp=F){
  curr_counts = counts
  curr_meta = meta
  
  # normalize counts
  if(is.null(tmm)){
    tmm = tmm_norm(curr_counts)
  }
  
  # run PCA
  curr_pca = prcomp(t(tmm),scale. = T,center = T)
  curr_pcax = curr_pca$x[,1:5]
  if(return_prcomp){
    return(curr_pca)
  }
  explained_var = summary(curr_pca)[["importance"]][3,5]
  # plot
  df = data.frame(
    PC1 = curr_pcax[,1],
    PC2 = curr_pcax[,2],
    randgroup = data_list[[tissue]]$meta[,group],
    sex = data_list[[tissue]]$meta[,sex],
    phase = data_list[[tissue]]$meta[,phase],
    stringsAsFactors = F
  )
  
  # put overlapping groups on the top
  df$randgroup = factor(df$randgroup, levels=c('acute_0.5h','acute_1h','control_0h','acute_4h','acute_7h','acute_24h','acute_48h','control_IPEh','control_0.5h','control_4h','control_7h','acute_IPEh'))
  df = df[rev(order(df$randgroup)),]
  p = ggplot(df) +
    geom_point(aes(x=PC1, y=PC2, fill=randgroup, shape=phase, colour=sex), size=3.5, stroke=1.5) +
    geom_point(data=df[df$randgroup %in% repeated_groups, ], size=0.2, aes(x=PC1, y=PC2)) +
    labs(title=tissue,
         x=sprintf("PC1 (%s%%)", round(summary(curr_pca)[["importance"]][2,1]*100, 2)),
         y=sprintf("PC2 (%s%%)", round(summary(curr_pca)[["importance"]][2,2]*100, 2))) +
    theme_classic() +
    guides(fill=guide_legend(ncol=2, override.aes = list(shape=21))) +
    scale_shape_manual(values=c(pass1a=21, pass1c=24)) +
    scale_colour_manual(values=sex_cols) +
    scale_fill_manual(values=group_cols,
                      limits=c('acute_IPEh','acute_0.5h','acute_1h','acute_4h','acute_7h','acute_24h','acute_48h','control_IPEh','control_0h','control_0.5h','control_4h','control_7h'))
  return(p)
}

# added functionality
myRemoveBatchEffect = function(x,batch=NULL,batch2=NULL,covariates=NULL,design=matrix(1,ncol(x),1),..., returnBeta=F)
  #  Remove batch effects from matrix of expression data
  #  Gordon Smyth and Carolyn de Graaf
  #  Created 1 Aug 2008. Last revised 1 June 2014.
{
  if(is.null(batch) && is.null(batch2) && is.null(covariates)) return(as.matrix(x))
  if(!is.null(batch)) {
    batch <- as.factor(batch)
    contrasts(batch) <- contr.sum(levels(batch))
    batch <- model.matrix(~batch)[,-1,drop=FALSE]
  }
  if(!is.null(batch2)) {
    batch2 <- as.factor(batch2)
    contrasts(batch2) <- contr.sum(levels(batch2))
    batch2 <- model.matrix(~batch2)[,-1,drop=FALSE]
  }
  if(!is.null(covariates)) covariates <- as.matrix(covariates)
  X.batch <- cbind(batch,batch2,covariates)
  fit <- lmFit(x,cbind(design,X.batch),...)
  beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
  beta[is.na(beta)] <- 0
  
  if(returnBeta){
    return(beta)
  }else{
    return(as.matrix(x) - beta %*% t(X.batch))
  }
}

# PASS1A
# requires gsutil access to the internal releases
pass1a_meta = dl_format_pass1a_pheno(tmp, gsutil)
pass1a_meta = pass1a_meta[,.(viallabel, bid, pid, sex, group)]
# select just Sinai PASS1A
pass1a_qc = dl_read_gcp("gs://motrpac-internal-release2-results/pass1a-06/transcriptomics/qa-qc/motrpac_pass1a-06_transcript-rna-seq_qa-qc-metrics.csv", sep=",")

# PASS1C
# just use shipment metadata. DMAQC pheno not yet available
# this actually includes metadaat from all PASS
pass1c_meta = fread("../merged-meta-ref-standards.csv")
pass1c_meta = pass1c_meta[Protocol == 'phase1c', .(vialLabel, bid, pid, sex, intervention, sacrificeTime)]
pass1c_meta[,group := paste0(intervention, '_', sacrificeTime)]
pass1c_meta[,c('intervention', 'sacrificeTime') := NULL]

# merge
colnames(pass1c_meta) = tolower(colnames(pass1c_meta))
pass1c_meta[,viallabel := as.character(viallabel)]
pass1a_meta[,phase := 'pass1a']
pass1c_meta[,phase := 'pass1c']
meta = rbindlist(list(pass1a_meta, pass1c_meta), use.names=T)

# standard tissue
meta[,bic_tissue_code := sapply(viallabel, function(x){
  return(paste0(c('T',unname(unlist(strsplit(x,'')))[8:9]),collapse=''))
})]
meta[,tissue_name_release := bic_animal_tissue_code$tissue_name_release[match(meta[,bic_tissue_code],
                                                                              bic_animal_tissue_code$bic_tissue_code)]]
merge_counts = function(tissue_name_release, # e.g. "t55-gastrocnemius"
                        pass1a_dir = '/projects/motrpac/PASS1A-1C-06/RNA/PASS1A',
                        pass1c_dir = '/projects/motrpac/PASS1A-1C-06/RNA/PASS1C'){
  
  pass1a_file = sprintf("%s/motrpac_pass1a-06_t%s_transcript-rna-seq_rsem-genes-count.txt",
                        pass1a_dir,
                        gsub('^t','',tissue_name_release))
  
  pass1c_file = sprintf("%s/motrpac_20210914_pass1c-06_T%s_rsem_genes_count.txt",
                        pass1c_dir,
                        gsub('^t','',tissue_name_release))
  
  if(!file.exists(pass1a_file) | !file.exists(pass1c_file)){
    message(sprintf("Gene expression for %s not yet available from both phases. Skipping.", tissue_name_release))
    return()
  }
  
  pass1a = fread(pass1a_file, header=T)
  # gastroc RNA-seq data from both sites for PASS1A
  # choose Stanford data
  if(tissue_name_release == "t55-gastrocnemius"){
    if(!exists("pass1a_qc")){
      pass1a_qc = dl_read_gcp("gs://motrpac-internal-release2-results/pass1a-06/transcriptomics/qa-qc/motrpac_pass1a-06_transcript-rna-seq_qa-qc-metrics.csv", sep=",")
    }
    # select Stanford data
    stanford_vl = pass1a_qc[GET_site=='Stanford', vial_label]
    cols = c('gene_id',colnames(pass1a)[colnames(pass1a) %in% as.character(stanford_vl)])
    pass1a = pass1a[, cols, with=F]
  }
  pass1c = fread(pass1c_file, header=T)
  counts = merge(pass1a, pass1c, by='gene_id')
  
  # convert to df
  counts = as.data.frame(counts, check.names=F)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  
  # filter
  min_cpm = 0.5
  min_num_samples = 2
  raw_dge = edgeR::DGEList(counts=counts)
  keep = rowSums(cpm(raw_dge) > min_cpm) >= min_num_samples
  filt_counts = counts[keep,]
  
  return(filt_counts)
}

data_list = list()
for(tissue in unique(meta[,tissue_name_release])){
  counts = merge_counts(tissue)
  if(!is.null(counts)){
    data_list[[tissue]] = list()
    # get meta
    curr_meta = meta[viallabel%in%colnames(counts)]
    data_list[[tissue]]$meta = curr_meta
    # get counts
    counts = counts[curr_meta[,viallabel]]
    data_list[[tissue]]$counts = counts
    # get tmm
    data_list[[tissue]]$norm = tmm_norm(counts)
  }
}
save(data_list, file="data_list.RData")
load("data_list.RData")

for(tissue in names(data_list)){
  print(pass1a_1c_pca(data_list[[tissue]]$counts, data_list[[tissue]]$meta, tissue, tmm=data_list[[tissue]]$norm))
}

nde_list = list()
for (tissue in names(data_list)){
  # initialize table to keep track of N DEGs
  nde = data.table(group=c('all',repeated_groups),
                   n_up_pass1a=NA_real_,
                   n_down_pass1a=NA_real_, 
                   tissue=tissue)
  
  # subset to overlapping groups
  curr_counts = data_list[[tissue]]$counts
  curr_meta = data_list[[tissue]]$meta
  #table(curr_meta[,phase], curr_meta[,group])
  curr_meta = curr_meta[group %in% repeated_groups]
  curr_counts = curr_counts[curr_meta[,viallabel]]
  
  if(tissue %in% c("t63-testes","t64-ovaries")){
    deseq_res = run_deseq(curr_counts, curr_meta, c('group'), 'phase', list(c('phase','pass1a','pass1c')), shrink=F, verbose=T)
  }else{
    deseq_res = run_deseq(curr_counts, curr_meta, c('sex','group'), 'phase', list(c('phase','pass1a','pass1c')), shrink=F, verbose=T)
  }
  lfc = deseq_res$res
  # 5% fdr
  lfc[, adj_p_value := p.adjust(pvalue, method='BH')]
  lfc[, colour := ifelse(adj_p_value < 0.05, tissue_cols[[tissue]], 'black') ]
  # volcano plot
  g = ggplot(lfc, aes(x=log2FoldChange, y=-log10(adj_p_value), colour=colour)) +
    geom_point() +
    theme_classic() +
    labs(y='p-value (-log10)', x='log2 fold change', title=sprintf('%s PASS1A versus PASS1C', tissue),
         subtitle = sprintf("N up: %s (%s%%); N down: %s (%s%%)",
                            nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0]),
                            round(nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0])/nrow(lfc)*100, 2),
                            nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0]),
                            round(nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0])/nrow(lfc)*100,), 2)) +
    geom_vline(xintercept = 0, linetype='dashed') +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') +
    scale_colour_identity(guide='none')
  print(g)
  
  nde[group == 'all', n_up_pass1a := nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0])]
  nde[group == 'all', n_down_pass1a := nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0])]
  
  # look at individual time points
  for (GROUP in repeated_groups ){
    curr_meta2 = curr_meta[group == GROUP]
    curr_counts2 = curr_counts[curr_meta[,viallabel]]
    if(tissue %in% c("t63-testes","t64-ovaries")){
      deseq_res = run_deseq(curr_counts2, curr_meta2, c(), 'phase', list(c('phase','pass1a','pass1c')), shrink=F, verbose=T)
    }else{
      deseq_res = run_deseq(curr_counts2, curr_meta2, c('sex'), 'phase', list(c('phase','pass1a','pass1c')), shrink=F, verbose=T)
    }
    lfc = deseq_res$res
    # 5% fdr
    lfc[,adj_p_value := p.adjust(pvalue, method='BH')]
    lfc[, colour := ifelse(adj_p_value < 0.05, tissue_cols[[tissue]], 'black') ]
    # volcano plot
    g = ggplot(lfc, aes(x=log2FoldChange, y=-log10(adj_p_value), colour=colour)) +
      geom_point() +
      theme_classic() +
      labs(y='p-value (-log10)', x='log2 fold change', title=sprintf('%s PASS1A versus PASS1C: %s', tissue, GROUP),
           subtitle = sprintf("N up: %s (%s%%); N down: %s (%s%%)",
                              nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0]),
                              round(nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0])/nrow(lfc)*100, 2),
                              nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0]),
                              round(nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0])/nrow(lfc)*100,), 2)) +
      geom_vline(xintercept = 0, linetype='dashed') +
      geom_hline(yintercept = -log10(0.05), linetype='dashed') +
      scale_colour_identity(guide='none')
    print(g)
    
    nde[group == GROUP, n_up_pass1a := nrow(lfc[adj_p_value < 0.05 & log2FoldChange > 0])]
    nde[group == GROUP, n_down_pass1a := nrow(lfc[adj_p_value < 0.05 & log2FoldChange < 0])]
  }
  nde_list[[tissue]] = nde
}

nde = rbindlist(nde_list)
nde_melt = melt(nde, id.vars=c('tissue','group'), measure.vars = c('n_up_pass1a','n_down_pass1a'))
nde_melt[variable == 'n_down_pass1a', value := -value]
nde_melt[,group := factor(group, levels=c('all',repeated_groups))]
ggplot(nde_melt, aes(x=tissue, y=value, fill=group)) +
  geom_bar(stat='identity', position='dodge') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, colour='black')) +
  geom_hline(yintercept=0) +
  labs(title="Number of genes differentially expressed between repeated groups (replicates)", y="N DEGs (- = down-reg; + = up-reg)", subtitle = "PASS1A versus PASS1C")


for (tissue in names(data_list)){
  
  tmm = data_list[[tissue]]$norm 
  meta = data_list[[tissue]]$meta
  
  # first, for each group
  betas = list()
  for(GROUP in repeated_groups){
    # subset to this group
    curr_m = meta[group==GROUP]
    curr_tmm = tmm[curr_m[,viallabel]]
    
    # get the batch factors using the same method as limma::removeBatchEffect
    beta = myRemoveBatchEffect(curr_tmm, curr_m[,phase],
                               design = model.matrix(~sex, data = curr_m),
                               returnBeta = T)
    colnames(beta) = GROUP
    betas[[GROUP]] = beta
  }
  
  # use all repeated groups 
  curr_m = meta[group %in% repeated_groups]
  curr_tmm = tmm[curr_m[,viallabel]]
  beta = myRemoveBatchEffect(curr_tmm, curr_m[,phase],
                             design = model.matrix(~sex + group, data = curr_m),
                             returnBeta = T)
  colnames(beta) = 'all_repeated'
  betas[['all_repeated']] = beta
  
  # this is identical to using only repeated groups 
  #   # use all samples 
  #   beta = myRemoveBatchEffect(data_list[[tissue]]$norm, data_list[[tissue]]$meta[,phase],
  #                                 design = model.matrix(~sex + group, data = data_list[[tissue]]$meta),
  #                                 returnBeta = T)
  #   colnames(beta) = 'all_samples'
  #   betas[['all_samples']] = beta
  
  betas_df = data.frame(do.call("cbind", betas))
  betas_df$gene_id = rownames(tmm)
  
  data_list[[tissue]][['batch_betas']] = betas
  
  # plot 'all' versus the rest
  beta_melt = reshape2::melt(betas_df, id.vars='gene_id', measure.vars=c(repeated_groups))
  beta_melt = merge(beta_melt, betas_df[,c('gene_id','all_repeated')], by='gene_id')
  
  min_x = betas_df[which.min(betas$all_repeated),'all_repeated']
  g = ggplot(beta_melt, aes(x=all_repeated, y=value, color=variable)) +
    geom_vline(xintercept = min_x, colour='gray') +
    geom_point(alpha=0.6, size=0.8) +
    theme_classic() +
    scale_colour_manual(values=c(group_cols)) +
    labs(title=sprintf("Estimated gene-level batch effects in %s", tissue), x='Batch effect estimated from all repeated groups',
         y='Batch effect from individual repeated groups') +
    geom_abline() +
    geom_point(data=beta_melt[beta_melt$all_repeated==min_x,], alpha=1, size=1.5)
  print(g)
  
}

for(tissue in names(data_list)){
  tmm = data_list[[tissue]]$norm
  meta = data_list[[tissue]]$meta
  
  assert(all(colnames(tmm) == meta[,viallabel]))
  
  # pca before
  g1 = pass1a_1c_pca(data_list[[tissue]]$counts, data_list[[tissue]]$meta, tissue, tmm=data_list[[tissue]]$norm)
  g1 = g1 + labs(subtitle = "Before batch correction") + theme(legend.position = 'none')
  
  tmm_rbe = myRemoveBatchEffect(tmm, meta[,phase],
                                design = model.matrix(~sex + group, data = meta))
  
  data_list[[tissue]]$norm_corrected = tmm_rbe
  
  # pca after
  g2 = pass1a_1c_pca(data_list[[tissue]]$counts, data_list[[tissue]]$meta, tissue, tmm=tmm_rbe)
  g2 = g2 + labs(subtitle = "After batch correction")
  
  grid.arrange(g1, g2, 
               layout_matrix=rbind(c(1,1,2,2,2),c(1,1,2,2,2)))
}

save(data_list, file="pass1a-1c-merge_data_list.RData")


pc_var_exp = function(prcomp_obj, meta, variable, tissue){
  
  if(!variable %in% colnames(meta)){
    error(sprintf("%s not in the columns of meta", variable))
  }
  
  # get PCs
  meta = data.table(cbind(meta, prcomp_obj$x))
  
  r2 <- c()
  for(pc in paste0('PC',1:10)){
    r2 <- c(r2, summary(lm(get(pc)~get(variable), meta))$r.squared)
  }
  
  df <- data.frame(pc=1:10, value=r2, x=1)
  
  cor <- ggplot(df, aes(x=x, y=rev(pc))) + 
    geom_tile(aes(fill = value), colour = "black") + 
    scale_fill_gradient(low = "white", high = tissue_cols[[tissue]],limits=c(0,1)) +
    theme_grey() + 
    theme(axis.title.y = element_blank(),
          legend.position='none',
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label=paste0(round(value,digits=2)*100,'%')), hjust=0.5) +
    labs(x=variable)
  
  return(cor) # correlation plot 
  
}

# get variance explained from first 10 PCs
var_explained <- function(prcomp_obj, tissue){
  
  eigs = prcomp_obj$sdev^2
  proportion = eigs/sum(eigs)
  
  df = data.frame(pc=1:10, value=proportion[1:10], x=1)
  
  var <- ggplot(df, aes(x=x, y=rev(pc))) + 
    geom_tile(aes(fill = value), colour = "black") + 
    scale_fill_gradient(low = "white", high = tissue_cols[[tissue]],limits=c(0,1)) +
    theme_grey() + 
    theme(axis.title.y = element_blank(),
          legend.position='none',
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label=paste0(round(value,digits=2)*100,'%')), hjust=0.5) +
    labs(x='%VE by PC')
  
  return(var)
}

pc_summary = function(prcomp_obj, meta, tissue, title="%%variance explained in PCs by phase, group, and sex"){
  label_df = data.frame(x=1,y=rev(1:10),label=paste0('PC',1:10))
  labels = ggplot(label_df, aes(x=x, y=y)) + 
    geom_tile(fill='white',colour='black') + 
    theme(axis.title.y = element_blank(),
          legend.position='none',
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label=label), hjust=0.5) +
    labs(x='')
  
  p0 = var_explained(prcomp_obj, tissue)
  p1 = pc_var_exp(prcomp_obj, meta, 'phase', tissue)
  p2 = pc_var_exp(prcomp_obj, meta, 'sex', tissue)
  p3 = pc_var_exp(prcomp_obj, meta, 'exercise', tissue)
  
  grid.arrange(labels, p0, p1, p2, p3, layout_matrix = rbind(c(1,2,2,3,3,4,4,5,5)),
               top = textGrob(title,gp=gpar(fontsize=12)))
}

for(tissue in names(data_list)){
  meta = data_list[[tissue]]$meta
  # make a new "exercise"
  meta[,exercise := NA_character_]
  meta[grepl("control", group), exercise := "control"]
  meta[group %in% c("acute_IPEh","acute_0.5h","acute_1h"), exercise := "early"]
  meta[group %in% c("acute_4h","acute_7h"), exercise := "mid"]
  meta[group %in% c("acute_24h","acute_48h"), exercise := "late"]
  
  # pca before
  pca_before = pass1a_1c_pca(data_list[[tissue]]$counts, 
                             data_list[[tissue]]$meta, 
                             tissue, 
                             tmm=data_list[[tissue]]$norm, 
                             return_prcomp = T)
  pc_summary(pca_before, meta, tissue, title=sprintf("%s: %%variance explained in PCs by phase, exercise, and sex, uncorrected", tissue))
  
  # pca after
  pca_after = pass1a_1c_pca(data_list[[tissue]]$counts, 
                            data_list[[tissue]]$meta, 
                            tissue, 
                            tmm=data_list[[tissue]]$norm_corrected, 
                            return_prcomp = T)
  pc_summary(pca_after, meta, tissue, title=sprintf("%s: %%variance explained in PCs by phase, exercise, and sex, corrected", tissue))
}
