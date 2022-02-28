
library(data.table)
library(qvalue)
library(ggplot2)
library(gprofiler2)
library(MotrpacBicQC)

#### Nicole's groundwork ####
# source("/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg/pass1b-06/tools/get_fx.R")
# 
# # meta-analysis
# longterm_muscle = dl_read_gcp("gs://mawg-data/pass1b-06/external-data/amar-human-metaanalysis/longterm_muscle.csv", sep=',')
# longterm_blood = dl_read_gcp("gs://mawg-data/pass1b-06/external-data/amar-human-metaanalysis/longterm_blood.csv", sep=',')
# 
# # David recommended only looking at hits from the base model 
# longterm_muscle = longterm_muscle[Group == "base_model"]
# longterm_blood = longterm_blood[Group == "base_model"]
# 
# # pass1b
# transcript_rna_seq = dl_load_gcp("gs://mawg-data/pass1b-06/transcript-rna-seq/dea/transcript_rna_seq_20210804.RData", 'transcript_rna_seq')
# 
# map = dl_read_gcp("gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt")
# head(map)
# table(unique(longterm_blood[,Entrez]) %in% map[,human_gene_symbol])
# table(unique(longterm_muscle[,Entrez]) %in% map[,human_gene_symbol])
# 
# merged_results = list()
# 
# merge_rat_human = function(pass1b, human_ma, map){
#   #human_ma = human_ma[Group %in% c('b0','base_model')]
#   human = merge(human_ma, map[,.(feature_ID, human_gene_symbol)], by.x='Entrez', by.y='human_gene_symbol')
#   rat_human = merge(pass1b, human, by='feature_ID')
#   rat_human = rat_human[,.(feature_ID, tissue, sex, comparison_group, logFC, p_value, selection_fdr, Entrez, Group, Effect, `p-value (-log10)`)]
#   setnames(rat_human, 
#            c("Entrez","Group","Effect","p-value (-log10)"),
#            c("human_gene_symbol","human_model","human_effect","human_p_value"))
#   rat_human[,human_p_value := 10^(-human_p_value)]
#   return(rat_human)
# }
# 
# # blood
# pass1b_blood = data.table(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue == "t30-blood-rna" ,])
# blood_merged = merge_rat_human(pass1b_blood, longterm_blood, map)
# 
# # vastus 
# pass1b_vl = data.table(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue == "t56-vastus-lateralis" ,])
# vl_merged = merge_rat_human(pass1b_vl, longterm_muscle, map)
# 
# # gastroc
# pass1b_gastroc = data.table(transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue == "t55-gastrocnemius" ,])
# gastroc_merged = merge_rat_human(pass1b_gastroc, longterm_muscle, map)
# 
# merged_results[["t30-blood-rna"]] = blood_merged
# merged_results[["t56-vastus-lateralis"]] = vl_merged
# merged_results[["t55-gastrocnemius"]] = gastroc_merged
# 
# calc_pi1 = function(p_values){
#   out = tryCatch(
#     {
#       pi1 = 1-qvalue(p_values)$pi0
#       return(pi1)
#     },
#     error=function(cond) {
#       pi1 = 1-qvalue(p_values, lambda = seq(0, max(p_values), 0.05))$pi0
#       return(pi1)
#     }
#   )    
#   return(out)
# }
# 
# pi1_list=list()
# for(t in names(merged_results)){
#   dt = merged_results[[t]]
#   pi1s = dt[,list(pi1 = calc_pi1(p_value)), by=.(sex, comparison_group)]
#   pi1s[,tissue := t]
#   pi1_list[[t]] = pi1s
# }
# all_pi1s = rbindlist(pi1_list)
# 
# ggplot(all_pi1s, aes(x=sex, y=pi1, fill=comparison_group)) +
#   geom_bar(stat='identity', position='dodge', colour='black') +
#   theme_classic() + 
#   facet_wrap(~tissue) +
#   scale_fill_manual(values=group_cols) +
#   theme(axis.text = element_text(colour='black'),
#         axis.title.x = element_blank()) +
#   labs(y='Replication (pi1)', title='Replication of human meta-analysis hits in PASS1B')
# 
# merged_results = lapply(merged_results, as.data.frame, stringsAsFactors=F)
# save(merged_results, file="/oak/stanford/groups/smontgom/shared/motrpac/shared_rdata/pass1b-human-metaanalysis__merged_results.RData")
# 

#### load Nicole's table ####
load("~/data/smontgom/pass1b-human-metaanalysis__merged_results.RData")
sign_match <- lapply(setNames(c("male", "female"), c("male", "female")), function(sex){
  sexmatch <- sapply(names(merged_results), function(tissue){
    sapply(paste0(c(2^(0:3)), "w"), function(time) {
      dat <- as.data.frame(merged_results[[tissue]])
      subdat <- dat[dat$comparison_group == time & dat$sex == sex,]
      round(mean(sign(subdat$logFC) == sign(subdat$human_effect)), 3)
      # length(unique(subdat$human_gene_symbol))
    })
  })
  colnames(sexmatch) <- MotrpacBicQC::bic_animal_tissue_code$abbreviation[match(colnames(sexmatch), MotrpacBicQC::bic_animal_tissue_code$tissue_name_release)]
  sexmatch
})

# sign_match <- lapply(setNames(c("male", "female"), c("male", "female")), 
#                               function(sex) round(sign_match[[sex]]*100, 1))

sign_match

#also read in the raw data?
muscle <- read.csv(file = "~/data/smontgom/longterm_muscle.csv")
muscle <- muscle[muscle$Group == "base_model",]
blood <- read.csv(file = "~/data/smontgom/longterm_blood.csv")
blood <- blood[blood$Group == "base_model",]

n_genes_per_tissue <- lapply(setNames(c("male", "female"), c("male", "female")), function(sex){
  sapply(names(merged_results), function(tissue){
    sapply(paste0(c(2^(0:3)), "w"), function(time) {
      dat <- as.data.frame(merged_results[[tissue]])
      subdat <- dat[dat$comparison_group == time & dat$sex == sex,]
      length(unique(subdat$human_gene_symbol))
    })
  })
})

pvalue_coercion_val <- 1E-5

all_log10_pvals <- log10(do.call(rbind, merged_results)$p_value)
all_log10_pvals[all_log10_pvals < log10(pvalue_coercion_val)] <- log10(pvalue_coercion_val)
log10_window_size <- mean(diff(quantile(all_log10_pvals, 1:10/10))) * 5
log_pval_windows <- seq(min(all_log10_pvals), max(all_log10_pvals) - log10_window_size, length.out = 100)
log_pval_windows <- cbind(log_pval_windows, log_pval_windows + log10_window_size)


windowed_props <- lapply(setNames(names(merged_results), names(merged_results)), function(tissue){
  lapply(setNames(c("male", "female"), c("male", "female")), function(sex){
    sapply(paste0(c(2^(0:3)), "w"), function(time) {
      sapply(1:nrow(log_pval_windows), function(w){
        dat <- as.data.frame(merged_results[[tissue]])
        dat$p_value[dat$p_value < log10(pvalue_coercion_val)] <- log10(pvalue_coercion_val) + 1E-15
        subdat <- dat[dat$comparison_group == time & dat$sex == sex & log10(dat$p_value) >= log_pval_windows[w,1] & log10(dat$p_value) <= log_pval_windows[w,2],]
        mean(sign(subdat$logFC) == sign(subdat$human_effect))
      })
    })
  })
})

windowed_props_n <- lapply(setNames(names(merged_results), names(merged_results)), function(tissue){
  lapply(setNames(c("male", "female"), c("male", "female")), function(sex){
    sapply(paste0(c(2^(0:3)), "w"), function(time) {
      sapply(1:nrow(log_pval_windows), function(w){
        dat <- as.data.frame(merged_results[[tissue]])
        dat$p_value[dat$p_value < log10(pvalue_coercion_val)] <- log10(pvalue_coercion_val) + 1E-15
        subdat <- dat[dat$comparison_group == time & dat$sex == sex & log10(dat$p_value) >= log_pval_windows[w,1] & log10(dat$p_value) <= log_pval_windows[w,2],]
        length(subdat$logFC)
      })
    })
  })
})
n_range <- range(unlist(windowed_props_n))
lwds <- log2(n_range[1]:(n_range[2]+1) + 1)
lwds <- lwds / max(lwds) * 3

rep2 <- function(x, n) c(sapply(x, function(i) rep(i, n)))
segments2 <- function(x0, x1, y0, y1, lwd1, lwd2, col = 1, ...){
  lt <- c(lwd1, lwd2) / 96 / par("pin")[2] * (par("usr")[4] - par("usr")[3])
  polygon(x = c(x0, x1, x1, x0), y = c(y0 + lt[1]/2, y1 + lt[2]/2, y1 - lt[2]/2, y0 - lt[1]/2), col = col, border = col)
}

lines2 <- function(x, y, lwds, col = 1){
  for(i in 1:(length(x)-1)){
    segments2(x0 = x[i], x1 = x[i+1], y0 = y[i], y1 = y[i+1],
              lwd1 = lwds[i], lwd2 = lwds[i+1], col = col)
  }
}

#### just the plotting ####

cairo_pdf("~/Documents/Documents - nikolai/pass1b_fig8_rat-man_sign_replication.pdf", width = 800 / 72, height = 500 / 72, family="Arial Unicode MS")
par(mfrow = c(2,3), mar = c(4,4,4,2))
for(sex in c("male", "female")){
  for(tissue in names(merged_results)){
    plot(100, 100, xlim = range(-apply(log_pval_windows, 1, mean)), cex.lab = 1.35, cex.axis = 1.3,
         ylim = c(0,1), xlab = latex2exp::TeX("$-log_{10}(pval)$"), ylab = "proportion of matching signs")
    
    #horizontal line at 0.5
    abline(h = 0.5, lwd = 2, col = "grey50", lty = 2)
    
    #tissue and sex labels
    text(tissue, x = mean(range(-apply(log_pval_windows, 1, mean))), y = 1.05, xpd = NA, 
         col = tissue_cols[tissue], pos = 3, cex = 2, font = 3)
    text(x = max(-apply(log_pval_windows, 1, mean)) - c(male = diff(par("usr")[1:2]) * 0.01, female = 0)[sex], y = 1.04, cex = 2.5, font = 2, pos = 1,
         labels = c(male = "\u2642", female = "\u2640")[sex], col = sex_cols[sex])
    
    
    for(time in paste0(c(2^(0:3)), "w")){
      # lines(x = -rep2(apply(log_pval_windows, 1, mean), 2), col = group_cols[[time]], lwd = 2,
      #       y = rep2(windowed_props[[tissue]][[sex]][,time], 2))
      lines2(x = -apply(log_pval_windows, 1, mean), col = group_cols[[time]],
             lwds = lwds[windowed_props_n[[tissue]][[sex]][,time]+1],
             y = windowed_props[[tissue]][[sex]][,time])
      
    }
    if(sex == "female" & tissue == names(merged_results)[length(merged_results)]){
      legend( x = par("usr")[2] - diff(par("usr")[1:2]) * 0.375, 
              y = par("usr")[3] + diff(par("usr")[3:4]) * 0.4, legend = c("50% match", paste0(c(2^(0:3)), "w")), 
             lwd = rep(2, 5), col = c("grey50", group_cols[paste0(c(2^(0:3)), "w")]), lty = c(2,1,1,1,1))
      
      for(lwdi in 1:(length(lwds) - 1)){
        segments
      }
      
    }
    
    #n genes in each tissue
    text(labels = paste0(n_genes_per_tissue[[sex]]["8w",tissue], " genes being considered"),
         x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4]) * 0.03, pos = 2, xpd = NA)
    
  }
}
dev.off()