library(data.table)
rdg_mapping <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
load("~/data/smontgom/meta_analysis_results.RData")

mouse_symb_to_human = tapply(rdg_mapping$HUMAN_ORTHOLOG_SYMBOL,
                             rdg_mapping$RAT_SYMBOL,
                             function(x)x[1])
human_entrez_to_symbol = tapply(rdg_mapping$HUMAN_ORTHOLOG_SYMBOL,
                                as.character(rdg_mapping$HUMAN_ORTHOLOG_NCBI_GENE_ID),
                                function(x)x[1])

meta_d = all_meta_analysis_res$`longterm,muscle`
meta_d_base_models = lapply(
  meta_d,
  function(x)x[[which(grepl("base",names(x)))[1]]]
)
# extract the simple base model, useful for I2 scores
meta_d_base_models_simple = lapply(meta_d,function(x)x[["simple:base_model"]])
# create a summary data frame of the human data
human_muscle_res = sapply(meta_d_base_models,
                          function(x)c(x$coeffs[1:2],x$mod_p,x$het_p,x$I2))
human_muscle_res = t(human_muscle_res)
colnames(human_muscle_res) = c("beta","beta_se","p","het_p","I2")
human_muscle_res = as.data.frame(human_muscle_res)
# extract the I2 scores from the simple models
human_muscle_res$I2 = sapply(meta_d_base_models_simple,
                             function(x){
                               if(is.null(x$I2) || is.na(x$I2)){return(0)}
                               return(x$I2)
                             }
)
# map entrez ids to gene symbols, will be used to merge with the 
# rat results
human_muscle_res$symbol = human_entrez_to_symbol[rownames(human_muscle_res)]
human_muscle_res = human_muscle_res[!is.na(human_muscle_res$symbol),]
rownames(human_muscle_res) = human_muscle_res$symbol


meta_d = all_meta_analysis_res$`longterm,blood`
meta_d_base_models = lapply(
  meta_d,
  function(x)x[[which(grepl("base",names(x)))[1]]]
)
# extract the simple base model, useful for I2 scores
meta_d_base_models_simple = lapply(meta_d,function(x)x[["simple:base_model"]])
# create a summary data frame of the human data
human_blood_res = sapply(meta_d_base_models,
                          function(x)c(x$coeffs[1:2],x$mod_p,x$het_p,x$I2))
human_blood_res = t(human_blood_res)
colnames(human_blood_res) = c("beta","beta_se","p","het_p","I2")
human_blood_res = as.data.frame(human_blood_res)
# extract the I2 scores from the simple models
human_blood_res$I2 = sapply(meta_d_base_models_simple,
                             function(x){
                               if(is.null(x$I2) || is.na(x$I2)){return(0)}
                               return(x$I2)
                             }
)
# map entrez ids to gene symbols, will be used to merge with the 
# rat results
human_blood_res$symbol = human_entrez_to_symbol[rownames(human_blood_res)]
human_blood_res = human_blood_res[!is.na(human_blood_res$symbol),]
rownames(human_blood_res) = human_blood_res$symbol

