library(ks)
library(data.table)
library(edgeR)
library(limma)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
library(MotrpacBicQC)
library(plotrix)
library(ggplot2)
library(testit)
library(circlize)
library(jpeg)
library(foreach)
library(doParallel)

####load deg eqtl data ####
# load MotrpacBicQC conversion tables 
bic_animal_tissue_code = as.data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code = bic_animal_tissue_code[tissue_name_release!='']
bic_animal_tissue_code[,my_tissue := tolower(gsub(" Powder", "", bic_tissue_name))]
bic_animal_tissue_code[,my_tissue := gsub(' ','_',my_tissue)]
bic_animal_tissue_code[my_tissue == 'blood_rna', my_tissue := 'paxgene_rna']
tissue_codes = bic_animal_tissue_code[, tissue_name_release]
tissue_codes = tissue_codes[tissue_codes != '']

## Define data paths

#download data to local drive
# system("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/shared/motrpac/mawg_data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt /Users/nikolai/data/smontgom/")

GTEx_logo <- readJPEG("~/Documents/GTEx_logo.jpg")
eqtl = '~/data/smontgom/GTEx_Analysis_v8_eQTL'
#deg = 'gs://mawg-data/pass1b-06/transcript-rna-seq/dea/'
deg = '~/data/smontgom/dea/'
#map = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/mapping/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t')
map = fread('~/data/smontgom/pass1b-06_transcript-rna-seq_feature-mapping_20201002.txt', sep='\t', header=T)
gwas = '~/data/smontgom/imputed_gwas_hg38_1.1'
coloc = '~/data/smontgom/results_enloc_priors'

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



if(!exists("rna_dea")){
  load('~/data/smontgom/rna_dea_20210114.RData')
  load('~/data/smontgom/dea/transcript_rna_seq_20210126.RData')
  old_rna_dea <- rna_dea
  rna_dea$training_dea <- as.data.table(transcript_rna_seq$training_dea)
  rna_dea$timewise_dea <- as.data.table(transcript_rna_seq$timewise_dea)
  new_rna_dea <- rna_dea
}

if(!exists("deg_eqtl_list")){
  deg_eqtl_list = list()
  for(motrpac_tissue in unique(rna_dea$timewise_dea$tissue)){
    
    if(!motrpac_tissue %in% names(motrpac_gtex_map)){next}
    
    cat(paste0(motrpac_tissue, "\n"))
    
    # read in eQTLs
    gtex_tissue = motrpac_gtex_map[[motrpac_tissue]]
    gtex_egene = fread(sprintf('%s/%s.v8.egenes.txt.gz',eqtl,gtex_tissue), sep='\t', header=T)
    gtex_egene[, human_ensembl_gene := gsub('\\..*','',gene_id)]
    
    # match human genes with rat ensembl genes 
    gtex_egene = merge(gtex_egene, map, by='human_ensembl_gene')
    gtex_motrpac = merge(gtex_egene, rna_dea$timewise_dea[tissue == motrpac_tissue], by='feature_ID')
    gtex_motrpac$abs_slope <- abs(gtex_motrpac$slope)
    
    deg_eqtl_list[[motrpac_tissue]] = gtex_motrpac

  }
}

cols = list(Tissue=tissue_cols[names(deg_eqtl_list)], 
            Time=group_cols[paste0(c(1,2,4,8), "w")],
            Sex=sex_cols[c('male','female')])
cols$Tissue[which(is.na(cols$Tissue))] <- '#C0C0C0'
names(cols$Tissue)[which(is.na(names(cols$Tissue)))] <- "t1000-gonads"

