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


files <- paste0("/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/", 
                trimws(GTEx_eQTLs[sapply(motrpac_gtex_map, function(x) grep(x, GTEx_eQTLs))]))
commands <- paste0("scp nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:", files," /Volumes/SSD500GB/GTEx_Analysis_v8_eQTL_all_associations/")
sink(file = "~/dload_eQTLs.sh")
for(i in commands){
  cat(paste0(i, "\n"))
}
sink()
