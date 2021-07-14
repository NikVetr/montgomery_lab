##########################
########################## Login to an intercactive node on Hoffman 
##########################
qrsh -l h_data=16G,h_rt=24:00:00

##########################
########################## Set directories
##########################
work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR
sum_stat_dir=$work_dir/data/GWAS_Sumstats
ldsc_seg_dir=$work_dir/results/ldsc_reg
ldsc_soft_dir=/u/home/b/bballiu/project-pajukant/software/ldsc

##########################
########################## Download ldsc
##########################
cd /u/home/b/bballiu/project-pajukant/software
git clone https://github.com/bulik/ldsc.git

##########################
########################## Create an environment with _ldsc_ dependencies and activate ldsc 
##########################
echo "Loading dependencies"
. /u/local/Modules/default/init/modules.sh
module load python/anaconda2

echo "Activating conda environment for ldsc"
# conda env create --file $ldsc_soft_dir/environment.yml
source activate ldsc

##########################
###### Download genotypes from 1000 Genomes data for European individuals and limit analyses to HapMap3 SNPs.
##########################
work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR
ldsc_seg_dir=$work_dir/results/ldsc_reg
cd $ldsc_seg_dir 

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz 
tar -xvzf 1000G_Phase3_plinkfiles.tgz
rm 1000G_Phase3_plinkfiles.tgz

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz 
tar -xvzf hapmap3_snps.tgz 
rm hapmap3_snps.tgz 


##########################
########################## Download baseline model LD scores, regression weights, and allele frequencies from the Finucane et al 2018 paper.
##########################
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz 
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
rm 1000G_Phase3_baseline_ldscores.tgz

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar -xvzf weights_hm3_no_hla.tgz
rm weights_hm3_no_hla.tgz

# Only needed for summary statistics files
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2

##########################
########################## Download and process summary statistics data for `ldsc` regression using `munge_sumstats.py`. 
##########################

##########################
########################## IR traits 
##########################

############ BMI from GIANT and UKBiob meta analysis
awk -v OFD='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14}' $sum_stat_dir/BMI_GIANT_UKBiob_Europeans_hg19.txt > $sum_stat_dir/tmp_BMI_GIANT_UKBiob_Europeans_hg19.txt

python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_BMI_GIANT_UKBiob_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/BMI_GIANT_UKBiob_Europeans_hg19 --a1-inc

############ WHR from GIANT and UKBiob meta analysis
awk -v OFD='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14}' $sum_stat_dir/WHR_GIANT_UKBiob_Europeans_hg19.txt > $sum_stat_dir/tmp_WHR_GIANT_UKBiob_Europeans_hg19.txt

python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_WHR_GIANT_UKBiob_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/WHR_GIANT_UKBiob_Europeans_hg19 --a1-inc

############ T2D 
## from Mahajan et al 2018 
awk -v OFD='\t' '{print $1,$2,$3,$4,$5,$6,$8,$9, $13, $14}' $sum_stat_dir/T2D_Mahajan_Europeans_hg19.txt > $sum_stat_dir/tmp_T2D_Mahajan_Europeans_hg19.txt
#change last column to N with 
sed -i 's/Neff/N/g' $sum_stat_dir/tmp_T2D_Mahajan_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_T2D_Mahajan_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/T2D_Mahajan_Europeans_hg19 --a1-inc

# UKBB - Something seems off, very low heritability. 
cp $sum_stat_dir/T2D_UKBB_Europeans_hg19.tsv $sum_stat_dir/tmp_T2D_UKBB_Europeans_hg19.tsv
awk '{print $1}' $sum_stat_dir/tmp_T2D_UKBB_Europeans_hg19.tsv | awk -F':' -v OFS='\t' '{print $3, $4}' > tmp.alleles
paste $sum_stat_dir/tmp_T2D_UKBB_Europeans_hg19.tsv tmp.alleles > $sum_stat_dir/tmp2_T2D_UKBB_Europeans_hg19.tsv
sed -i 's/nCompleteSamples/N/g' $sum_stat_dir/tmp2_T2D_UKBB_Europeans_hg19.tsv
#change allele column names with vi effect_allele and non_effect_allele
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp2_T2D_UKBB_Europeans_hg19.tsv --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/T2D_UKBB_Europeans_hg19 --a1-inc

############ Fasting insulin 
### From Scott et al. 2012 - No longer used because it only has few SNPs
awk '{print $0, "108000"}' $sum_stat_dir/FastInsu_adjBMI_MAGIC_Scott_et_al_Europeans_hg19.txt > $sum_stat_dir/tmp_FastInsu_adjBMI_MAGIC_Scott_et_al_Europeans_hg19.txt
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_FastInsu_adjBMI_MAGIC_Scott_et_al_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/FastInsu_adjBMI_MAGIC_Scott_et_al_Europeans_hg19 --a1-inc

#### From Manning et al 2012
awk '{print $0, "51750"}' FastInsu_MAGIC_Manning_et_al_Europeans_hg19.txt > tmp_FastInsu_MAGIC_Manning_et_al_Europeans_hg19.txt
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_FastInsu_MAGIC_Manning_et_al_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/FastInsu_MAGIC_Manning_et_al_Europeans_hg19 --a1-inc

############ Fasting glucose 
### From Scott et al. 2012 - No longer used because it only has few SNPs
awk '{print $0, "108000"}' FastGlu_MAGIC_Scott_et_al_Europeans_hg19.txt > tmp_FastGlu_MAGIC_Scott_et_al_Europeans_hg19.txt
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_FastGlu_MAGIC_Scott_et_al_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/FastGlu_MAGIC_Scott_et_al_Europeans_hg19 --a1-inc

### From Manning et al 2012 
awk '{print $0, "58074"}' FastGlu_MAGIC_Manning_et_al_Europeans_hg19.txt > tmp_FastGlu_MAGIC_Manning_et_al_Europeans_hg19.txt
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_FastGlu_MAGIC_Manning_et_al_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/FastGlu_MAGIC_Manning_et_al_Europeans_hg19 --a1-inc
rm $sum_stat_dir/tmp_*

############ Surrogate Insulin Sensitivity from MAGIC 
awk '{print $0, "16000"}' $sum_stat_dir/ISI_MAGIC_AdjAgeSexBMI_Europeans_hg19.txt > $sum_stat_dir/tmp_ISI_MAGIC_AdjAgeSexBMI_Europeans_hg19.txt 
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_ISI_MAGIC_AdjAgeSexBMI_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_MAGIC_AdjAgeSexBMI_Europeans_hg19 --a1-inc

awk '{print $0, "16000"}' $sum_stat_dir/ISI_MAGIC_AdjAgeSex_Europeans_hg19.txt > $sum_stat_dir/tmp_ISI_MAGIC_AdjAgeSex_Europeans_hg19.txt 
# replace header with N using vi
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_ISI_MAGIC_AdjAgeSex_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_MAGIC_AdjAgeSex_Europeans_hg19 --a1-inc

############ Insulin Sensitivity 
## From GENESIS + GUARDIAN
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ISI_adjBMI_GENESIS_GUARDIAN_Mixed_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_adjBMI_GENESIS_GUARDIAN_Mixed_hg19 --a1-inc
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ISI_GENESIS_GUARDIAN_Mixed_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_GENESIS_GUARDIAN_Mixed_hg19 --a1-inc

## Insulin Sensitivity from GENESIS
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ISI_GENESIS_adjAgeSex_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_GENESIS_adjAgeSex_Europeans_hg19 --a1-inc
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ISI_GENESIS_adjAgeSexBMI_Europeans_hg19.txt --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ISI_GENESIS_adjAgeSexBMI_Europeans_hg19 --a1-inc

############ HDL 
# from meta analysis of GLGC and MVP (No longer used because of mixed ancestry)
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/HDL_GLGC_MVP_Mixed_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/HDL_GLGC_MVP_Mixed_hg19 --a1-inc
# from Teslovich_et_al
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/HDL_Teslovich_et_al_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/HDL_Teslovich_et_al_Europeans_hg19 --a1-inc

############ TG 
#from meta analysis of GLGC and MVP (No longer used because of mixed ancestry)
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/TG_GLGC_MVP_Mixed_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/TG_GLGC_MVP_Mixed_hg19 --a1-inc
# from Teslovich_et_al
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Triglycerides_Teslovich_et_al_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Triglycerides_Teslovich_et_al_Europeans_hg19 --a1-inc

############ Extreme height 
wget http://portals.broadinstitute.org/collaboration/giant/images/b/bd/GIANT_EXTREME_HEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt.gz
mv GIANT_EXTREME_HEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt Height_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Height_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Height_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19 --a1-inc

############ Extreme BMI
wget http://portals.broadinstitute.org/collaboration/giant/images/9/98/GIANT_EXTREME_BMI_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt.gz
mv GIANT_EXTREME_BMI_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt BMI_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/BMI_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/BMI_Extreme_GIANT_Stage1_Berndt2013_Europeans_hg19 --a1-inc


##########################
########################## psychiatric traits 
##########################

############ Schizophrenia from PGC
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Schizophrenia_SWGPGC_Europeans_hg18.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Schizophrenia_SWGPGC_Europeans_hg18 --a1-inc

#No longer used because of mixed ancestry
awk '{print $0, "150064"}' $sum_stat_dir/SCZ_SWGPGC_Mixed_hg19.txt > $sum_stat_dir/tmp_SCZ_SWGPGC_Mixed_hg19.txt
# replace header with N using vi
mv $sum_stat_dir/tmp_SCZ_SWGPGC_Mixed_hg19.txt $sum_stat_dir/SCZ_SWGPGC_Mixed_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/SCZ_SWGPGC_Mixed_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/SCZ_SWGPGC_Mixed_hg19 --a1-inc

############  MDD 
# from PGC 
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/MDD_PGC_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/MDD_PGC_Europeans_hg19 --a1-inc
# From PGC and PGC+UKBB - No longer used, smaller heritability estimate
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/MDD_PGC_UKBB_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/MDD_PGC_UKBB_Europeans_hg19 --a1-inc

############ Bipolar from PGC
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Bipolar_PGC_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Bipolar_PGC_Europeans_hg19 --a1-inc

############ Alzheimer's Disease
awk '{print $0, "54162"}' $sum_stat_dir/Alzheimer_IGAP_Europeans_hg19.txt > $sum_stat_dir/tmp_Alzheimer_IGAP_Europeans_hg19.txt
# replace header with N using vi
mv $sum_stat_dir/tmp_Alzheimer_IGAP_Europeans_hg19.txt $sum_stat_dir/Alzheimer_IGAP_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Alzheimer_IGAP_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Alzheimer_IGAP_Europeans_hg19 --a1-inc

############ Anorexia Nervosa
mv Anorexia_PGC.tsv Anorexia_PGC_Europeans_hg19.tsv
awk '{print $0, "72517"}' $sum_stat_dir/Anorexia_PGC_Europeans_hg19.tsv > $sum_stat_dir/tmp_Anorexia_PGC_Europeans_hg19.tsv
# replace header with N using vi
mv $sum_stat_dir/tmp_Anorexia_PGC_Europeans_hg19.tsv $sum_stat_dir/Anorexia_PGC_Europeans_hg19.tsv
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Anorexia_PGC_Europeans_hg19.tsv  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Anorexia_PGC_Europeans_hg19 --a1-inc

############ Autism 
wget https://storage.googleapis.com/broad-alkesgroup-public/sumstats_formatted/PASS_Autism.sumstats
mv PASS_Autism.sumstats Autism_PGC_Europeans_hg19.txt
mv Autism_PGC_Europeans_hg19.txt Autism_PGC_Europeans_hg19.sumstats.gz

##########################
########################## autoimmune traits 
##########################

############ Rheumatoid Arthritis
awk '{print $0, "58284"}' $sum_stat_dir/RA_Okada_Europeans_hg19.txt > $sum_stat_dir/tmp_RA_Okada_Europeans_hg19.txt
# replace header with N using vi
mv $sum_stat_dir/tmp_RA_Okada_Europeans_hg19.txt $sum_stat_dir/RA_Okada_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/RA_Okada_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/RA_Okada_Europeans_hg19 --a1-inc


############ Crohn's disease
awk '{print $0, "20883"}' $sum_stat_dir/Crohns_IIDBGC_Europeans_hg19.txt > $sum_stat_dir/tmp_Crohns_IIDBGC_Europeans_hg19.txt
# replace header with N using vi
mv $sum_stat_dir/tmp_Crohns_IIDBGC_Europeans_hg19.txt $sum_stat_dir/Crohns_IIDBGC_Europeans_hg19.txt
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Crohns_IIDBGC_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Crohns_IIDBGC_Europeans_hg19 --a1-inc

############ Primary Biliary Cirrhosis
wget https://storage.googleapis.com/broad-alkesgroup-public/sumstats_formatted/PASS_Primary_biliary_cirrhosis.sumstats
mv PASS_Primary_biliary_cirrhosis.sumstats Primary_biliary_cirrhosis_Europeans_hg19.txt
mv Primary_biliary_cirrhosis_Europeans_hg19.txt Primary_biliary_cirrhosis_Europeans_hg19.sumstats.gz

# # Inflammatory Bowel Disease
# wget ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz
# # Crohn's disease
# wget ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/cd-meta.txt.gz
# # Ulcerative colitis
# wget ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/ucmeta-sumstats.txt.gz

############ Allergy
# UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz
mv disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ALLERGY_ECZEMA_DIAGNOSED --a1-inc

############ Asthma
# UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/disease_ASTHMA_DIAGNOSED.sumstats.gz
mv disease_ASTHMA_DIAGNOSED.sumstats.gz ASTHMA_DIAGNOSED.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/ASTHMA_DIAGNOSED.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/ASTHMA_DIAGNOSED --a1-inc

##########################
########################## reproductive traits 
##########################

############ Age at Menarche 
# 1000 Genomes meta-analysis from Day et al (Nature Genetics 2017) 
wget https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip
unzip Menarche_1KG_NatGen2017_WebsiteUpload.zip 
rm Menarche_1KG_NatGen2017_WebsiteUpload.zip
mv Menarche_1KG_NatGen2017_WebsiteUpload.txt Menarche_NatGen2017_ReproGen_Europeans_hg19.txt
# In R add N=370000 and save as $sum_stat_dir/tmp_Menarche_NatGen2017_ReproGen_Europeans_hg19.txt
mv $sum_stat_dir/tmp_Menarche_NatGen2017_ReproGen_Europeans_hg19.txt $sum_stat_dir/Menarche_NatGen2017_ReproGen_Europeans_hg19.txt
# This will not work because I need rsids! 
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Menarche_NatGen2017_ReproGen_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Menarche_NatGen2017_ReproGen_Europeans_hg19 --a1-inc


# HapMap 2 GWAS meta-analysis results from Perry et al (Nature 2014) 
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/repro_MENARCHE_AGE.sumstats.gz
wget https://www.reprogen.org/Menarche_Nature2014_GWASMetaResults_17122014.zip
unzip Menarche_Nature2014_GWASMetaResults_17122014.zip 
mv Menarche_Nature2014_GWASMetaResults_17122014.txt Menarche_Nature2014_ReproGen_Europeans_hg19.txt
# In R add N=182416 and save as $sum_stat_dir/tmp_Menarche_Nature2014_ReproGen_Europeans_hg19.txt
mv $sum_stat_dir/tmp_Menarche_Nature2014_ReproGen_Europeans_hg19.txt $sum_stat_dir/tmp_Menarche_Nature2014_ReproGen_Europeans_hg19.txt
# This will not work because I need rsids! 
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/tmp_Menarche_Nature2014_ReproGen_Europeans_hg19.txt  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/tmp_Menarche_Nature2014_ReproGen_Europeans_hg19 --a1-inc

##### UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/repro_MENARCHE_AGE.sumstats.gz
mv repro_MENARCHE_AGE.sumstats.gz Menarche_Age_UKBB_Europeans_hg19.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Menarche_Age_UKBB_Europeans_hg19.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Menarche_Age_UKBB_Europeans_hg19 --a1-inc

############ Age at Menopause 
##### UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/repro_MENOPAUSE_AGE.sumstats.gz
mv repro_MENOPAUSE_AGE.sumstats.gz Menopause_Age_UKBB_Europeans_hg19.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/Menopause_Age_UKBB_Europeans_hg19.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/Menopause_Age_UKBB_Europeans_hg19 --a1-inc

##########################
########################## haemotological traits 
##########################

# Diastolic blood pressure UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/bp_DIASTOLICadjMEDz.sumstats.gz
mv bp_DIASTOLICadjMEDz.sumstats.gz blood_pressure_diastolic_UKBB_Europeans_hg19.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_pressure_diastolic_UKBB_Europeans_hg19.sumstats.gz   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_pressure_diastolic_UKBB_Europeans_hg19  --a1-inc

# Systolic blood pressure UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/bp_SYSTOLICadjMEDz.sumstats.gz
mv bp_SYSTOLICadjMEDz.sumstats.gz blood_pressure_systolic_UKBB_Europeans_hg19.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_pressure_systolic_UKBB_Europeans_hg19.sumstats.gz   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_pressure_systolic_UKBB_Europeans_hg19  --a1-inc

# Blood phenotypes
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_EOSINOPHIL_COUNT.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_EOSINOPHIL_COUNT.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_EOSINOPHIL_COUNT --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_LYMPHOCYTE_COUNT.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_LYMPHOCYTE_COUNT.sumstats.gz  --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_LYMPHOCYTE_COUNT --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_MONOCYTE_COUNT.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_MONOCYTE_COUNT.sumstats.gz   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_MONOCYTE_COUNT  --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_MEAN_PLATELET_VOL.sumstats.gz
gunzip blood_MEAN_PLATELET_VOL.sumstats.gz
# manually cap the p values at 1.0E-300.
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_MEAN_PLATELET_VOL.sumstats
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_MEAN_PLATELET_VOL.sumstats
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_MEAN_PLATELET_VOL.sumstats   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_MEAN_PLATELET_VOL  --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_PLATELET_COUNT.sumstats.gz
gunzip blood_PLATELET_COUNT.sumstats.gz
# manually cap the p values at 1.0E-300.
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_PLATELET_COUNT.sumstats
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_PLATELET_COUNT.sumstats
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_PLATELET_COUNT.sumstats   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_PLATELET_COUNT  --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_PLATELET_DISTRIB_WIDTH.sumstats.gz
# manually cap the p values at 1.0E-300.
gunzip blood_PLATELET_DISTRIB_WIDTH.sumstats.gz
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_PLATELET_DISTRIB_WIDTH.sumstats
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_PLATELET_DISTRIB_WIDTH.sumstats
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_PLATELET_DISTRIB_WIDTH.sumstats   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_PLATELET_DISTRIB_WIDTH  --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_RED_COUNT.sumstats.gz
# manually cap the p values at 1.0E-300.
gunzip blood_RED_COUNT.sumstats.gz
sed -i 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_RED_COUNT.sumstats
sed -i 's/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/' $sum_stat_dir/blood_RED_COUNT.sumstats
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_RED_COUNT.sumstats   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_RED_COUNT  --a1-inc

wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/blood_RBC_DISTRIB_WIDTH.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/blood_RBC_DISTRIB_WIDTH.sumstats.gz   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/blood_RBC_DISTRIB_WIDTH  --a1-inc

##########################
########################## cardiometabolic traits 
##########################
# UKBB from Price group
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/disease_CARDIOVASCULAR.sumstats.gz
mv disease_CARDIOVASCULAR.sumstats.gz CARDIOVASCULAR.sumstats.gz
python $ldsc_soft_dir/munge_sumstats.py --sumstats $sum_stat_dir/CARDIOVASCULAR.sumstats.gz   --merge-alleles $ldsc_seg_dir/w_hm3.snplist.txt --out $sum_stat_dir/CARDIOVASCULAR  --a1-inc

##########################
########################## Create annotation files for DE and expressed genes
##########################
# To compute annotation-specific LD scores, we need an annot file, with extension `.annot` or `.annot.gz`. 
# An annot file typically consists of CHR, BP, SNP, and CM columns, followed by one column per annotation, with the value of the annotation for each SNP (0/1 for binary categories, arbitrary numbers for continuous annotations). 
# To create annotation files, we first extract genes DE in each perturbation and cell at 5% FDR. 
# We also extract all genes tested to create annotation of all genes. 

mkdir $ldsc_seg_dir/DE_by_perturb_genes
cat $work_dir/results/tmp_cells | while read cell ; do
cat $work_dir/results/tmp_perturbations | while read perturb ; do
awk -v cell="$cell" -v perturb="$perturb" '{if( $1==cell && $2 == perturb && $24 == "TRUE") print $25}' $work_dir/results/DE_sum_stat_all_treatments_cells_ensg.tsv > $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.tsv
echo done $cell and $perturb
done
done
 
# set of all genes
awk '{print $1}' $work_dir/results/gene_coordinates_hg19.txt > $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.tsv 


# Then we run `make_annot.py` to create annotation files and `ldsc.py` to compute partitioned LD scores based on these annotations. 
# This takes a long time and it is done in parallel on Hoffman2.
