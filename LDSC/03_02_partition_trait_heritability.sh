#!/bin/bash

part_method=$1

ldsc_soft_dir=/u/home/b/bballiu/project-pajukant/software/ldsc
work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR
ldsc_seg_dir=$work_dir/results/ldsc_reg
sum_stat_dir=$work_dir/data/GWAS_Sumstats

echo "Loading dependencies"
. /u/local/Modules/default/init/modules.sh
module load python/anaconda2

echo "Activating conda environment for ldsc"
source activate ldsc

if ([ $part_method -eq 1 ]); then
  cell=$2
  perturb=$3
  
  cat $work_dir/results/tmp_traits | while read trait ; do
    trait_name=$(echo $trait | tr "." " " | awk '{print $1}')
    echo Computing % of $trait_name heritability explained by genes DE in $cell and $perturb vs baseline and all genes

    # This is conditional but I will have to compute the p-values for the Z-score myself
    python $ldsc_soft_dir/ldsc.py --h2 ${sum_stat_dir}/${trait} \
    --w-ld-chr $ldsc_seg_dir/weights_hm3_no_hla/weights. \
    --ref-ld-chr $ldsc_seg_dir/1000G_EUR_Phase3_baseline/baseline.,$ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr,$ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr \
    --overlap-annot \
    --frqfile-chr $ldsc_seg_dir/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients \
    --out $ldsc_seg_dir/hsq_enrich_res/${trait_name}_${cell}_${perturb}

    echo Completed computation of % of $trait heritability explained by genes DE in $cell and $perturb
  done

fi

if ([ $part_method -eq 2 ]); then

  # This directly / only gives the p-values for the Z-score 
   echo Partitioning heritability of $trait and testing for enrichment of each perturbation-related annotation vs annotation of all genes and baseline annotation
  trait=$2

  python $ldsc_soft_dir/ldsc.py --h2-cts ${sum_stat_dir}/${trait} \
  --w-ld-chr $ldsc_seg_dir/weights_hm3_no_hla/weights. \
  --ref-ld-chr $ldsc_seg_dir/1000G_EUR_Phase3_baseline/baseline. \
  --ref-ld-chr-cts $ldsc_seg_dir/test_DEbyPerturb_vs_all_genes.txt \
  --overlap-annot \
  --frqfile-chr $ldsc_seg_dir/1000G_Phase3_frq/1000G.EUR.QC. \
  --print-coefficients \
  --out $ldsc_seg_dir/hsq_enrich_res/${trait}
  
  echo Completed computation of % of $trait heritability explained by each perturbation-related annotation vs annotation of all genes and baseline annotation
fi


if ([ $part_method -eq 3 ]); then
  cell=$2
  
  cat $work_dir/results/tmp_traits | while read trait ; do
    trait_name=$(echo $trait | tr "." " " | awk '{print $1}')
    echo Computing % of $trait_name heritability explained by genes DE in $cell conditional on baseline annotations and annotations of all genes

    # This is conditional but I will have to compute the p-values for the Z-score myself
    python $ldsc_soft_dir/ldsc.py --h2 ${sum_stat_dir}/${trait} \
    --w-ld-chr $ldsc_seg_dir/weights_hm3_no_hla/weights. \
    --ref-ld-chr $ldsc_seg_dir/1000G_EUR_Phase3_baseline/baseline.,$ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr,$ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_${cell}.chr \
    --overlap-annot \
    --frqfile-chr $ldsc_seg_dir/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients \
    --out $ldsc_seg_dir/hsq_enrich_res/${trait_name}_${cell}_sharedVSspecific

    echo Completed computation of % of $trait heritability explained by shared versus perturbation-specific DE genes in $cell
  done

fi