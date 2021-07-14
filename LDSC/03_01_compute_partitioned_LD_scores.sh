#!/bin/bash

cell=$1
perturb=$2
perturb_annot=$3


ldsc_soft_dir=/u/home/b/bballiu/project-pajukant/software/ldsc
work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR
ldsc_seg_dir=$work_dir/results/ldsc_reg
    
echo "Loading dependencies"
. /u/local/Modules/default/init/modules.sh
module load python/anaconda2
    
echo "Activating conda environment for ldsc"
#conda env create --file $ldsc_soft_dir/environment.yml
source activate ldsc
    

# Computation of paritioned LD scores for DE genes
if ([ $perturb_annot -eq 1 ]); then
for chr in $(seq 1 1 22) ; do
    echo Starting computation of paritioned LD scores for DE genes in $cell, $perturb, and chr$chr
    
    echo "Building annotation"
    python $ldsc_soft_dir/make_annot.py \
    --gene-set-file $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.tsv \
    --gene-coord-file $work_dir/results/gene_coordinates_hg19.txt \
    --windowsize 100000 \
    --bimfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file $ldsc_seg_dir/DE_by_perturb_genes/tmp_DE_geneset_${cell}_${perturb}.chr${chr}.annot
    
    echo "Adding annotation name to file"
    sed -i -- "s/\<ANNOT\>/DE_${cell}_${perturb}/g" $ldsc_seg_dir/DE_by_perturb_genes/tmp_DE_geneset_${cell}_${perturb}.chr${chr}.annot
    
    echo "Adding chr, SNP, CM, and BP info"
    paste $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.with_header.bim $ldsc_seg_dir/DE_by_perturb_genes/tmp_DE_geneset_${cell}_${perturb}.chr${chr}.annot > $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.annot
    rm $ldsc_seg_dir/DE_by_perturb_genes/tmp_DE_geneset_${cell}_${perturb}.chr${chr}.annot
    
    gzip $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.annot
    
     echo "Computing paritioned LD scores, limiting to sites present in baseline partition from Finucane et al 2018"
     python $ldsc_soft_dir/ldsc.py --l2 --ld-wind-cm 1 \
     --bfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
     --annot $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.annot.gz \
     --out $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr} \
     --print-snps $ldsc_seg_dir/1000G_EUR_Phase3_baseline/print_snps.txt #$hm3_dir/hm.${chr}.snp

     echo "Add annotation name to file"
     gunzip $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.l2.ldscore.gz
     sed -i -- "s/\<L2\>/DE_${cell}_${perturb}_L2/g" $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.l2.ldscore
     gzip $ldsc_seg_dir/DE_by_perturb_genes/DE_geneset_${cell}_${perturb}.chr${chr}.l2.ldscore
    
    echo Paritioned LD scores done for DE genes in $perturb, $cell, and chromosome $chr
done
fi

# Computation of paritioned LD scores for expressed genes 
if ([ $perturb_annot -eq 2 ]); then
    chr=$4

    echo Starting computation of paritioned LD scores for all genes in chr$chr
    
    echo "Building annotation"
    python $ldsc_soft_dir/make_annot.py \
    --gene-set-file $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.tsv \
    --gene-coord-file $work_dir/results/gene_coordinates_hg19.txt \
    --windowsize 100000 \
    --bimfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/tmp_all_genes.chr${chr}.annot
    
    echo "Adding annotation name to file"
    sed -i -- "s/\<ANNOT\>/all_genes/g" $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/tmp_all_genes.chr${chr}.annot
    
    echo "Adding chr, SNP, CM, and BP info"
    paste $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.with_header.bim $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/tmp_all_genes.chr${chr}.annot > $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.annot
    rm $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/tmp_all_genes.chr${chr}.annot
    
    gzip $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.annot
  
    echo Starting computation of paritioned LD scores for set of all genes
    echo "Computing paritioned LD scores, limiting to sites present in baseline partition from Finucane et al 2018"
    python $ldsc_soft_dir/ldsc.py --l2 --ld-wind-cm 1 \
    --bfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --annot $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.annot.gz \
    --out $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr} \
    --print-snps $ldsc_seg_dir/1000G_EUR_Phase3_baseline/print_snps.txt
    
    echo "Add annotation name to file"
    gunzip $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.l2.ldscore.gz
    sed -i -- "s/\<L2\>/all_genes_L2/g" $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.l2.ldscore
    gzip $ldsc_seg_dir/1000G_EUR_Phase3_all_genes/all_genes.chr${chr}.l2.ldscore
    
    echo Paritioned LD scores done for all genes in chromosome $chr
fi


# Computation of paritioned LD scores for shared versus perturbation specific DE genes 
if ([ $perturb_annot -eq 3 ]); then
    chr=$4

    echo Computing paritioned LD scores for shared versus perturbation specific DE genes 
    
    awk '{print $1,$2,$3,$4,$5}' $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene.chr${chr}.annot > $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SGBS.chr${chr}.annot
    awk '{print $1,$2,$3,$4,$6}' $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene.chr${chr}.annot > $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SKMC.chr${chr}.annot
    awk '{print $1,$2,$3,$4,$7}' $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene.chr${chr}.annot > $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_HEPG2.chr${chr}.annot


    gzip $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_*.chr${chr}.annot

    python $ldsc_soft_dir/ldsc.py --l2 --ld-wind-cm 1 \
    --bfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --annot $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SGBS.chr${chr}.annot.gz \
    --out $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SGBS.chr${chr} \
    --print-snps $ldsc_seg_dir/1000G_EUR_Phase3_baseline/print_snps.txt
  
      python $ldsc_soft_dir/ldsc.py --l2 --ld-wind-cm 1 \
    --bfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --annot $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SKMC.chr${chr}.annot.gz \
    --out $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_SKMC.chr${chr} \
    --print-snps $ldsc_seg_dir/1000G_EUR_Phase3_baseline/print_snps.txt
    
        python $ldsc_soft_dir/ldsc.py --l2 --ld-wind-cm 1 \
    --bfile $ldsc_seg_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --annot $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_HEPG2.chr${chr}.annot.gz \
    --out $ldsc_seg_dir/Nr_Pert_By_Gene/Nr_Pert_By_Gene_HEPG2.chr${chr} \
    --print-snps $ldsc_seg_dir/1000G_EUR_Phase3_baseline/print_snps.txt
    echo Paritioned LD scores done for all genes in chromosome $chr
fi


