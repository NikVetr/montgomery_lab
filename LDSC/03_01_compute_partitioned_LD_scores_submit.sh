#!/bin/bash

work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR

perturb_annot=3

# Computation of paritioned LD scores for DE genes
if ([ $perturb_annot -eq 1 ]); then
    cat $work_dir/results/tmp_cells | while read cell ; do
        cat $work_dir/results/tmp_perturbations | while read perturb ; do
            qsub -N ${cell}_${perturb} -o ${work_dir}/logfiles/ldsc_${cell}_${perturb}.o -e ${work_dir}/logfiles/ldsc_${cell}_${perturb}.e -l h_data=16G,h_rt=24:00:00 compute_partitioned_LD_scores.sh $cell $perturb $perturb_annot
        done
    done
fi

# Computation of paritioned LD scores for expressed genes 
if ([ $perturb_annot -eq 2 ]); then
    cell=NULL
    perturb=NULL
    for chr in $(seq 1 1 22) ; do
      qsub -N all_genes_${chr} -o ${work_dir}/logfiles/ldsc_all_genes_${chr}.o -e ${work_dir}/logfiles/ldsc_all_genes_${chr}.e -l h_data=16G,h_rt=24:00:00 $work_dir/scripts/compute_partitioned_LD_scores.sh $cell $perturb $perturb_annot $chr
    done
fi

# Computation of paritioned LD scores for shared versus perturbation specific DE genes 
if ([ $perturb_annot -eq 3 ]); then
    cell=NULL
    perturb=NULL
    for chr in $(seq 1 1 22) ; do
      qsub -N all_genes_${chr} -o ${work_dir}/logfiles/ldsc_cont_annot_${chr}.o -e ${work_dir}/logfiles/ldsc_cont_annot_${chr}.e -l h_data=16G,h_rt=24:00:00 $work_dir/scripts/compute_partitioned_LD_scores.sh $cell $perturb $perturb_annot $chr
    done
fi
