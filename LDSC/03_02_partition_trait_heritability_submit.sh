#!/bin/bash

work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR/

part_method=1

# part_method=1 is conditional but I will have to compute the p-values for the Z-score myself
# part_method=2 is conditional and returns the p-values for the Z-score but does not return the paritioned hsq estimates
# part_method=3 is for partition based on number of perturbatios for each DE gene

if ([ $part_method -eq 1 ]); then
  cat $work_dir/results/tmp_cells | while read cell ; do
    cat $work_dir/results/tmp_perturbations | while read perturb ; do
      qsub -N ${cell}_${perturb} -o ${work_dir}/logfiles/part_her_${cell}_${perturb}.o -e ${work_dir}/logfiles/part_her_${cell}_${perturb}.e -l h_data=16G,h_rt=24:00:00 $work_dir/scripts/partition_trait_heritability.sh $part_method $cell $perturb
    done
  done
fi

if ([ $part_method -eq 2 ]); then
  cat $work_dir/results/tmp_traits | while read trait ; do
    trait_name=$(echo $trait | tr "." " " | awk '{print $1}')
    qsub -N ${trait_name} -o ${work_dir}/logfiles/part_her_${trait_name}.o -e ${work_dir}/logfiles/part_her_${trait_name}.e -l h_data=16G,h_rt=24:00:00 $work_dir/scripts/partition_trait_heritability.sh $part_method $trait 
  done
fi


if ([ $part_method -eq 3 ]); then
  cat $work_dir/results/tmp_cells | while read cell ; do
      qsub -N ${cell} -o ${work_dir}/logfiles/part_her_${cell}.o -e ${work_dir}/logfiles/part_her_${cell}.e -l h_data=16G,h_rt=24:00:00 $work_dir/scripts/partition_trait_heritability.sh $part_method $cell
  done
fi
