##########################
########################## Login to an intercactive node on Hoffman 
##########################
qrsh -l h_data=16G,h_rt=24:00:00
date=20201120

##########################
########################## Set directories
##########################
ldsc_soft_dir=/u/home/b/bballiu/project-pajukant/software/ldsc
work_dir=/u/home/b/bballiu/project-pajukant/MoMeIR
ldsc_seg_dir=$work_dir/results/ldsc_reg
sum_stat_dir=$work_dir/data/IR_GWAS_Sumstats

# See partition_trait_heritability.sh and partition_trait_heritability_submit.sh for how the analysese were performed.

##########################
##########################  Extract total observed heritability for each trait
##########################
grep 'Total Observed scale h2' $ldsc_seg_dir/hsq_enrich_res/*.log > $ldsc_seg_dir/hsq_enrich_res/total_observed_scale_h2_$date.tsv

##########################
########################## Extract proportion of heritability and importance of each annotation  
##########################
head -n 1 $ldsc_seg_dir/hsq_enrich_res/WHR_GIANT_UKBiob_Europeans_hg19_SKMC_DEXA.results | awk -v OFS='\t' '{print "Trait", $0}' > $ldsc_seg_dir/hsq_enrich_res/hsq_partition_by_perturb_$date.tsv

cat $work_dir/results/tmp_cells | while read cell ; do
cat $work_dir/results/tmp_perturbations | while read perturb ; do
cat $work_dir/results/tmp_traits | while read trait ; do

trait_name=$(echo $trait | tr "." " " | awk '{print $1}')

tail -n 1 $ldsc_seg_dir/hsq_enrich_res/${trait_name}_${cell}_${perturb}.results | awk -v OFS='\t' -v var="$trait_name" '{print var, $0}' >> $ldsc_seg_dir/hsq_enrich_res/hsq_partition_by_perturb_$date.tsv
done
done
done


##########################
##########################  Extract importance of each annotation and p-value 
##########################
head -n 1 $ldsc_seg_dir/hsq_enrich_res/BMI_GIANT_UKBiob_Europeans_hg19.sumstats.gz.cell_type_results.txt  | awk -v OFS='\t' '{print "Trait", $0}' > $ldsc_seg_dir/hsq_enrich_res/hsq_partition_by_perturb_enrich_pval_$date.tsv

cat $work_dir/results/tmp_traits | while read trait ; do
trait_name=$(echo $trait | tr "." " " | awk '{print $1}')
tail -n +2 $ldsc_seg_dir/hsq_enrich_res/${trait_name}.sumstats.gz.cell_type_results.txt | awk -v OFS='\t' -v var="$trait_name" '{print var, $0}' >> $ldsc_seg_dir/hsq_enrich_res/hsq_partition_by_perturb_enrich_pval_$date.tsv
done



##########################
########################## Extract coefficients and z-scores of continuous annotation 
##########################
head -n 1 $ldsc_seg_dir/hsq_enrich_res/WHR_GIANT_UKBiob_Europeans_hg19_SKMC_sharedVSspecific.results | awk -v OFS='\t' '{print "Trait.Cell", $0}' > $ldsc_seg_dir/hsq_enrich_res/hsq_partition_sharedVSspecific.tsv

cat $work_dir/results/tmp_traits | while read trait ; do
cat $work_dir/results/tmp_cells | while read cell ; do
trait_name=$(echo $trait | tr "." " " | awk '{print $1}')
tail -n 1 $ldsc_seg_dir/hsq_enrich_res/${trait_name}_${cell}_sharedVSspecific.results  | awk -v OFS='\t' -v var="$trait_name.$cell" '{print var, $0}' >> $ldsc_seg_dir/hsq_enrich_res/hsq_partition_sharedVSspecific.tsv
done
done


scp bballiu@hoffman2.idre.ucla.edu:/u/home/b/bballiu/MoMeIR/results/ldsc_reg/hsq_enrich_res/total_observed_scale_h2* .
scp bballiu@hoffman2.idre.ucla.edu:/u/home/b/bballiu/MoMeIR/results/ldsc_reg/hsq_enrich_res/hsq_partition_by_perturb_* .
