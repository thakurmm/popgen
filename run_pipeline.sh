#!/usr/bin/env bash


./create_all_admixed_chromosomes.sh
output_11="pops_data/admixed/CEU_YRI_admixed_10admixed_300pure_HOMOLOGOUS.vcf"
output_12="pops_data/admixed/CEU_YRI_admixed_10admixed_300pure_ALLELE_vcf.txt"
output_13="pops_data/admixed/CEU_YRI_admixed_10admixed_300pure_proportions.txt"

output_14="pops_data/admixture/CEU_YRI_admixed_10admixed_300pure.2.Q"

python bin/pythonScripts/generate_test_set.py
output_21=test_input/CEU_YRI_true_population.csv
output_22=test_input/CEU_YRI_test_SNPs_ALLELE_vcf.txt


python3 bin/pythonScripts/Test_Local_Ancestry_Inference_ASW_from_3Pops.py \
    --reference_filename $output_12 \
    --all_admix_filename $output_14 \
    --chrom_admix_filename $output_14 \
    --num_test 8 \
    --test_filename $output_22
output_31="results/source=CEU_YRI_admixed_10admixed_300pure_ALLELE_vcf.txt_test=CEU_YRI_test_SNPs_ALLELE_vcf.txt.txt"


python3 bin/pythonScripts/evaluate_inferences.py \
    --inferences_filename $output_31 \
    --true_ancestry_filename $output_21 \
    --save_plot_filename results/plot_5chroms_try2point5.png

output_41=results/plot_5chroms_try2point5.png
