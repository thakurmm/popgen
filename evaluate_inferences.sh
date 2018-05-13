#!/bin/bash

set -x # echo command before running

python evaluate_inferences.py --inferences_filename results/source=CEU_YRI_admixed_200admixed_40pure_homologous.vcf_test=CEU_YRI_test_cases_chrom_try2.csv.txt --true_ancestry_filename test_input/CEU_YRI_test_cases_true_ancestry_try2.csv --save_plot_filename results/plot_5chroms_try2point5.png