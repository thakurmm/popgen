#!/bin/bash

set -x # echo command before running

python Test_Local_Ancestry_Inference_ASW_from_3Pops.py --reference_filename data/admixed/CEU_YRI_admixed_200admixed_40pure_homologous.vcf --all_admix_filename data/admixture/new/admix/CEU_YRI_admixed_200admixed_40pure.2.Q --chrom_admix_filename data/admixture/new/admix/CEU_YRI_admixed_200admixed_40pure.2.Q --num_test 15 --test_filename test_input/CEU_YRI_test_cases_chrom_try2.csv