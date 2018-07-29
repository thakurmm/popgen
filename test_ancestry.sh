#!/bin/bash

set -x # echo command before running

# Reading the python script help, it appears all_admix_filename is meant to be the same as chrom_admix_filename, except
# for ALL chromosomes (chr1, chr2 ..... ). However, since we are only working with Chr22 - the two files are
# the same below.

# num_test applies to number of entries to be used from the "test_filename" below, unless you specify num_test = -1
# which means, use ALL the entries from "test_filename"
# num_test must be <= number of entries in test_filename OR just -1

python bin/pythonScripts/Test_Local_Ancestry_Inference_ASW_from_3Pops.py \
    --reference_filename pops_data/admixed/CEU_YRI_admixed_10admixed_300pure_ALLELE_vcf.txt \
    --all_admix_filename pops_data/admixture/CEU_YRI_admixed_10admixed_300pure.2.Q \
    --chrom_admix_filename pops_data/admixture/CEU_YRI_admixed_10admixed_300pure.2.Q \
    --num_test 8 \
    --test_filename test_input/CEU_YRI_test_SNPs_ALLELE_vcf.txt