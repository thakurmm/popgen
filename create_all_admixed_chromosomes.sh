#!/bin/bash
set -x
# @ToDo This script should be moved to bin/bashScripts/new
# This script should be run from popgen folder

#for i in 100,20 200,40 500,100 100,40 200,80 500,200;
# for i in 500,20 500,40;
# for i in 400,20 400,40;
# for i in 300,40;
# for i in 100,30;
# for i in 200,20;
# for i in 250,20 300,20, 350,20;
# for i in 200,30 300,30, 400,30;
mkdir -p pops_data/admixed > /dev/null 2> /dev/null
mkdir -p pops_data/admixture/bed > /dev/null 2> /dev/null
for i in 10,300;
do
#    IFS=',' read num_admixed num_anchor <<< ${i};
    num_admixed=$(echo $i | cut -f1 -d',')
    num_anchor=$(echo $i | cut -f2 -d',')
	python bin/pythonScripts/create_admixed_chromosomes.py --num_admixed_chromosomes $num_admixed --num_anchor $num_anchor --source_pops CEU YRI;

    admixed_output_homologous_file="pops_data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure_HOMOLOGOUS.vcf"
    admixed_output_allele_file="pops_data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure_ALLELE_vcf.txt"
    admixed_proportions_file="pops_data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure_proportions.txt"
    
    plink_output_file_prefix="pops_data/admixture/bed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure"
    # First we prune the file to eliminate SNPs that are linked
	plink2 --noweb --vcf ${admixed_output_homologous_file}  --indep-pairwise 50 5 0.5 --out tmp_plink1_out
	# Now lets use the prune.in file to generate the bed file
	plink2 --noweb --vcf ${admixed_output_homologous_file}  --extract tmp_plink1_out.prune.in --make-bed --thin 0.9 --out ${plink_output_file_prefix};
	/bin/rm tmp_plink1_out*

	admixture --cv ${plink_output_file_prefix}.bed 2;
	mv CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure.2.* pops_data/admixture	
done

# The python script below plots the proportions file output by create_admxed_chromosomes.py AND the .Q file created by the admixture command above
# We had done this using gnuplot - and thus are currently skipping analyzing this program for now.
# We looked through the proportions and .Q file and determined that the files results are a "decent" match.

#python plot_admix_results.py

