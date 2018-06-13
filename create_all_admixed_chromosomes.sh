#!/bin/bash

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
for i in 10,300;
do
    IFS=',' read num_admixed num_anchor <<< ${i};
	python create_admixed_chromosomes.py --num_admixed_chromosomes $num_admixed --num_anchor $num_anchor --source_pops CEU YRI;
	plink --vcf data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure.vcf --make-bed --out data/admixture/new/bed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure;
	admixture --cv data/admixture/new/bed/CEU_YRI_admixed_${num_admixed}admixed_${num_anchor}pure.bed 2;
done

mv CEU_YRI_admixed_* data/admixture/new/admix;

python plot_admix_results.py

