#!/bin/bash

for i in 100,20;
do
    IFS=',' read num_admixed num_pure <<< ${i};
    vcf_filename=data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_pure}pure_homologous.vcf
    ancestry_filename=data/admixed/CEU_YRI_admixed_${num_admixed}admixed_${num_pure}pure_proportions.txt
	python generate_rfmix_input.py $vcf_filename $ancestry_filename
done