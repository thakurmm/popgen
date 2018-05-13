#!/bin/bash/

# for i in {1..22}
for (( i=22; i<=22; i++ )); do
	echo "###### Analyzing Chr$i ######"
	in_folder_loc=/Users/gjohnston9/Documents/popgen
	out_folder_loc=/Users/gjohnston9/Documents/popgen/ASW_CEU_YRI_Data
	mkdir -p "$out_folder_loc"/Chr"$i"/
	# genotype_folder_loc=/home/aberens3/Documents/1000_Genomes_Project/Genotypes
	vcftools --gzvcf "$in_folder_loc"/data/vcf/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep "$in_folder_loc"/SampleIDs/ASW_CEU_YRI_IDs.txt --min-alleles 2 --max-alleles 2 --non-ref-ac 2 --remove-indels --maf 0.05 --recode --out "$out_folder_loc"/Chr"$i"/chr"$i".phase3.ASW_CEU_YRI.SNPs
done
