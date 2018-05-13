#!/bin/bash/

train="ASW"
test1="CEU"
test2="YRI"

pops="$train_$test1_$test2"

# for i in {1..22}
for (( i=22; i<=22; i++ )); do
	echo "###### Analyzing Chr$i ######"
	in_folder_loc=/Users/gjohnston9/Documents/popgen
	out_folder_loc=/Users/gjohnston9/Documents/popgen/"$pops"_Data
	mkdir -p "$out_folder_loc"/Chr"$i"/
	vcftools --gzvcf "$in_folder_loc"/data/vcf/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep "$in_folder_loc"/SampleIDs/"$pops"_IDs.txt --min-alleles 2 --max-alleles 2 --non-ref-ac 2 --remove-indels --maf 0.05 --recode --out "$out_folder_loc"/Chr"$i"/chr"$i".phase3."$pops".SNPs
done
