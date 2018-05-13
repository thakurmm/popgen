#!/bin/bash

# This script must be executed from the popgen folder.

# @ToDo: Why do we have the variable names as train, test1 and test2. Why not pop1, pop2 and pop3 ?
# @ToDo: Handle the linkage of pop input to the SampleIDs files so it is not always ASW_CEU_YRI
train="ASW"
test1="CEU"
test2="YRI"

pops="${train}_${test1}_${test2}"

# sample_ids_file="SampleIDs/${pops}_IDs.txt"
# vcf_folder="data/vcf"
# @ToDo Merge the output/ and pops_data symlink
pops_folder_base="pops_data/${pops}_Data"

# for i in {1..22}
for (( i=22; i<=22; i++ )); do
	echo "###### Analyzing Chr$i ######"
	chr_folder_loc="${pops_folder_base}/Chr${i}"

	tmp_folder=${chr_folder_loc}/tmp
	if [ ! -d $tmp_folder ]; then
	    mkdir -p $tmp_folder
	fi

	# cd $chr_folder_loc

	orig_chr_recoded_vcf="${chr_folder_loc}/chr${i}.phase3.${pops}.SNPs.recode.vcf"
	working_chr_recoded_vcf="${tmp_folder}/chr${i}.phase3.${pops}.SNPs.recode.vcf"
	tmp_chr_header_file="${tmp_folder}/chr${i}.phase3.${pops}.header.txt"

    # @ToDo: What is the need for removing the # from the CHROM entry in the vcf file?
    # This command copies the (vcf file) to (vcf file).bak before removing the "#" from "#CHROM" in (vcf file)
    # @ToDo: Do we need to do the two lines below in the select_populations script?
	sed 's/#CHROM/CHROM/g' $orig_chr_recoded_vcf > $working_chr_recoded_vcf
	grep '#' $working_chr_recoded_vcf > $tmp_chr_header_file

	#### Split phased chromosomes
	# Rscript /Users/gjohnston9/Documents/popgen/split_homologous_chr.R "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf
	python3 bin/pythonScripts/split_homologous_chr.py $working_chr_recoded_vcf

	allele_filename=$(echo $working_chr_recoded_vcf | sed 's?recode.vcf?homologous.txt?g')
	homologous_filename=$(echo $working_chr_recoded_vcf | sed 's?recode.vcf?homologous.vcf?g')
	echo "finished running Python script"
exit
	cat "$tmp_folder_loc"/chr"$i"."$pops".header.txt "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.homologous.vcf > "$tmp_folder_loc"/chr"$i"."$pops".SNPs.homologous.header.vcf
	echo "finished something else"
	sed 's/CHROM/#CHROM/g' "$tmp_folder_loc"/chr"$i"."$pops".SNPs.homologous.header.vcf > "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.header.vcf
	echo "###### Finished splitting phased chromosomes ######"

	vcftools --vcf "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.header.vcf --plink-tped --out "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.header2
	echo "finished vcftools"
	plink --tfile "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.header2 --indep-pairwise 50 5 0.5 --out "$tmp_folder_loc"/chr"$i"."$pops".SNPs.homologous.LD
	echo "finished plink part1"
	plink --tfile "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.header2 --extract "$tmp_folder_loc"/chr"$i"."$pops".SNPs.homologous.LD.prune.in --thin-count 25000 --make-bed --out "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.25k
	echo "###### Selected 25k Independent SNPs ######"

	admixture --cv "$chr_folder_loc"/chr"$i"."$pops".SNPs.homologous.25k.bed 2 | tee log2.out
	echo "###### Completed ADMIXTURE estimates ######"

	rm -rf "$tmp_folder_loc"
done