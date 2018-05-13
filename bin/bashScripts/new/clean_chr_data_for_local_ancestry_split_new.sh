#!/bin/bash/

train="ASW"
test1="CEU"
test2="YRI"

pops="$train_$test1_$test2"

for (( i=22; i<=22; i++ )); do
	echo "###### Analyzing Chr$i ######"

	chr_folder_loc=/Users/gjohnston9/Documents/popgen/"$pops"_Data/Chr"$i"
	mkdir "$chr_folder_loc"/tmp
	tmp_folder_loc="$chr_folder_loc"/tmp
	cd "$chr_folder_loc"

	#### Select samples (from relevant populations) with MAF > 0.05, biallelic, no singletons, thin by 10k
	# vcftools --vcf "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf --min-alleles 2 --max-alleles 2 --non-ref-ac 2 --remove-indels --maf 0.05 --recode --out "$tmp_folder_loc"/chr"$i"."$pops".SNPs ### this is done in previous script
	sed -i.bak 's/#CHROM/CHROM/g' "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf
	grep '#' "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf > "$tmp_folder_loc"/chr"$i"."$pops".header.txt
	echo "###### Completed filtering ######"

	#### Split phased chromosomes
	# Rscript /Users/gjohnston9/Documents/popgen/split_homologous_chr.R "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf
	python3 /Users/gjohnston9/Documents/popgen/split_homologous_chr.py "$chr_folder_loc"/chr"$i".phase3."$pops".SNPs.recode.vcf
	echo "finished running Python script"
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