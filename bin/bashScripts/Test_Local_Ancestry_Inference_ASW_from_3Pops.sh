#!/bin/bash/

for (( i=22; i<=22; i++ )); do
	chr_folder_loc=/Users/gjohnston9/Documents/popgen/ASW_CEU_YRI_Data/Chr"$i"
	bin_folder_loc=/Users/gjohnston9/Documents/popgen/bin
	admix_folder_loc=/Users/gjohnston9/Documents/popgen/data/admixture
	sample_ids_folder_loc=/Users/gjohnston9/Documents/popgen/SampleIDs

	echo "###### Analyzing Chr$i ######"
	Rscript "$bin_folder_loc"/RScripts/Test_Local_Ancestry_Inference_ASW_from_3Pops.R \
	"$chr_folder_loc"/chr"$i".phase3.ASW_CEU_YRI.SNPs.homologous.txt \
	"$admix_folder_loc"/all.chr.ASW_CEU_YRI.SNPs.homologous.25k.2.Q \
	"$chr_folder_loc"/chr"$i".ASW_CEU_YRI.SNPs.homologous.25k.2.Q \
	"$sample_ids_folder_loc"/ASW_Sample_IDs_Haploid.txt \
	"$chr_folder_loc"/chr"$i".Test.LA.ASW_Only.txt \
	5
done