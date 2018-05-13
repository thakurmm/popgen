#!/bin/bash/

for i in {1..22}
do
folder_loc=/home/aberens3/Documents/LocalAncestryMethod/ASW_CEU_YRI_Data/Chr"$i"
all_folder_loc=/home/aberens3/Documents/LocalAncestryMethod/ASW_CEU_YRI_Data/AllChr
bin_folder_loc=/home/aberens3/Documents/LocalAncestryMethod/ASW_CEU_YRI_Data/bin
echo "###### Analyzing Chr$i ######"
Rscript "$bin_folder_loc"/Test_LAI_ASW_from_3Pops_SlidingWindows.R \
"$folder_loc"/chr"$i".phase3.ASW_CEU_YRI.SNPs.homologous.txt \
"$all_folder_loc"/all.chr.ASW_CEU_YRI.SNPs.homologous.25k.2.Q \
"$folder_loc"/chr"$i".ASW_CEU_YRI.SNPs.homologous.25k.2.Q \
"$bin_folder_loc"/ASW_Sample_IDs_Haploid.txt \
"$folder_loc"/chr"$i".Test.LAI.ASW_Only.SlidingWindows.txt \
5
done
