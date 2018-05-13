#!/bin/bash
echo "Script ${0} Start Time: $(date)"
# This script is expected to be run from the top "popgen" folder.

# The input data files for the vcftools command can be downloaded from http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
# or ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# @ToDo: Why do we have the variable names as train, test1 and test2. Why not pop1, pop2 and pop3 ?
# @ToDo: Handle the linkage of pop input to the SampleIDs files so it is not always ASW_CEU_YRI
train="ASW"
test1="CEU"
test2="YRI"

pops="${train}_${test1}_${test2}"

vcf_folder="data/vcf"
out_folder_base="output/${pops}_Data"
sample_ids_file="SampleIDs/${pops}_IDs.txt"

# for i in {1..22}
for (( i=22; i<=22; i++ )); do
	echo "###### Analyzing Chr$i ######"

	out_folder_loc=${out_folder_base}/Chr${i}
	if [ ! -d $out_folder_loc ]; then
	    mkdir -p $out_folder_loc
	fi

    input_vcf_filename="ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

    output_snps_filename="chr${i}.phase3.${pops}.SNPs"

    # If we need the "Info" data from the input vcf in the output vcf file, see https://sourceforge.net/p/vcftools/mailman/message/32285447/
    # You add --recode-INFO-all to the options below
    # Select samples (from relevant populations) with MAF > 0.05, biallelic, no singletons, thin by 10k
    vcftools  --gzvcf "${vcf_folder}/${input_vcf_filename}" \
                --keep ${sample_ids_file} \
                --min-alleles 2 --max-alleles 2 --non-ref-ac 2 --remove-indels --maf 0.05 --recode \
                --out "${out_folder_loc}/${output_snps_filename}"
done
echo "Script ${0} End Time: $(date)"
