#!/usr/bin/env python3

from sys import argv
import os

if len(argv) <= 4:
	print("Usage: {0} <start_chr> <stop_chr> <pop_1> [<population 2> .....]".format(argv[0]))
	exit()

start_chr = int(argv[1])	# 1
stop_chr = int(argv[2])		# 22
pops = argv[3:]				# ASW CEU YRI

input_vcf_folder = "data/vcf"
sample_ids_file = "SampleIDs/{0}_IDs.txt".format('_'.join(pops))

for chromosome in range(start_chr, stop_chr+1):

	print("###### Analyzing Chr{0} ######".format(chromosome))

	# This is the shortened symlink
	input_vcf_file="chr{0}.vcf.gz".format(chromosome)

	# The folder and file where the recoded SNPs will go
	output_snps_folder = "pops_data/{0}_Data/Chr{1}".format('_'.join(pops), chromosome)
	if not os.path.exists(output_snps_folder):
		os.makedirs(output_snps_folder)
	output_snps_file = "chr{0}.phase3.{1}.SNPs".format(chromosome, '_'.join(pops))

	# Only run vcftools if it hasn't been run before
	if not os.path.isfile("{}/{}.recode.vcf".format(output_snps_folder, output_snps_file)):
		# If we need the "Info" data from the input vcf in the output vcf file, see https://sourceforge.net/p/vcftools/mailman/message/32285447/
	    # You add --recode-INFO-all to the options below. 
	    # Select samples (from relevant populations) with MAF > 0.05, biallelic, no singletons, thin by 10k
		os.system(\
			'vcftools \
				--gzvcf "{0}/{1}" \
				--keep {2} \
				--min-alleles 2 \
				--max-alleles 2 \
				--non-ref-ac 2 \
				--remove-indels \
				--maf 0.05 \
				--recode \
				--out "{3}/{4}"'.format(input_vcf_folder, input_vcf_file, sample_ids_file, output_snps_folder, output_snps_file) \
				)
	else:
		print("Chr{0} already recoded".format(chromosome))
