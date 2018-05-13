#!/usr/bin/env python3

### % purple? % correct? time taken?
### ^ calculate this for each set of test chromosomes with {4, 6, 8} switches
####################### effects of modifying recomb rate?
### maybe use this information to change window size (kmer)
### and look for areas on chromosome that are frequently incorrect, or ambiguous

### TO DO:
### 1) exclude test chromosomes from training params for HMM (select_informative_SNPs, creation of emission matrix and transition matrix)
### 2) large, permanent test set?
### 3) implement Viterbi???

### TO ASK ABOUT:
### 1) problem of finding locations with highest error rate
###		- weight by specific window size (i.e. 5 informative SNPs that were far from each other -> large area -> more likely to have some error)
###		- weight by percent of error?? (this makes the most sense but it's the hardest to implement...)
###		- window of windows, find find error rate for each [window of windows], much more variability this way (num_chromosomes * outer_window_size)

### keep track of ambiguous vs incorrect calls at each SNP

### November 2017 changes:
# - training set is obtained by loading .vcf file produced by create_admixed_chromosomes; slightly different file format
#		(first line is a comment, and header has #CHROM instead of CHROM)
# - train_ids are no longer relevant, since each individual in the training file is used

### 2018 TODO:
# - accept data in RFMix format
# - fix HMM warnings
# - do some more profiling
# - rename files to sensible titles

import Chr_kmer_Informative_SNPs_Haploid_HMM as HMM

import numpy as np
import pandas as pd

import argparse
import itertools
import os
import pdb
import random
import time

pd.set_option('display.float_format', lambda x: '%.4f' % x)

# homologous_filename = argv[1]
# all_admixture_filename = argv[2]
# admixture_filename = argv[3]
# ASW_ids_filename = argv[4]
# out_filename = argv[5]
# kmer = int(argv[6])


def test_chromosomes_sliding_windows(kmer, seed, num_windows, num_test, train_filename, test_filename,
	all_admixture_filename, admixture_filename, out_filename):

	np.random.seed(seed)
	print('performing inference with window size: {}'.format(kmer))

	#### Read in strands ####
	print('loading data')
	train_chr_strands = pd.read_csv(train_filename, sep='\t')
	train_chr_strands.drop_duplicates(inplace=True)
	test_chr_strands_source = pd.read_csv(test_filename, sep='\t', header=0, comment='#')

	identifiers_to_use = list(test_chr_strands_source.columns[9:])
	print('all ids before cutoff ({} total): {}'.format(len(identifiers_to_use), identifiers_to_use))
	if num_test != -1:
		assert num_test <= len(identifiers_to_use)
		identifiers_to_use = identifiers_to_use[:num_test]
	print('ids after cutoff ({} total): {}'.format(len(identifiers_to_use), identifiers_to_use))

	test_chr_strands_list = [test_chr_strands_source[identifier].tolist() for identifier in identifiers_to_use]
	print('loaded {} test chromosomes'.format(num_test))

	overall_start_time = time.time()

	#### Read in Admixture results for all chromosomes and ith chromosome
	print('loading Admixture')
	all_chromosomes_ADMIX = pd.read_csv(all_admixture_filename, header=None, delim_whitespace=True)
	single_chromosome_ADMIX_unordered = pd.read_csv(admixture_filename, header=None, delim_whitespace=True)
	
	### Order columns of all chromosomes ADMIX in order of increasing population frequency
	all_chromosomes_ADMIX_ordered = all_chromosomes_ADMIX.reindex_axis(all_chromosomes_ADMIX.mean().sort_values(ascending=False).index, axis=1)
	all_chromosomes_ADMIX_ordered.columns = ['pop{}'.format(i) for i in range(1, len(all_chromosomes_ADMIX_ordered.columns) + 1)]

	### Subtract the sum of each column from the sum of pop1 proportions in all chromosomes ADMIX
	### TODO: examine how this will work when using multiple chromosomes
	admix_chr_diff_all = single_chromosome_ADMIX_unordered.apply(lambda col: all_chromosomes_ADMIX_ordered['pop1'] - col, axis=0).sum()


	### Sort ADMIX columns (each corresponding to a different inferred ancestral population) based on how much they differ from all_chromosomes_ADMIX pop1 values
	ADMIX = single_chromosome_ADMIX_unordered.reindex_axis(admix_chr_diff_all.sort_values(ascending=True).index, axis=1)

	### columns: [0, 1, 2] -> ['pop1', 'pop2', 'pop3'], where freq(pop1) < freq(pop2) < freq(pop3)
	ADMIX.columns = ['pop{}'.format(i) for i in range(1, len(ADMIX.columns) + 1)]

	### add labels to each row (sample chromosome)
	ADMIX.rename(index={i : colname for i, colname in enumerate(train_chr_strands.columns[9:])}, inplace=True)

	### calculate overall proportion of each population
	prob_pops = ADMIX.sum() / len(ADMIX)

	train_chr_strands_with_test = train_chr_strands.copy()
	for identifier, SNPs in zip(identifiers_to_use, test_chr_strands_list):
		train_chr_strands_with_test[identifier] = SNPs

	train_chr_strands_top = HMM.select_informative_SNPs(train_chr_strands, prob_pops, ADMIX, diff_quantile=0.1)
	train_chr_strands_with_test_top = train_chr_strands_with_test.ix[train_chr_strands_top.index, :]

	all_windows = dict()
	all_windows_prob = dict()

	for window in range(num_windows):
		print('\nperforming inference for window {} of {}'.format(window + 1, num_windows))

		new_strands_with_test_top = train_chr_strands_with_test_top.iloc[window:] # drop first $window rows (do nothing for window=0)

		### group SNPs (combine rows) into k-size windows
		chr_strand_with_test_substrings = HMM.create_kmer_haplotypes(new_strands_with_test_top, kmer)
		# chr_strand_substrings = HMM.create_kmer_haplotypes(new_strands_top, kmer) #new

		### generate all k-length bit strings. 
		all_haplotypes = HMM.all_kmer_haplotypes(kmer)

		ADMIX_with_test = ADMIX.copy()
		for identifier in identifiers_to_use:
			ADMIX_with_test.loc[identifier] = (ADMIX.sum() / len(ADMIX))

		log2_emission_matrices = HMM.create_emission_matrix(chr_strand_with_test_substrings, ADMIX_with_test, prob_pops, all_haplotypes)
		# log2_emission_matrices = HMM.create_emission_matrix(chr_strand_substrings, ADMIX, prob_pops, all_haplotypes)
		### getting empty entries in emission matrices when leaving out test chrs

		### local ancestry inference ###
		inferences = HMM.both_directions_local_ancestry_prob(chr_strand_with_test_substrings, ADMIX, identifiers_to_use, log2_emission_matrices, kmer, recomb_rate=0.001)

		### combine forward ancestry with backward ancestry, with a value of 'unknown' for positions that differ between the forward and backward inferences
		all_windows[window] = inferences['both_ancestry']
		all_windows_prob[window] = inferences['both_probs']

	all_prob = pd.concat(all_windows_prob.values())
	### The next two lines are used for filling in uninformative SNPs before the first informative SNP and after the last one.
	### Using values of 1 means the resulting confidence for these areas will be cut in half, but the ancestry calls will not be changed.
	all_prob.loc[-1] = 1
	all_prob.loc[train_chr_strands_with_test.index[-1] + 1] = 1
	all_prob.sort_index(inplace=True)


	complete_ancestry_probs = pd.DataFrame(train_chr_strands_with_test['POS'])
	### add columns for each person.
	# TODO: dtypes for these columns will be 'object' after adding 'pop1'/'pop2' (dtype is float64 upon creation) - not too efficient.
	#       Can't use fixed-length strings (https://stackoverflow.com/a/34884078/5377941) so only other option is to represent population with a number
	complete_ancestry_probs = complete_ancestry_probs.reindex(columns=['POS'] + list(all_prob.columns[2:]))
	### this will hold population inferences, based on contents of complete_ancestry_probs
	complete_ancestry = complete_ancestry_probs.copy()

	### fill in informative SNPs in complete_ancestry_probs
	### this is working (and takes about two minutes with 10 samples on chromosome 22). TODO: profile
	### skip [-1] window and start with [-1, all_prob.index[0]] window
	print('\n\naveraging windows for informative SNPs')
	inf_start = time.time()
	for index, windows in enumerate(itertools.islice(HMM.get_sliding_windows(all_prob.index, num_windows), 1, None)):
		if index % 1000 == 0:
			print('on window {} of ~{}'.format(index, all_prob.shape[0]))
		# avg = all_prob.loc[windows].iloc[:, 2:].mean() ### get the relevant contiguous windows in all_prob, then discard 'POS_start/end' columns
		avg = all_prob.loc[windows].iloc[:, 2:].prod() ### averaging fix
		complete_ancestry_probs.ix[windows[-1], 1:] = avg
	inf_end = time.time()
	print('\ninformative SNPs took {:.1f} seconds'.format(inf_end - inf_start))

	### fill in uninformative SNPs by averaging the relevant windows
	### skip [-1] window and start with [-1, all_prob.index[0]] window
	print('\n\naveraging windows for uninformative SNPs')
	uninf_start = time.time()
	for index, windows in enumerate(itertools.islice(HMM.get_sliding_windows(all_prob.index, num_windows - 1), 1, None)):
		if index % 1000 == 0:
			print('on window {} of ~{}'.format(index, all_prob.shape[0]))
		# avg = all_prob.loc[windows].iloc[:, 2:].mean()
		avg = all_prob.loc[windows].iloc[:, 2:].prod() ### averaging fix
		### Seems like I can't broadcast one row to multiple rows.
		### Some more efficient methods are listed here: https://stackoverflow.com/questions/18771963/pandas-efficient-dataframe-set-row
		### TODO: check performance for this section
		### TODO: try:
		###			complete_ancestry_probs.ix[windows[-2] + 1 : windows[-1], 1:] = itertools.repeat(avg, windows[-1] - (windows[-2] + 1))
		for row in range(windows[-2] + 1, windows[-1]):
			complete_ancestry_probs.ix[row, 1:] = avg
	uninf_end = time.time()
	print('\nuninformative SNPs took {:.1f} seconds'.format(uninf_end - uninf_start))

	### drop final row, added in line above: `all_prob.loc[train_chr_strands_with_test.index[-1] + 1] = 1`
	complete_ancestry_probs.drop(complete_ancestry_probs.index[-1], inplace=True)

	for column in complete_ancestry_probs.columns[1:]: ### for each sample:
		complete_ancestry[column] = np.where(complete_ancestry_probs[column] > 1, 'pop2', 'pop1')

	print('saving inferences to {}.txt'.format(out_filename))
	complete_ancestry.to_csv(out_filename + '.txt', sep='\t', index=False)

	overall_end_time = time.time()
	time_diff = round(overall_end_time - overall_start_time, 1)
	print('computed ancestry for {} chromosomes in {} seconds'.format(len(identifiers_to_use), time_diff))

	# complete_ancestry_top = complete_ancestry.ix[chr_strands_with_test_top.index]


if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--kmer', help='window size - this is also used as the number of sliding windows', type=int, default=5)
	parser.add_argument('--num_windows', help='number of sliding windows to use. If not specified, the '
		'value for kmer is used', default=None)
	parser.add_argument('--seed', help='random seed to use.', type=int, default=0)
	parser.add_argument('--num_test', help='number of test chromosomes to evaluate (default all)', type=int, default=-1)

	parser.add_argument('--reference_filename', help='path to file that contains training set', type=str, required=True)
	parser.add_argument('--all_admix_filename', help='path to overall admix filename', type=str, required=True)
	parser.add_argument('--chrom_admix_filename', help='path to admix filename for this chromosome', type=str, required=True)
	parser.add_argument('--test_filename', help='path to file that contains chromosomes for which ancestry should be inferred', type=str, required=True)

	args = parser.parse_args()

	kmer = args.kmer
	seed = args.seed
	num_test = args.num_test
	num_windows = args.num_windows or kmer

	train_filename = args.reference_filename
	test_filename = args.test_filename
	all_admixture_filename = args.all_admix_filename
	admixture_filename = args.chrom_admix_filename

	train_basename = os.path.basename(train_filename)
	test_basename = os.path.basename(test_filename)

	out_filename = '/home/greg/School/popgen/results/source={}_test={}'.format(train_basename, test_basename)

	test_chromosomes_sliding_windows(
		kmer,
		seed,
		num_windows,
		num_test,
		train_filename,
		test_filename,
		all_admixture_filename,
		admixture_filename,
		out_filename)
