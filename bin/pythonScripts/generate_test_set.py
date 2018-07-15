import numpy as np
import pandas as pd

import itertools
import pdb

#### Function to create test chromosomes ####
def create_test_chromosomes(chr_strands, pops, ref_IDs, num_switch):
	test_pos_chromosome = chr_strands['POS'].sample(num_switch, replace=False) ### sample from SNPs
	assert(len(test_pos_chromosome) == num_switch)

	test_pos_chromosome.sort_values(inplace=True)
	test_pos_chromosome.reset_index(drop=True, inplace=True)

	test_chr_idx_stop = pd.Series(range(num_switch)).map(lambda x: np.where(chr_strands['POS'] == test_pos_chromosome[x])[0][0]+1)
	test_chr_idx_start = pd.concat([pd.Series(0), test_chr_idx_stop], ignore_index=True)
	test_chr_idx_stop = pd.concat([test_chr_idx_stop, pd.Series(len(chr_strands))], ignore_index=True)
	start = test_chr_idx_start
	stop = test_chr_idx_stop

	test_chr = np.zeros(len(chr_strands), dtype=int)
	true_chr = np.array([''] * len(chr_strands), dtype='<U4') ### need to specify max string length
	for j in range(len(pops)):
		for i in range(j, len(test_chr_idx_start), len(pops)):
			test_chr[start[i] : stop[i]] = chr_strands[ref_IDs[j]][start[i] : stop[i]]
			true_chr[start[i] : stop[i]] = pops[j]

	return [test_chr, true_chr]


def create_test_chromosome_set(chr_strands, pops_dict, num_chromosomes, num_switch=8):
	### pops_dictionary should map population name to list of haploid IDs corresponding to that population
	### if num_chromosomes == -1, as many chromosomes as possible will be created

	ID_pairs = []
	### copy dict values (lists) since they will be modified; two copies of each haplotype, so ignore one
	### remove IDs not in chr_strands
	# pops_dict = {key: list(filter(lambda v: v[-2:] != '_2' and v in chr_strands.columns, val)) for key, val in pops_dictionary.items()}

	pops_ordered = ['pop1', 'pop2']
	for i in itertools.count():
		### keep getting a random ID from each population, until either one of the ID lists is empty, or num_chromosomes is exceeded
		if [] in pops_dict.values() or (num_chromosomes != -1 and i >= num_chromosomes):
			break
		ID_pairs.append([pops_dict[pop].pop(np.random.randint(0, len(pops_dict[pop]))) for pop in pops_ordered]) ### get a random ID from each pop

	return [create_test_chromosomes(chr_strands, pops_ordered, pair, num_switch) for pair in ID_pairs], ID_pairs


def get_pop_IDs(filename):
	ret = []
	with open(filename, 'r') as f:
		for line in f:
			line = line.strip('\n')
			if line != '':
				ret.append(line)
	return ret


if __name__ == '__main__':
	np.random.seed(0)
	
	sourcepop1 = "CEU"
	sourcepop2 = "YRI"
	sourcepop1_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop1)
	sourcepop2_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop2) 
	
	populations = "ASW_CEU_YRI"
	chr_number = 22
	allele_filename = 'pops_data/{}_Data/Chr{}/tmp/chr22.phase3.{}.SNPs.allele.vcf'.format(populations, chr_number, populations)
	# test1_ids_filename = '/home/greg/School/popgen/SampleIDs/CEU_Sample_IDs_haploid.txt'
	# test2_ids_filename = '/home/greg/School/popgen/SampleIDs/YRI_Sample_IDs_haploid.txt'
	# test_source_filename = '/home/greg/School/popgen/data/chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt'

	sourcepop1_ids = get_pop_IDs(sourcepop1_ids_filename)
	sourcepop2_ids = get_pop_IDs(sourcepop2_ids_filename)

	num_admixed_chromosomes = 10

	all_chr_strands = pd.read_csv(allele_filename, sep='\t', header=0, comment='#')
#	source_df = pd.read_csv(allele_filename, sep='\t', header=0, comment='#')
	all_chr_strands.drop_duplicates(inplace=True)

	# The function below ensures thats the sourcepop1_ids and pop2_ids consist ONLY of individuals (N12345_0, N12345_1 etc) from the allele.vcf file
	# that are present in sourcepop1_ids_filename and sourcepop1_ids_filename (the haploid.txt files for the two populations)
	sourcepop1_ids, sourcepop2_ids = map(lambda ids: list(filter(lambda ID: ID in all_chr_strands.columns, ids)),[sourcepop1_ids, sourcepop2_ids])
	
	test_chr = create_test_chromosome_set(all_chr_strands, {'pop1' : sourcepop1_ids, 'pop2' : sourcepop2_ids}, num_admixed_chromosomes)

	ancestry_df = pd.DataFrame(all_chr_strands[all_chr_strands.columns[:9]])
	chrom_df = pd.DataFrame(all_chr_strands[all_chr_strands.columns[:9]])
	
	for (test, true), id_pair in zip(*test_chr):
		name = '-'.join(id_pair).replace('_', '-')
		ancestry_df[name] = true
		chrom_df[name] = test

	ancestry_df.to_csv('test_input/CEU_YRI_test_cases_true_ancestry_try2_mm.csv', sep='\t', index=False)
	chrom_df.to_csv('test_input/CEU_YRI_test_cases_chrom_try2_mm.csv', sep='\t', index=False)