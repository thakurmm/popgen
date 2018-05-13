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


def create_test_chromosome_set(chr_strands, pops_dictionary, num_chromosomes, num_switch=8):
	### pops_dictionary should map population name to list of haploid IDs corresponding to that population
	### if num_chromosomes == -1, as many chromosomes as possible will be created

	ID_pairs = []
	### copy dict values (lists) since they will be modified; two copies of each haplotype, so ignore one
	### remove IDs not in chr_strands
	pops_dict = {key: list(filter(lambda v: v[-2:] != '_2' and v in chr_strands.columns, val)) for key, val in pops_dictionary.items()}
	pops_ordered = sorted(list(pops_dictionary.keys()))
	assert len(pops_ordered) == 2
	for i in itertools.count():
		### keep getting a random ID from each population, until either one of the ID lists is empty, or num_chromosomes is exceeded
		if [] in pops_dict.values() or (num_chromosomes != -1 and i >= num_chromosomes):
			break
		ID_pairs.append([pops_dict[pop].pop(np.random.randint(0, len(pops_dict[pop]))) for pop in pops_ordered]) ### get a random ID from each pop

	return [create_test_chromosomes(chr_strands, pops_ordered, pair, num_switch) for pair in ID_pairs], ID_pairs


def read_ids(filename):
	with open(filename, 'r') as f:
		return map(str.strip, f.readlines())

if __name__ == '__main__':
	np.random.seed(0)
	test1_ids_filename = '/home/greg/School/popgen/SampleIDs/CEU_Sample_IDs_haploid.txt'
	test2_ids_filename = '/home/greg/School/popgen/SampleIDs/YRI_Sample_IDs_haploid.txt'
	test_source_filename = '/home/greg/School/popgen/data/chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt'

	test1_ids = read_ids(test1_ids_filename)
	test2_ids = read_ids(test2_ids_filename)

	num_test = 10

	source_df = pd.read_csv(test_source_filename, sep='\t', header=0, comment='#')
	test_chr = create_test_chromosome_set(source_df, {'pop1' : test1_ids, 'pop2' : test2_ids}, num_test)

	ancestry_df = pd.DataFrame(source_df[source_df.columns[:9]])
	chrom_df = pd.DataFrame(source_df[source_df.columns[:9]])
	
	for (test, true), id_pair in zip(*test_chr):
		name = '-'.join(id_pair).replace('_', '-')
		ancestry_df[name] = true
		chrom_df[name] = test

	ancestry_df.to_csv('test_input/CEU_YRI_test_cases_true_ancestry_try2.csv', sep='\t', index=False)
	chrom_df.to_csv('test_input/CEU_YRI_test_cases_chrom_try2.csv', sep='\t', index=False)