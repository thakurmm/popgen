import ggplot ### conda install -c conda-forge ggplot
import numpy as np
import pandas as pd
import pymongo

from collections import deque
import itertools
import pdb
import sys

def select_informative_SNPs(chr_strands: pd.DataFrame, prob_pops: pd.Series, ADMIX2: pd.DataFrame, diff_quantile = 0.1) -> pd.DataFrame:
	### look at using FST instead of absolute difference in conditional prob
	### unevenly distributed across genome
	### could pick the SNP with highest FST in each window (maybe 8,000 - 10,000 bp)
	### or divide into windows, find SNPs with FST above certain percentile (90, 95, 99?)

	### look at using P file (for each SNP, gives ancestry proportions) from ADMIX to do the same calculations before, slightly less accurate but much faster once you have it
	print('Finding informative SNPs')

	alleles = [0, 1]

	conditional_prob_matrices = np.ones(shape=(len(prob_pops), len(alleles), len(chr_strands)), )

	for i in range(len(alleles)):
		A = alleles[i]

		allele_df = pd.DataFrame({
			### TODO: when A is 1, below is just rowsums(chr_strands); when A is 0, below is len(chr_strands) - rowsums(chr_strands)
			'row_sums': (chr_strands.ix[:,9:] == A).astype(int).sum(axis=1).values,
			}, index=chr_strands['ID'])

		allele_df['prob_allele'] = allele_df['row_sums'] / len(ADMIX2)

		for j in range(len(prob_pops)):
			## sum_over_individuals(#Allele[i] per Individual * %Pop[j] per Individual)

			### look at squaring ADMIX matrix (50/50 individual -> 25/25, only counts for half a person, where as 100/0 -> 100/0, higher weighting for individuals with unequal proportions)
			pop_sum = (chr_strands.ix[:,9:] == A).mul(ADMIX2.ix[:, j], axis=1).sum(axis=1) ### http://stackoverflow.com/questions/23117756/python-numpy-or-pandas-equivalent-of-the-r-function-sweep

			### http://stackoverflow.com/questions/35666272/equivalent-of-r-ifelse-in-python-pandas-compare-string-columns
			# pop_weighted_mean = np.where((allele_df['row_sums'] == 0), 0, pop_sum.reset_index(drop=True)/allele_df['row_sums'].reset_index(drop=True))
			pop_weighted_mean = np.where((allele_df['row_sums'] == 0), 0, pop_sum.values/allele_df['row_sums'].values)

			print('\npop_weighted_mean:')
			print(pop_weighted_mean)

			conditional_prob_matrices[j, i, :] = (pop_weighted_mean * allele_df['prob_allele']) / prob_pops[j]

	conditional_prob_differences = abs(conditional_prob_matrices[0, 0, :] - conditional_prob_matrices[1, 0, :])
	diff_quantile_value = np.percentile(conditional_prob_differences, max(diff_quantile, 1 - diff_quantile) * 100) ### http://stackoverflow.com/questions/18580461/eliminating-all-data-over-a-given-percentile
	top_conditional_prob_diff = np.where(abs(conditional_prob_differences) >= diff_quantile_value)[0]

	# print('\ninformative SNPs:')
	# print(chr_strands.ix[top_conditional_prob_diff, :])

	print('\ndiff quantile value is:')
	print(diff_quantile_value)

	print('number of SNPs selected:')
	print(len(chr_strands.ix[top_conditional_prob_diff, :]))

	return chr_strands.ix[top_conditional_prob_diff, :]


def all_kmer_haplotypes(k):
	# Uncomment the below line to include "N" as a SNP possibility in the haplotypes
	# return list(map(lambda x: ''.join(x), itertools.product('01N', repeat=k)))
	return list(map(lambda x: ''.join(x), itertools.product('01', repeat=k)))


def split_string_logic(text: str, k: int) -> list:
	### split text string into substrings of length k
	extended_text = text + text[ : (k - (len(text) % k)) % k] ### repeat characters from the beginning of text until its length is divisible by k
	return (extended_text[i : i+k] for i in range(0, len(extended_text), k))

def get_sliding_windows(iterable, n):
	### given ['a', 'b', 'c', 'd', 'e'] and n=3, yields
	### ['a']
	### ['a', 'b']
	### ['a', 'b', 'c']
	### ['b', 'c', 'd']
	### ['c', 'd', 'e']
	it = iter(iterable)
	window = deque()
	for _ in range(n-1):
		window.append(next(it))
		yield list(window)

	while len(window) > 0:
		next_item = next(it)
		window.append(next_item)		
		yield list(window)

		window.popleft()

# @profile
def create_kmer_haplotypes(chr_strands: pd.DataFrame, kmer: int) -> pd.DataFrame:
	print('creating {}-mer haplotypes'.format(kmer))

	if kmer != 1:
		remainder = len(chr_strands) % kmer
		if remainder != 0:
			chr_strands.drop(chr_strands.tail(remainder).index, inplace=True) ### drop last $remainder rows

		chr_strand_substrings = chr_strands.iloc[np.arange(0, len(chr_strands), kmer), :].copy()
		chr_strand_substrings.rename(columns={'POS' : 'POS_start'}, inplace=True)
		chr_strand_substrings['POS_end'] = chr_strands['POS'].iloc[np.arange(kmer-1, len(chr_strands), kmer)].values
		cols = chr_strand_substrings.columns.tolist()
		cols = cols[:2] + [cols[-1]] + cols[2:-1]
		chr_strand_substrings = chr_strand_substrings[cols]

		# all_cols = pd.concat(col[1] for col in chr_strands.iloc[:, 9:].iteritems()).map(str).str.cat() ### <--- no errors but it changes the output...
		### Rather than casting each column to string and then concatenating all the values (in the old version, this
		### happens once for each column in the for loop below), this first concatenates every column together, casts
		### to string, and concatenates.

		for i in range(10, len(chr_strand_substrings.columns)):
			chr_strands_as_strings = chr_strands.iloc[:, i-1].map(str).str.cat()
			# start = (i - 10) * len(chr_strands)
			# stop = start + len(chr_strands)
			# chr_strands_as_strings = all_cols[start:stop]
			chr_strand_substrings.iloc[:, i] = pd.Series(split_string_logic(chr_strands_as_strings, kmer)).values ### <-- 74% of time spent here

	else:
		chr_strand_substrings = chr_strands
		chr_strand_substrings.rename(columns={'POS' : 'POS_start'}, inplace=True) ### may give SettingWithCopyWarning
		chr_strand_substrings['POS_end'] = chr_strand_substrings['POS_start'].copy()
		cols = chr_strand_substrings.columns.tolist()
		cols = cols[:2] + [cols[-1]] + cols[2:-1] ### move last column (POS_end column) to third column's position and shift everything else over
		chr_strand_substrings = chr_strand_substrings[cols]		

	return chr_strand_substrings


def create_emission_matrix(chr_strand_substrings, ADMIX2, prob_pops, kmer_haplotypes):
	'''
	- each 'page' (determined by first index) in emission_matrix refers to a specific population
	- each row (determined by second index) refers to a haplotype
	- each column (third index) refers to a SNP
	'''
	print('creating emission matrix')

	emission_matrices = np.ones(shape=(len(ADMIX2.columns), len(kmer_haplotypes), len(chr_strand_substrings)))

	for i in range(len(kmer_haplotypes)):
		hap = kmer_haplotypes[i]
		matches = (chr_strand_substrings.ix[:, 10:] == hap).astype(int) ### <-- 70% of time in function is spent here
		haplotype_df = pd.DataFrame({
			'row_sums': matches.sum(axis=1).values,
			}, index=chr_strand_substrings['ID'])
		haplotype_df['prob_haplotype'] = haplotype_df['row_sums'] / len(ADMIX2)

		for j in range(len(ADMIX2.columns)):
			pop_sum = matches.mul(ADMIX2.ix[:, j], axis=1).sum(axis=1)
			pop_weighted_mean = np.where((haplotype_df['row_sums'] == 0), 0, pop_sum.values / haplotype_df['row_sums'].values)
			emission_matrices[j, i, :] = (pop_weighted_mean * haplotype_df['prob_haplotype']) / prob_pops[j]

	return np.log2(emission_matrices)


def viterbi_init(next_haplotype, init_prob, emission_matrix, kmer_haplotypes):
	haplotype_idx = np.where(np.array(kmer_haplotypes) == next_haplotype)[0][0]

	next_prob = init_prob
	for pop_idx in range(len(next_prob)):
		next_prob.iloc[pop_idx] = emission_matrix[pop_idx, haplotype_idx] + init_prob.iloc[pop_idx]

	return next_prob

def viterbi(next_haplotype, prev_prob, emission_matrix, transition_matrix, kmer_haplotypes):
	haplotype_idx = np.where(kmer_haplotypes == next_haplotype)[0][0]

	next_prob = prev_prob.copy()
	for pop_idx in range(len(next_prob)):
		next_prob.iloc[pop_idx] = emission_matrix[pop_idx, haplotype_idx] + max(prev_prob + transition_matrix[pop_idx, :])

	return next_prob		


def HMM_prob_one_direction(chr_strand_substrings_array, init_prob, log2_emission_matrices, log2_transition_matrix, kmer_haplotypes, direction):
	chr_states = chr_strand_substrings_array.copy() ### copied to avoid SettingWithCopyWarning for 'chr_states[0] = np.argmax(next_prob)' line
	chr_states_prob = chr_strand_substrings_array.copy()

	### determine which direction to go when iterating over chr_strand_substrings_array
	if direction == 'forward':
		val_range = range(0, len(chr_states), 1)
	elif direction == 'backward':
		val_range = range(len(chr_states)-1, -1, -1)
	else:
		raise Exception('direction should be either "forward" or "backward"')

	print('direction is {}'.format(direction))

	start = val_range[0]

	next_prob = viterbi_init(chr_strand_substrings_array.iloc[start], np.log2(init_prob), log2_emission_matrices[:,:,start], kmer_haplotypes)

	chr_states_prob.iloc[start] = next_prob[0] / next_prob[1]
	chr_states.iloc[start] = np.argmax(next_prob)
	prev_prob = next_prob

	for j in val_range[1:]:
		next_prob = viterbi(chr_strand_substrings_array.iloc[j], prev_prob, log2_emission_matrices[:, :, j], log2_transition_matrix, kmer_haplotypes)
		chr_states_prob.iloc[j] = next_prob[0] / next_prob[1]
		if np.isnan(chr_states_prob.iloc[j]):
			print('nan at j={}'.format(j))
		chr_states.iloc[j] = np.argmax(next_prob)
		prev_prob = next_prob

	print('finished {} direction'.format(direction))
	return {'ancestry' : chr_states, 'probs' : chr_states_prob}


def both_directions_local_ancestry_prob(chr_strand_substrings, ADMIX, individual_IDs, log2_emission_matrices, kmer, recombination_rate):
	'''
	individual_IDs is a Python list of IDs of individuals to perform LAI for
	'''

	all_haplotypes = np.array(all_kmer_haplotypes(kmer))

	### initially contains column of NaN for each id in individual_IDs
	chr_states_forward = pd.concat([chr_strand_substrings.loc[:, ['POS_start', 'POS_end']]] +
		[pd.Series().rename(name) for name in individual_IDs], axis=1)

	chr_states_forward_prob = chr_states_forward.copy()
	chr_states_backward = chr_states_forward.copy()
	chr_states_backward_prob = chr_states_forward.copy()
	chr_states_both_directions = chr_states_forward.copy()
	chr_states_both_directions_prob = chr_states_forward.copy()

	for name in individual_IDs:
		print('\nrunning Viterbi for {}'.format(name))
		# init_prob = ADMIX.ix[name, :]
		# init_prob = pd.Series([0.5, 0.5])
		init_prob = ADMIX.sum() / len(ADMIX)

		### transition matrix for each individual
		transition_matrix = np.zeros(shape=(len(ADMIX.columns), len(ADMIX.columns)))

		npops = len(ADMIX.columns)

		for pop_index in range(npops):

			### change nth column, except for element on diagonal

			row_ix = [True] * npops
			row_ix[pop_index] = False

			col_ix = [False] * npops
			col_ix[pop_index] = True

			transition_matrix[np.ix_(row_ix, col_ix)] = recombination_rate * init_prob[pop_index]

		np.fill_diagonal(transition_matrix, 1 - np.sum(transition_matrix, axis=1))
		log2_transition_matrix = np.log2(transition_matrix)

		forward_HMM_results = HMM_prob_one_direction(chr_strand_substrings[name], init_prob, log2_emission_matrices, log2_transition_matrix, all_haplotypes, 'forward')

		chr_states_forward[name] = forward_HMM_results['ancestry']
		chr_states_forward_prob[name] = forward_HMM_results['probs']

		backward_HMM_results = HMM_prob_one_direction(chr_strand_substrings[name], init_prob, log2_emission_matrices, log2_transition_matrix, all_haplotypes, 'backward')

		chr_states_backward[name] = backward_HMM_results['ancestry']
		chr_states_backward_prob[name] = backward_HMM_results['probs']

		chr_states_both_directions[name] = np.where(chr_states_forward[name] == chr_states_backward[name], chr_states_backward[name], 'unknown')

		combined_prob_results = pd.concat([forward_HMM_results['probs'], backward_HMM_results['probs']], axis=1)
		# combined_prob_results_averages = combined_prob_results.mean(axis=1)
		combined_prob_results_averages = combined_prob_results.prod(axis=1) ### averaging fix
		chr_states_both_directions_prob[name] = combined_prob_results ### NOTE: value > 1 -> prob(pop2) > prob(pop1)

	return({
		'both_ancestry' : chr_states_both_directions,
		'both_probs' : chr_states_both_directions_prob,
		'forward_ancestry' : chr_states_forward,
		'forward_probs' : chr_states_forward_prob,
		'backward_ancestry' : chr_states_backward,
		'backward_probs' : chr_states_backward_prob,
		})


def evaluate_inferences_full(LA_df, LA_df_top, test_to_true, rounding=5):
	ret = {
		'overall' : {
			'acc' : {},
		},

		'informative' : {
			'acc' : {},
		},
	}

	overall_accuracies = []
	for test, true in test_to_true.items():
		acc = (LA_df[test] == LA_df[true]).mean()
		overall_accuracies.append(acc)

	ret['overall']['acc']['avg'] = round(np.mean(overall_accuracies), rounding)
	ret['overall']['acc']['min'] = round(min(overall_accuracies), rounding)
	ret['overall']['acc']['max'] = round(max(overall_accuracies), rounding)
	ret['overall']['acc']['var'] = round(np.var(overall_accuracies), rounding)

	informative_accuracies = []
	for test, true in test_to_true.items():
		acc = (LA_df_top[test] == LA_df_top[true]).mean()
		informative_accuracies.append(acc)

	ret['informative']['acc']['avg'] = round(np.mean(informative_accuracies), rounding)
	ret['informative']['acc']['min'] = round(min(informative_accuracies), rounding)
	ret['informative']['acc']['max'] = round(max(informative_accuracies), rounding)
	ret['informative']['acc']['var'] = round(np.var(informative_accuracies), rounding)

	return ret




def write_to_db(
		evaluation: dict,
		time: int,
		window_size: int,
		sliding_window: bool,
		random_seed: int,
		batch_size: int):
	client = pymongo.MongoClient()
	db = client.lai
	results = db.results

	results.insert_one(
		{
			'overall_acc' : evaluation['overall']['acc'],
			'informative_SNPs_acc' : evaluation['informative']['acc'],
			'time' : time,
			'window_size' : window_size,
			'sliding_window' : sliding_window,
			'random_seed' : random_seed,
			'batch_size' : batch_size
		})


def evaluate_inferences(LA_df, test_to_true, miss_percentile=95, rounding=5):
	### - LA_df: contains inferred and actual ancestries for each chromosome
	### - test_to_true: map test IDs to true IDs, i.e.
	### 	{'test_1' : 'true_1', 'test_2' : 'true_2', ...}
	### - miss_percentile: after calculating ratio of misses for each position on the chromosome, only return positions
	###		with a miss ratio above this percentile
	### - rounding: number of digits after the decimal to round to
	ret = {
		'acc' : {}, ### statistics based purely on checking whether inferred values are equal to correct values
		'unknown' : {}, ### statistics related to the frequency of ambiguous calls (not based on correct values at all)
		'incorrect' : {}, ### statistics based on incorrect calls that were unambiguous
		'pos' : {}, ### statistics based on position, i.e. for which positions was it difficult to predict ancestry
	}

	accuracies = []
	for test, true in test_to_true.items():
		acc = len(np.where(LA_df[test] == LA_df[true])[0]) / len(LA_df)
		accuracies.append(acc)
	ret['acc']['avg'] = sum(accuracies) / len(accuracies)
	ret['acc']['min'] = min(accuracies)
	ret['acc']['max'] = max(accuracies)
	ret['acc']['var'] = np.var(accuracies)

	unknowns = []
	for test in test_to_true:
		unk = LA_df[test].value_counts()['unknown'] / len(LA_df)
		unknowns.append(unk)
	ret['unknown']['avg'] = sum(unknowns) / len(unknowns)
	ret['unknown']['min'] = min(unknowns)
	ret['unknown']['max'] = max(unknowns)
	ret['unknown']['var'] = np.var(unknowns)

	incorrect = []
	for acc, unk in zip(accuracies, unknowns):
		inc = 1 - acc - unk
		incorrect.append(inc)
	ret['incorrect']['avg'] = sum(incorrect) / len(incorrect)
	ret['incorrect']['min'] = min(incorrect)
	ret['incorrect']['max'] = max(incorrect)
	ret['incorrect']['var'] = np.var(incorrect)

	matches = LA_df[['POS_start', 'POS_end']].copy()
	for index, (test, true) in enumerate(test_to_true.items()):
		matches['match_{}'.format(index)] = np.where(LA_df[test] == LA_df[true], 0, 1) ### 1 indicates a mismatch

	percent_misses = matches.filter(regex=r'match_\d').mean(axis=1) ### mean of each row
	### could do a lot of different things here
	### normalize by (POS_end - POS_start)? (more likely to make incorrect call on larger window)
	percentile = np.percentile(percent_misses, miss_percentile)
	misses_idx = np.where(percent_misses >= percentile)[0]
	top_misses = LA_df[['POS_start', 'POS_end']].iloc[misses_idx] ### get positions with the most misses
	top_misses['percent_misses'] = percent_misses.ix[top_misses.index] ### include miss percentage for each position

	ret['pos']['top_misses'] = top_misses

	return ret


def plot_and_save_local_ancestry_top(df, kmer, image_filename, num_chromosomes, y_scale=0.5):
	return _plot_and_save_local_ancestry(df, kmer, image_filename, num_chromosomes, ['POS_start', 'POS_end'], 'POS_start', y_scale)


def plot_and_save_local_ancestry_full(df, kmer, image_filename, num_chromosomes, y_scale=0.5):
	return _plot_and_save_local_ancestry(df, kmer, image_filename, num_chromosomes, ['POS'], 'POS', y_scale)


def _plot_and_save_local_ancestry(df, kmer, image_filename, num_chromosomes, id_vars, x_axis, y_scale):
	print('saving plot as: {}'.format(image_filename))
	var_name='chromosome'

	local_ancestry_df_long = pd.melt(df, id_vars=id_vars, var_name=var_name, value_name='estimated_ancestry')

	new_names = {}
	for i in range(1, num_chromosomes + 1):
		new_names['test_{}'.format(i)] = 2*i - 2 * y_scale
		new_names['true_{}'.format(i)] = 2*i - 1 * y_scale

	for key, value in new_names.items():
		local_ancestry_df_long.replace(key, value, inplace=True)

	plot = ggplot.ggplot(ggplot.aes(x=x_axis, y=var_name, color='estimated_ancestry'), data=local_ancestry_df_long) \
		+ ggplot.geom_point() \
		+ ggplot.scale_y_continuous(labels=list(new_names.keys()), breaks=list(new_names.values())) \
		+ ggplot.scale_color_manual(values=['#FF0000', '#0000FF', '#73008C']) \
		+ ggplot.theme(plot_margin={'top':0.7, 'bottom':0.3}) ### TODO: this should depend on scale

	plot.save(image_filename)


# def report_evaluation_full(evaluation: dict, time: int, window_size: int, sliding: bool) -> str:
# 	overall_acc = evaluation['overall']['acc']
# 	informative_acc = evaluation['informative']['acc']
# 	ret = ''
# 	ret += 'window size: {}\n'.format(window_size)
# 	ret += 'used sliding windows: {}\n'.format(sliding)
# 	ret += 'time taken: {} seconds\n'.format(time)
# 	ret += '\n'
# 	ret += '=============\n'
# 	ret += 'overall accuracy\n'
# 	ret += '=============\n'
# 	ret += 'avg: {:.2f}%\n'.format(overall_acc['avg'] * 100)
# 	ret += 'var: {:.2f}%\n'.format(overall_acc['var'] * 100)
# 	ret += 'min: {:.2f}%\n'.format(overall_acc['min'] * 100)
# 	ret += 'max: {:.2f}%\n'.format(overall_acc['max'] * 100)
# 	ret += '\n'
# 	ret += 'informative SNP accuracy\n'
# 	ret += '===============\n'
# 	ret += 'avg: {:.2f}%\n'.format(informative_acc['avg'] * 100)
# 	ret += 'var: {:.2f}%\n'.format(informative_acc['var'] * 100)
# 	ret += 'min: {:.2f}%\n'.format(informative_acc['min'] * 100)
# 	ret += 'max: {:.2f}%\n'.format(informative_acc['max'] * 100)
# 	return ret


# def report_evaluation(evaluation):
# 	### take in the output from evaluate_inferences
# 	acc = evaluation['acc']
# 	unknown = evaluation['unknown']
# 	incorrect = evaluation['incorrect']
# 	ret = ''
# 	ret += 'correct calls\n'
# 	ret += '=============\n'
# 	ret += 'avg: {:.2f}%\n'.format(acc['avg'] * 100)
# 	ret += 'var: {:.2f}%\n'.format(acc['var'] * 100)
# 	ret += 'min: {:.2f}%\n'.format(acc['min'] * 100)
# 	ret += 'max: {:.2f}%\n'.format(acc['max'] * 100)
# 	ret += '\n'
# 	ret += 'ambiguous calls\n'
# 	ret += '===============\n'
# 	ret += 'avg: {:.2f}%\n'.format(unknown['avg'] * 100)
# 	ret += 'var: {:.2f}%\n'.format(unknown['var'] * 100)
# 	ret += 'min: {:.2f}%\n'.format(unknown['min'] * 100)
# 	ret += 'max: {:.2f}%\n'.format(unknown['max'] * 100)
# 	ret += '\n'
# 	ret += 'incorrect, unambiguous calls\n'
# 	ret += '============================\n'
# 	ret += 'avg: {:.2f}%\n'.format(incorrect['avg'] * 100)
# 	ret += 'var: {:.2f}%\n'.format(incorrect['var'] * 100)
# 	ret += 'min: {:.2f}%\n'.format(incorrect['min'] * 100)
# 	ret += 'max: {:.2f}%\n'.format(incorrect['max'] * 100)
# 	return ret

def evaluate(test, true, rounding=5):
   acc = {}

   columns_to_ignore = ['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
   ancestry_cols = list(filter(lambda x: x not in columns_to_ignore, test.columns))

   overall_accuracies = []

   for col_name in ancestry_cols:
       if col_name not in true:
           raise KeyError('true ancestry dataframe is missing ancestry for id: {}'.format(col_name))
       mean_acc = (test[col_name] == true[col_name]).mean()
       overall_accuracies.append(mean_acc)

   acc['avg'] = round(np.mean(overall_accuracies), rounding)
   acc['min'] = round(min(overall_accuracies), rounding)
   acc['max'] = round(max(overall_accuracies), rounding)
   acc['var'] = round(np.var(overall_accuracies), rounding)

   return acc


def new_plot_ancestry_with_correct_results(test, true, y_scale=0.5, image_filename=None):
   columns_to_ignore = ['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] ### we want only 'POS' and ancestry columns
   ancestry_cols = list(filter(lambda x: x not in columns_to_ignore, test.columns))

   merged = pd.DataFrame(test['POS'])
   for col_name in ancestry_cols:
       if col_name not in true:
           raise KeyError('true ancestry dataframe is missing ancestry for id: {}'.format(col_name))
       merged[col_name+'_test'] = test[col_name]
       merged[col_name+'_true'] = true[col_name]

   melted = pd.melt(merged, id_vars=['POS'], var_name='chromosome', value_name='ancestry')
   # the above takes merged from something like this:
   ###
   ### columns: POS sample1_test sample1_true sample2_test sample2_true
   ###                      111      pop1         pop1         pop2         pop1
   ###          124      pop1         pop1         pop2         pop1
   ###
   # to this: (spaces between rows added for clarity)
   ###
   ### columns: POS   chromosome    ancestry
   #            111   sample1_test    pop1
   #            124   sample1_test    pop1
   #
   #            111   sample1_true    pop1
   #            124   sample1_true    pop1
   #
   #            111   sample2_test    pop2
   #            124   sample2_test    pop2
   #
   #            111   sample2_true    pop1
   #            124   sample2_true    pop1

   spacing = {}
   for i, col_name in enumerate(ancestry_cols):
       spacing[col_name+'_test'] = 2*i - 2 * y_scale
       spacing[col_name+'_true'] = 2*i - 1 * y_scale

   # taks above example to something like:
   ###
   ### columns: POS  chromosome  ancestry
   #            111       0        pop1
   #            124       0        pop1
   #
   #            111       1        pop1
   #            124       1        pop1
   #
   #            111       2        pop2
   #            124       2        pop2
   #
   #            111       3        pop1
   #            124       3        pop1

   for col_name, spacing_val in spacing.items():
       melted.replace(col_name, spacing_val, inplace=True)

   plot = ggplot.ggplot(ggplot.aes(x='POS', y='chromosome', color='ancestry'), data=melted) \
       + ggplot.geom_point() \
       + ggplot.scale_y_continuous(labels=list(spacing.keys()), breaks=list(spacing.values())) \
       + ggplot.scale_color_manual(values=['#FF0000', '#0000FF', '#73008C']) \
       + ggplot.theme(plot_margin={'top':0.7, 'bottom':0.3}) ### TODO: this should depend on scale

   if image_filename is not None:
       plot.save(image_filename)
   else:
       plot.show()
