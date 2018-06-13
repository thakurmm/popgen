#!/usr/bin/env python3

'''
parses the output of `split_homologous_chr.py` and creates a number of admixed genomes from the two given populations
- run the script:
    mkdir -p data/admixed
    python create_admixed_chromosomes.py --num_admixed_chromosomes 10 --num_anchor 300 --source_pops CEU YRI

(outputs two files in data/admixed/ directory: a .vcf file containing admixed chromosomes, and a proportions.txt file
    containing the ancestry proportions of each of the generated chromosome, for comparison with the output
    of admixture program)

next steps:
- use plink to generate .bed file:
    mkdir -p data/admixture/new/bed
    plink --vcf data/admixed/CEU_YRI_admixed_n\=200.vcf --make-bed --out data/admixture/new/bed/CEU_YRI_admixed_n\=200

- run admixture on .bed file:
    admixture --cv data/admixture/new/bed/CEU_YRI_admixed_n\=100.bed 2
    mkdir -p data/admixture/new/admix
    mv CEU_YRI_admixed_n\=100.2.* data/admixture/new/admix
'''
import sys
import numpy as np
import pandas

import argparse
from collections import defaultdict
#@ToDo Commented the two imports below as they are unused
#    import itertools
#    import pdb


def create_test_chromosome(chr_strands, pops, ref_IDs):
    ### sample from SNPs
    recombination_spots = (
        chr_strands['POS']
        .sample(NUM_RECOMBINATIONS, replace=False)
        .sort_values()
        .reset_index(drop=True)
    )

    stop = pandas.Series(range(NUM_RECOMBINATIONS)).map(lambda x: np.where(chr_strands['POS'] == recombination_spots[x])[0][0]+1)
    start = pandas.concat([pandas.Series(0), stop], ignore_index=True)
    stop = pandas.concat([stop, pandas.Series(len(chr_strands))], ignore_index=True)

    test_chr = np.zeros(len(chr_strands), dtype=np.uint8)
    true_chr = np.array([0] * len(chr_strands), dtype=np.uint8)
    for j in range(2):
        for i in range(j, len(start), 2):
            test_chr[start[i] : stop[i]] = chr_strands[ref_IDs[j]][start[i] : stop[i]]
            true_chr[start[i] : stop[i]] = j

    ancestry_proportions = true_chr.sum() / true_chr.shape[0]
    return test_chr, ancestry_proportions # proportion of SNPs selected from second ancestry


def create_pure_chromosomes_from_all_pops(chr_strands, n, pop1_ids, pop2_ids):
    return [create_pure_chromosomes_from_one_pop(chr_strands, int(n/2.0), popIDs) for popIDs in [pop1_ids, pop2_ids]]


def create_pure_chromosomes_from_one_pop(chr_strands, n, popIDs):
    assert n <= len(popIDs)

    IDs_to_use = np.random.choice(popIDs, size=n, replace=False)
    chroms = [chr_strands[popID] for popID in IDs_to_use]
    return [chroms, IDs_to_use]


def create_test_chromosome_set(chr_strands, pops_dict, num_chromosomes):
    ### pops_dict should map population name to list of haploid IDs corresponding to that population

    ID_pairs = []
    ### copy dict values (lists) since they will be modified; two copies of each haplotype, so ignore one
    ### remove IDs not in chr_strands
    pops_ordered = ['pop1', 'pop2']
    for i in range(num_admixed_chromosomes):
        if i % 20 == 0:
            print('on chromosome {}'.format(i))
        ID_pairs.append([pops_dict[x][np.random.randint(0, len(pops_dict[x]))] for x in pops_ordered]) ### get a random ID from each pop

    test_chromosomes_and_proportions = [create_test_chromosome(chr_strands, pops_ordered, pair) for pair in ID_pairs]
    test_chromosomes, second_ancestry_proportions = zip(*test_chromosomes_and_proportions)
    ancestry_proportions = [[1 - proportion, proportion] for proportion in second_ancestry_proportions]

    return test_chromosomes, ancestry_proportions, ID_pairs


def get_pop_IDs(filename):
    ret = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip('\n')
            if line != '':
                ret.append(line)
    return ret


if __name__ == '__main__':
    NUM_RECOMBINATIONS = 8 ### when creating offspring chromosomes

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--seed', help='random seed to use.', type=int, default=0)
    parser.add_argument('--num_admixed_chromosomes', help='number of admixed chromosomes to create.', type=int)
    parser.add_argument('--num_anchor', help='TOTAL number of pure chrosomes to create', type=int)
    parser.add_argument('--source_pops',
        help='population codes (e.g. "YRI CEU") for the two populations that will be used in creating the simulated test chromosomes',
        nargs=2, required=True)

    args = parser.parse_args()

    seed = args.seed
    num_admixed_chromosomes = args.num_admixed_chromosomes
    num_anchor = args.num_anchor
    sourcepop1, sourcepop2 = args.source_pops

    np.random.seed(seed)

    # homologous_filename = '/home/greg/School/popgen/data/chr22.phase3.{}.SNPs.homologous.txt'.format('ASW_CEU_YRI') ### TODO: change this
    populations = 'ASW_CEU_YRI'
    chr_number=22
    #@ToDo the chr22.phase3 .... file should be renamed to just chr.phase3.... , because this is in the Chr21, Chr22 folder and so on
    allele_filename = 'pops_data/{}_Data/Chr{}/tmp/chr22.phase3.{}.SNPs.allele.vcf'.format(populations, chr_number, populations)

    #sourcepop1_ids_filename = '/home/greg/School/popgen/SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop1)
    #sourcepop2_ids_filename = '/home/greg/School/popgen/SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop2)
    sourcepop1_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop1)
    sourcepop2_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop2)

    #out_basename = '/home/greg/School/popgen/data/admixed/{}_{}_admixed_{}admixed_{}pure'.format(sourcepop1, sourcepop2, num_admixed_chromosomes, num_anchor)
    out_basename = 'pops_data/admixed/{}_{}_admixed_{}admixed_{}pure'.format(sourcepop1, sourcepop2, num_admixed_chromosomes, num_anchor)
    # out_filename = '/home/greg/School/popgen/data/admixed/{}_{}_admixed_n={}.vcf'.format(sourcepop1, sourcepop2, num_admixed_chromosomes)
    #out_filename = out_basename + '.vcf' ### 0|0 0|0 1|1 0|0
    #out_filename_homologous = out_basename + '_homologous.vcf' ### 0 0 1 0
    out_filename_homologous = out_basename + '_HOMOLOGOUS.vcf' ### 0|0 0|0 1|1 0|0
    out_filename_allele = out_basename + '_ALLELE.vcf' ### 0 0 1 0
    #proportions_out_filename_homologous = '/home/greg/School/popgen/data/admixed/{}_{}_admixed_{}admixed_{}pure_proportions.txt'.format(sourcepop1, sourcepop2, num_admixed_chromosomes, num_anchor)
    proportions_out_filename = 'pops_data/admixed/{}_{}_admixed_{}admixed_{}pure_proportions.txt'.format(sourcepop1, sourcepop2, num_admixed_chromosomes, num_anchor)

    sourcepop1_ids = get_pop_IDs(sourcepop1_ids_filename)
    sourcepop2_ids = get_pop_IDs(sourcepop2_ids_filename)

    #all_chr_strands = pandas.read_csv(homologous_filename, sep='\t', header=0, comment='#')
    #@ToDo the allele_filename is from the tmp/ folder, with no "hash" on the CHROM line. Can we instead use the file with the hash in the parent folder?
    all_chr_strands = pandas.read_csv(allele_filename, sep='\t', header=0, comment='#')
    all_chr_strands.drop_duplicates(inplace=True)

    sourcepop1_ids, sourcepop2_ids = map(lambda ids: list(filter(lambda ID: ID in all_chr_strands.columns, ids)), [sourcepop1_ids, sourcepop2_ids])

    # The compicated lambda command is doing what the two (for ... for) loops are doing
    # resultant_pop1_ids = []
    # for id in sourcepop1_ids:
    #     for chr_id in all_chr_strands.columns:
    #         if id == chr_id:
    #             resultant_pop1_ids.append(id)
    # sourcepop1_ids = resultant_pop1_ids
    #
    #
    # resultant_pop2_ids = []
    # for id in sourcepop2_ids:
    #     for chr_id in all_chr_strands.columns:
    #         if id == chr_id:
    #             resultant_pop2_ids.append(id)
    # sourcepop2_ids = resultant_pop2_ids

    chromosomes, ancestry_proportions, ID_pairs = create_test_chromosome_set(all_chr_strands, {'pop1' : sourcepop1_ids, 'pop2' : sourcepop2_ids}, num_admixed_chromosomes)
    ### ancestry proportions is a list of [pop1_proportion, pop2_proportion] for each test chromosome
    print('created admixed chromosomes; creating anchor chromosomes')

    pure = create_pure_chromosomes_from_all_pops(all_chr_strands, num_anchor, sourcepop1_ids, sourcepop2_ids)
    source1pure, source2pure = pure

    num_anchor_per_ancestry = int(num_anchor/2.0)
    ancestry_proportions.extend([[1,0]] * num_anchor_per_ancestry)
    ancestry_proportions.extend([[0,1]] * num_anchor_per_ancestry)

    print('created anchor chromosomes; dropping unneeded columns')

    all_chr_strands.drop(all_chr_strands.columns[9:], axis=1, inplace=True)
    to_write = all_chr_strands
    to_write_homologous = to_write.copy()

    print('finished dropping unneeded columns')

    print('building dataframe, starting with admixed chromosomes')

    used_names = defaultdict(int)
    for index, (chrom, id_pair) in enumerate(zip(chromosomes, ID_pairs)):
        if index % 20 == 0:
            print('on column: {}'.format(index))

        name = '-'.join(id_pair).replace('_', '-')
        idx = used_names[name]
        used_names[name] += 1
        name = '{}-{}'.format(name, idx)

        to_write[name] = np.fromiter(('{0}|{0}'.format(snp) for snp in chrom), dtype='<U3')
        to_write_homologous[name] = chrom

    if (num_anchor > 0):
        print('adding pure chromosomes to dataframe')
    else:
        print('not adding any pure chromosomes to dataframe')

    for pure_chroms, IDs in [source1pure, source2pure]:
        for column, (chrom, ID) in enumerate(zip(pure_chroms, IDs), len(to_write.columns)):
            if column % 20 == 0:
                print('on column: {}'.format(column))

            to_write[ID] = np.fromiter(('{0}|{0}'.format(snp) for snp in chrom), dtype='<U3')
            to_write_homologous[ID] = chrom

    #print('finished building dataframe. Writing to .vcf file: {}'.format(out_filename))
    print('finished building dataframe. Writing to .vcf file: {}'.format(out_filename_homologous))

    #with open(out_filename, 'w') as f:
    with open(out_filename_homologous, 'w') as f:
        f.write('##fileformat=VCFv4.2\n#')
        # f.write('##fileformat=VCFv4.2\n')
        # to_write.rename(columns={'CHROM' : '#CHROM'}, inplace=True)
        to_write.to_csv(f, sep='\t', index=False, mode='a')

    #with open(out_filename_homologous, 'w') as f_homologous:
    with open(out_filename_allele, 'w') as f_homologous:
        to_write_homologous.to_csv(f_homologous, sep='\t', index=False)

    print('finished writing to .vcf file. writing ancestry proportions: {}'.format(proportions_out_filename))

    with open(proportions_out_filename, 'w') as f:
        for ancestry_zero, ancestry_one in ancestry_proportions:
            f.write('{} {}\n'.format(ancestry_zero, ancestry_one))
