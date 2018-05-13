from sys import argv
import time

start_time = time.time()
print('starting split_homologous_chr')

in_filename = argv[1]
allele_filename = in_filename.replace('recode.vcf', 'homologous.txt')
homologous_filename = in_filename.replace('recode.vcf', 'homologous.vcf')

with open(in_filename, 'r') as in_file, open(allele_filename, 'w') as allele_file, open(homologous_filename, 'w') as homologous_file:
	header_found = False
	for line in in_file:
		if line[0] != '#': ### skip comment lines
			if not header_found: ### look for first line that isn't a comment
				header_found = True ### mark header as found
				line = line.split() ### header: NA001 NA002
				line[9:] = sum([[col + '_1', col + '_2'] for col in line[9:]], [])
				allele_file.write('\t'.join(line) + '\n') ### header: NA001_1 NA000_2
				homologous_file.write('\t'.join(line) + '\n') ### header: NA001_1 NA000_2
			else:
				line = line.split()	###  0|0   0|1
				line[9:] = sum([col.split('|') for col in line[9:]], [])
				allele_file.write('\t'.join(line) + '\n') ### 0 0 0 1
				line[9:] = ['{}|{}'.format(col, col) for col in line[9:]]
				homologous_file.write('\t'.join(line) + '\n') ### 0|0 0|0 0|0 1|1

elapsed_time = time.time() - start_time
print('finished')
print('time elapsed: {}'.format(elapsed_time))
