import pandas as pd

import pdb
from sys import argv

if len(argv) <= 2:
    print("This script should be run in the parent folder of SampleIDs folder, where the outputs will be created")
    print("Usage: {0} <path to samples file> <population_1> [<population 2> .....]".format(argv[0]))
    exit()

input_samples_filename=argv[1]
pops = argv[2:]

filename = '_'.join(pops) + '_IDs.txt'

# pd.read_csv() return a pandas DataFrame object
input_df = pd.read_csv(input_samples_filename, sep='\t')

# The line below will filter the items read above, so as to only include entries where Population code column
# is one of the codes entered as "pops" argument above.
# See learning_samples/pandas_dataframe.py to see how this works
filtered_df = input_df[input_df['Population code'].isin(pops)]
filtered_sample_name_df = filtered_df['Sample name']

filepath = 'SampleIDs/' + filename
print('saving IDs belonging to any of ( {0} ) to {1}'.format(', '.join(pops), filepath))
# write IDs (from the Sample name column) from all the input populations in the filepath defined above
# index=False means the DataFrame index column is omitted. Default index is True
filtered_sample_name_df.sort_values().to_csv(filepath, index=False)

for pop in pops:
    pop_df = input_df[input_df['Population code'] == pop]
    pop_sample_name_df = pop_df['Sample name']

    # For each Sample name, we append _1 and _2 to create the two haploids
    haploid_df = pd.concat([pop_sample_name_df.apply(lambda str: str + '_1'), pop_sample_name_df.apply(lambda str: str + '_2')]).sort_index()

    # Path to the output filename for each population
    haploid_filepath = 'SampleIDs/{0}_Sample_IDs_haploid.txt'.format(pop)
    print('saving haploid IDs for {0} to {1}'.format(pop,haploid_filepath))

    # write haploid IDs (id_1, id_2) from one population
    haploid_df.sort_values().to_csv(haploid_filepath, index=False)
