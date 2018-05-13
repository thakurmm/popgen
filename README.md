#### Dependent Programs ####
Need to install [admixture](https://www.genetics.ucla.edu/software/admixture/), [plink2](https://www.cog-genomics.org/plink2) and [vcftools](https://vcftools.github.io/). Add these programs to your PATH.

#### Prepare data ####
Data downloaded from the [1000 Genomes Project](http://www.internationalgenome.org/data) - specifically, this ftp directory: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

1. Extract sample IDs.
	- Run the python script `SampleIDs/write_sampleIDs.py ASW CEU YRI` from this directory
		- Input data:
			- In SampleIDs: `pops_for_sample_IDs.tsv`
		- Output (in SampleIDs folder):
			- `ASW_CEU_YRI_IDs.txt`
			- `ASW_Sample_IDs_haploid.txt`
			- `CEU_Sample_IDs_haploid.txt`
			- `YRI_Sample_IDs_haploid.txt`

2. Select out populations of interest: ASW, CEU, and YRI.
	- Use the bash script (in bin folder): `Select_ASW_CEU_YRI_Seqs.sh`
		- 2018 update: Instead, use `select_populations.sh` from `bin/bashScripts/new/`. This does the same thing but is part of an effort to make scripts more flexible with regards to which populations they are looking at.
		- Calls:
			- vcftools
	- Input data:
		- In vcf folder: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
		- In SampleIDs: `ASW_CEU_YRI_IDs.txt`
	- Output: `chr22.phase3.ASW_CEU_YRI.SNPs.recode.vcf`, `chr22.phase3.ASW_CEU_YRI.SNPs.log`

3. Split the phased chromosomes into separate chromosomes
	- Use the bash script (in bin folder): `clean_chr_data_for_local_ancestry_split.sh`
	- 2018 update: Like with the previous step, use `bin/bashScripts/clean_chr_data_for_local_ancestry_split_new.sh` instead.
		- Calls:
			- Python script (in bin folder): `split_homologous_chr.py`
			- admixture
			- plink2
			- vcftools
	- Input data (in vcf folder): `chr22.phase3.ASW_CEU_YRI.SNPs.recode.vcf`
	- Output: 
		- `chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.Q`
		- `chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.P`
		- `chr22.phase3.ASW_CEU_YRI.SNPs.homologous.vcf`
		- `chr22.ASW_CEU_YRI.SNPs.homologous.header.vcf`
		- `chr22.ASW_CEU_YRI.SNPs.homologous.25k.bed`
		- `chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt`

4. Run `create_admixed_chromosomes` to parse the output of `split_homologous_chr.py` (which is called as part of the previous step) and create an admixed training set.
	- Input data:
		- In `SampleIDs/` directory: two lists of sample IDs, one for each training population.
		- In `data/` directory: source of genetic data, e.g. `chr22.phase3.ASW_CEU_YRI.SNPs.homologous.txt`
	- Output data:
		- In `data/admixed/` directory:
			- admixed training set, in two forms:
				- `CEU_YRI_admixed_40admixed_100pure.vcf`
				- `CEU_YRI_admixed_40admixed_100pure_homologous.vcf`
				- These have the same column names, but when the first vcf has `0|0`, the homologous vcf has `1|1`, and when the first vcf has `1|1`, the homologous vcf has `1`.
			- True ancestry proportions for the admixed training set
				- e.g. `CEU_YRI_admixed_40admixed_100pure_proportions.txt`
				- **not** used by STRUCTUREpainter (we only have this information here because of our artificially created training set). Used in `plot_admix_results.py` script to see how well ADMIXTURE does in inferring ancestry proportions.

5. Run `generate_test_set.py`.
	- Input data:
		- In `SampleIDs/` directory: two lists of sample IDs, one for each of the populations that we want to admix together.
		- In `data/` directory: genetic information (formatted like a VCF) of the individuals referenced by the sample IDs.
	- Output data:
		- In `test_input/` directory (test input to STRUCTUREpainter):
			- csv containing ancestry information
			- csv containing genetic information
	- Right now the parameters are hard-coded, you can change them at the bottom of the file. Everything under `if __name__ == "__main__"` is only run when you directly call `python generate_test_set.py`, but not when you import *from* this file. The same thing is done in the other Python scripts, so that they can import things from each other without side effects, i.e. extra code being run.

6. After you've followed these steps, `test_ancestry.sh` will run STRUCTUREpainter with appropriate parameters.

7. `evaluate_inferences.py` will plot the results of STRUCTUREpainter against the true ancestry values for the test set of chromosomes. See `evaluate_inferences.sh` for sample parameters. The Python script only requires two files, containing true and inferred ancestry, as well as a path for the plot that will be produced and saved.

#### 2018 TODOs ####
- Right now the transition matrix is based on population-level ancestry proportions, so the transition matrix is the same for each individual. It would probably be more accurate if the transition matrix is calculated for each individual based on that individual's estimated ancestry proportions (estimated using ADMIXTURE). Currently, STRUCTUREpainter only expects the results of ADMIXTURE on the training set, not on the evaluation set ("canvas" set?), which means that this is not a trivial change.

#### Local Ancestry Inference ####
- Run the local ancestry method which selects ancestry informative SNPs, estimates the transition + emission matrices, and iterates through the Hidden Markov Model (HMM).
- Optional input parameters:
	- `--kmer` (int, default = 5): window size to use
	- `--num_windows` (int, default = same as kmer): number of ways to split each chromosome into windows. If not specified, the value for `kmer` will be used. If specified, should be less than or equal to `kmer`.
		- Another way to think about what this parameter means: for each SNP, how many calls do you want to make and then average together to determine ancestry? `num_windows` exactly determines this.
		- Using values other than `kmer` is mostly untested right now.
		- Note that approximately, runtime should scale linearly with each of `kmer` and `num_windows`.
	- `--seed` (int, default = 0): random seed to use. Since the default value is an int (not Python's `None`), results will be exactly the same from run to run. Set to `None` for irreproducible results (varying randomly from run to run), or to a different integer value to see a different set of results that will be the same from run to run.
	- `--num_test` (int, default = -1): Number of test chromosomes (from the supplied test set, detailed below) to paint. If -1 (the default value), paints all test chromosomes from the test set.
- Required input parameters:
	- `--reference_filename` (string): path to file that contains training set
	- `--all_admix_filename` (string): path to overall admix filename
	- `--chrom_admix_filename` (string): path to admix filename for this chromosome
	- `--test_filename` (string): path to file that contains chromosomes for which ancestry should be inferred


#### Data Sources ####
- The file popgen/pops_for_sample_IDs.tsv can be downloaded from http://www.internationalgenome.org/data-portal/sample . Select the populations of interest (YRI, CEU and ASW) from the checkboxes, and then select the option to "Download the list" at the top of the page.
	- Rename the file as pops_for_sample_IDs.tsv
- The input data files for the vcftools command can be downloaded from http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ or ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
    - As mentioned on this page, the VCF files produced by the final phase of the 1000 Genomes Project (phase 3) are phased (http://www.internationalgenome.org/faq/are-1000-genomes-variant-calls-phased/)


#### Tools used ####
- vcftools
    - git clone https://github.com/vcftools/vcftools.git
    - cd vcftools
    - ./autogen.sh
    - ./configure
    - make
    - make install
