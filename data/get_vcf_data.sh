#!/usr/bin/env bash
# This script is expected to be run from the popgen folder.
# This script will download very large amounts of data and will take a LONG time to run.

URL="http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3"
mkdir -p data/vcf
cd data/vcf
for i in {1..22}
do
    wget $URL/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done
