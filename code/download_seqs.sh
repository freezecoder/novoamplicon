#!/bin/sh


dir="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/"

mkdir -p $dir

cd $dir 

echo Downloading chr10
wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr10.fa.gz  -O chr10.fa.gz
echo Downloading 22
wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz  -O chr22.fa.gz

zcat chr10.fa.gz chr22.fa.gz > genome.fa

echo genome is made in $dir/genome.fa
ls -l $dir/genome.fa




