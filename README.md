Amplicon Analysis for Illumina TruSeq
=======================================





Dependencies
----------------------------------------

Requires Novoalign v3.02.10 or higher. Please download from www.novocraft.com and request a free license key if you're not a licensed user as yet.

* Bedtools (2.15 or higher) https://github.com/arq5x/bedtools2
* Samtools (0.1.18 or higher)
* bcftools (latest) See http://samtools.github.io/bcftools
* tabix /bgzip  

All the dependencies must be in the UNIX $PATH

Protocol
------------------------------------

Run from the UNIX command line.


Works for paired-end reads only
Requires amplicon definitions

Run as the following

```sh
bash novoamplicon/code/runApp.sh 1111 NA12878-AFP1_S17_L001_R1_001.fastq.gz NA12878-AFP1_S17_L001_R2_001.fastq.gz genome.fa novoamplicon/t/truseq_amplicon_cancer_panel_manifest_afp1_pn15032433_b.txt
```

```sh
bash novoamplicon/code/runApp.sh  <sid> <read1.fastq> <read2.fastq> <reference.fasta> <manifest_file> 

```

Output
----------
An output directory named by sample ID


