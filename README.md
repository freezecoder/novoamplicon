Amplicon Analysis for Illumina TruSeq
=======================================



Requires Novoalign v3.02.10 or higher 

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


