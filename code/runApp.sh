#!/bin/sh

code=`dirname $0`

usage="$0  <sampleID> <r1.fq> <r2.fq> <hg19.fasta> <manifest.txt>"

if [ $# -lt 5 ];then
	echo $usage
	echo "Code repo is $code"
	exit
fi


if [ ! -d $code ];then
	echo No code repo found for scripts
	echo exiting
	exit
fi


#inputs
#sample ID
sid="$1"
#paired FASTQ
fastqs="$2 $3 "
#ref genome in fasta
REFGENOME="$4"
#manifest file
manifest="$5"

outdir="$sid.outdir"
variantdir="$outdir/variants"

#created by App
GENOMEIDX="hg19.fasta.nix"
ampliconbed="$outdir/manifest2bed/ampliconNovoalign.bed"
amplicon_noprimer="$outdir/manifest2bed/ampliconNoPrimerRegion.bed"
variants="$outdir/variants.vcf.gz"
bam="$outdir/alignments.bam"

mkdir -p $outdir $variantdir

echo `date` Indexing the genome
novoindex $GENOMEIDX $REFGENOME 

#ManifestConversion
perl $code/manifestToAmplicon.pl $manifest -outdir $outdir/manifest2bed/


#Alignment
echo `date` Align and Sorting
novoalign -d $GENOMEIDX -f $fastqs --amplicons $ampliconbed -oSAM  2>novoalign.std_err.log |  samtools view -uS - > $outdir/novoalign.raw.sam
novosort $outdir/novoalign.raw.sam -o $bam -i $outdir/novoalign.raw.sam  2>log.txt 
rm  $outdir/novoalign.raw.sam

#Coverage Analysis
intersectBed -abam $bam  -b  $amplicon_noprimer  -f 0.99999 -bed -wa -wb |sort -k 10,10 | groupBy -g 10 -c 1 -o count > $outdir/novo_coverage.tsv
#Plot graph
#Rscript coverage.rscript

#Variants calling
mkdir -p $variantdir
perl $code/callAmpliconVariants.pl -percent 0.99999 -cov 50 $amplicon_noprimer -bam $outdir/novoalign.bam -genome $REFGENOME -out $variantdir/ | parallel -j+0
bcftools merge --force-samples -o $variants -Oz $variantdir/*gz
tabix -p vcf $variants









