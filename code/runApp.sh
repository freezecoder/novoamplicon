#!/bin/sh

LIM=""
LIM="-#25k"

code=`dirname $0`
usage="$0  <sampleID> <r1.fq> <r2.fq> <hg19.fasta> <manifest.txt>"
if [ $# -lt 5 ];then
	echo $usage
	echo "Code repo is $code"
	echo "Debug:Lim=$LIM,Host=`hostname`,Version=`novoalign version`"
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

outdir="$sid"
variantdir="$outdir/variants"

#created by App
GENOMEIDX="$outdir/genome.nix"
ampliconbed="$outdir/manifest2bed/ampliconNovoalign.bed"
amplicon_noprimer="$outdir/manifest2bed/ampliconNoPrimerRegion.bed"
variants="$outdir/variants.vcf.gz"
bam="$outdir/alignments.bam"


novoversion=`novoalign version`
mkdir -p $outdir $variantdir
echo "Novoalign $novoversion"
echo `date` Pipeline started on `hostname`
echo "Input: SID: $sid, Genome=$REFGENOME, fastqs=$fastqs, manifest=$manifest"

if [ ! -e $GENOMEIDX ]; then
	echo `date` Indexing the genome
	novoindex $GENOMEIDX $REFGENOME 
fi


if [ ! -e $ampliconbed ];then
#ManifestConversion
perl $code/manifestToAmplicon.pl $manifest -outdir $outdir/manifest2bed/
fi

#check for existence of amplicon BED file
if [ ! -e $ampliconbed ];then
	echo "Error: No Amplicon bed file found in $ampliconbed"
	echo "exiting"
	exit 
fi

#Alignment and Bam-sorting
if [ ! -e $bam ];then
	echo `date` Aliging reads with Novoalign $novoversion on `hostname`
	novoalign $LIM -d $GENOMEIDX -f $fastqs --amplicons $ampliconbed -oSAM -rR -k  2> $outdir/novoalign_log.txt |  samtools view -uS - > $outdir/novoalign.raw.bam
	echo `date` Sorting alignments 
	novosort $outdir/novoalign.raw.bam -o $bam -i  2>$outdir/novosort_log.txt 
	rm  $outdir/novoalign.raw.bam
fi

rm $GENOMEIDX


#Coverage Analysis

if [ ! -e "$outdir/novo_coverage.tsv" ];then
echo `date` Calculating Coverage
namecol=16
echo -e "TargetAmplicon\tCount" > $outdir/novo_coverage.tsv
intersectBed -abam $bam  -b  $amplicon_noprimer  -f 0.99999 -bed -wa -wb |sort -k $namecol,$namecol | groupBy -g $namecol -c 1 -o count >> $outdir/novo_coverage.tsv
#Plot graph
#Rscript coverage.rscript
fi

echo `date` Pipeline calling variants
#Variants calling
mkdir -p $variantdir
perl $code/callAmpliconVariants.pl \
 -percent 0.99999 \
 -cov 10  \
  $amplicon_noprimer \
 -bam $bam  \
 -genome $REFGENOME \
 -out $variantdir/ | parallel -j+0
#index files
ls $variantdir/*.gz | parallel -j+0 tabix -p vcf  {}

echo `date` Pipeline merging variant files 
#merge and index with Tabix
bcftools merge --force-samples -o $variants -Oz $variantdir/*gz
tabix -p vcf $variants

echo `date` Pipeline ended  on `hostname`









