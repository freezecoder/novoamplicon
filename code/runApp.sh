#!/bin/sh

LIM=""
LIM="-#25k"
export SHELL="/bin/sh"
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
echo "################### `date` Pipeline started on `hostname` ####################"
echo "Output dir: $outdir"
echo "Input: SID: $sid, Genome=$REFGENOME, fastqs=$fastqs, manifest=$manifest"

if [ ! -e $GENOMEIDX ]; then
	echo `date` Indexing the genome
	novoindex $GENOMEIDX $REFGENOME 
fi

echo Amplicon Bed file: $ampliconbed 
echo Converting:  perl $code/manifestToAmplicon.pl $manifest -outdir $outdir/manifest2bed/
    perl $code/manifestToAmplicon.pl $manifest -outdir $outdir/manifest2bed/

##### Test code for chr22 ####
# checks that amplicon and reference genome have same IDs
#ensure ref genome and amplicon bed have same chr list
cp $ampliconbed original.bed
cut -f 1  $REFGENOME.fai |sed 's/^/\^/' > chrlist.txt
grep -f chrlist.txt $ampliconbed > tmpfile && mv tmpfile $ampliconbed
grep -f chrlist.txt $amplicon_noprimer > tmpfile &&  mv tmpfile $amplicon_noprimer  


#check for existence of amplicon BED file
if [ ! -e $ampliconbed ];then
	echo "Error: No Amplicon bed file found in $ampliconbed"
	echo "exiting"
	exit 
fi

echo "################################################"
echo "################################################"
echo "######  novoalign $LIM -d $GENOMEIDX -f $fastqs --amplicons $ampliconbed -oSAM -rR -k  "
echo "################################################"
#Alignment and Bam-sorting
if [ ! -e $bam ];then
	echo `date` Aliging reads with Novoalign $novoversion on `hostname`
	novoalign $LIM -d $GENOMEIDX -f $fastqs --amplicons $ampliconbed -oSAM -rR -k  2> $outdir/novoalign_log.txt |  samtools view -uS - > $outdir/novoalign.raw.bam
	echo `date` Sorting alignments 
	novosort $outdir/novoalign.raw.bam -o $bam -i  2>$outdir/novosort_log.txt 
	rm  $outdir/novoalign.raw.bam
fi

doneflag=`grep "Done at" $outdir/novoalign_log.txt  |wc -l`

if [ $doneflag -lt 1 ];then
    echo "`date` Alignment Failed on `hostname`  see $outdir/novoalign_log.txt for more details"
    exit    
fi
echo "###############" `date` Novoalign done "###################:"



#Coverage Analysis

if [ ! -e "$outdir/novo_coverage.tsv" ];then
echo "################################################"
echo "################################################"
echo "`date` Calculating Coverage  for $outdir"
namecol=16
echo -e "TargetAmplicon\tCount" > $outdir/novo_coverage.tsv
intersectBed -abam $bam  -b  $amplicon_noprimer  -f 0.99999 -bed -wa -wb |sort -k $namecol,$namecol | groupBy -g $namecol -c 1 -o count >> $outdir/novo_coverage.tsv
Rscript $code/coveragePlot.R $outdir/novo_coverage.tsv $outdir/amplicon_coverage.png
#fi
#Plot graph
#Rscript coverage.rscript
fi

echo "################################################"
echo "################################################"
echo `date` Pipeline calling variants
echo "################################################"
echo "################################################"

echo Temp is $outdir/tmp ,our is $variantdir
#Variants calling
mkdir -p $variantdir 
mkdir -p "$outdir/tmp"
perl $code/callAmpliconVariants.pl \
 -percent 0.99999 \
 -cov 10  \
 -temp  $outdir/tmp \
 -bam $bam  \
 -genome $REFGENOME \
 -out $variantdir \
 -bed  $amplicon_noprimer 
ls $variantdir/*.gz | parallel -j+0 tabix -f -p vcf  {}
echo "`date` Variant calling complete"
echo "################################################"
echo "################################################"

echo `date` Pipeline merging variant files 
#merge and index with Tabix
bcftools merge --force-samples -o $variants -Oz $variantdir/*gz
tabix -f  -p vcf $variants
bcftools stats $variants > $variants.stats.txt


#Summary files
echo "Generating stat files"
perl -lane 's/:\s+/,/g  and s/#\s+//  and print  if /Paired Reads/ .. /CPU Time/' $outdir/novoalign_log.txt > $outdir/alignment_summary.csv
perl -lane 's/#\s+|Done.+//g and s/\t/,/g and print  if /#\s+Amplicon/ .. /Done/' $outdir/novoalign_log.txt > $outdir/amplicon_novoalign_count.csv
grep "^SN"  $variants.stats.txt |sed 's/samples/Amplicons/' |cut -f 3-100 | sed 's/:/,/' > $outdir/variantsummary.csv
grep "^TSTV"  $variants.stats.txt | perl -lane 'shift @F;  foreach(@F){$c+=$_;$e++};printf "TiTv,%.3f\n",$c/$e' >>  $outdir/variantsummary.csv

echo `date` Pipeline ended  on `hostname`

echo "Cleaning up index "
if [ -e $GENOMEIDX ];then
    rm $GENOMEIDX
fi
echo "################### `date` Pipeline Done ####################"






