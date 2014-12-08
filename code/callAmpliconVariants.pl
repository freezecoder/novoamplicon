#!/usr/bin/perl

my $genome="";
my $bam="";
my $bcftools="bcftools";
my $samtools="samtools";
use Getopt::Long;
use File::Temp qw /tempfile/;

my $expectedcov=1000;
my $outdir="outvcf";
my $percent="0.8";
&GetOptions(
	'percent=s'=>\$percent,
	'bam=s'=>\$bam,
	'samtools=s'=>\$samtools,
	'bcftools=s'=>\$bcftools,
	'cov=s'=>\$expectedcov,
	'genome=s'=>\$genome,	
	'outdir=s'=>\$outdir
);
my $version=0.1;

my $usage="#Amplicon variant caller $version\n$0  -percent 0.9  -cov 1000 regions.bed  -bam file.bam -genome genome.fasta -out <output folder>  | parallel -j+0\n";
$usage.="Requires\n\tsamtools[$samtools] and\n\tbcftools[$bcftools]\n";
if (scalar @ARGV==0) {
	print $usage;
	exit 1;
}

my $bed=shift;

print "### Input BED regions: $bed\n";

mkdir("./bedA");
mkdir("$outdir") unless -e $outdir;
open(IN,$bed) or die "$!";

while(<IN>) {
	next if /^$|^#/;
	my @F=split;
	my($chr,$start,$end,$name,$score,$strand)=@F;
	#print "$chr\t$start\t$end\n";
	#my($fh,$fn)=tempfile(DIR=>"./bedA",SUFFIX=>".bed");
	open(BD,"+>bedA/$name.bed") or die "!";
	#print "$fn\t$genome\t$bam\n";
	print BD  "$_";
	close BD ;
	
	my $vcf="$outdir/$name.ampl.vcf";
	interCall("bedA/$name.bed",$bam,$genome,"$vcf",$name);
}
close IN;

sub interCall {

	my $bed=shift;
	my $bam=shift;
	my $genome=shift;
	my $out=shift;
	my $name=shift;
	#my $cmd="intersectBed  -abam $bam -b $bed  -wa -f $percent | $samtools mpileup -d 99999 -B -A  -uf $genome - | $bcftools call -vm -Oz ";
	my $cmd="intersectBed  -abam $bam -b $bed  -wa -f $percent | $samtools mpileup -d 99999 -B -A  -uf $genome - | $bcftools call -vm -Oz -o $out.gz";
	#$cmd.="| $bcftools filter -g 10  -s lowDPQ -i 'DP>=$expectedcov | QUAL>30' -Oz -o $out.gz";
	#my $cmd2="intersectBed  -abam $bam -b $bed  -wa -f 0.8 > Avcf/$name.bam";

	print "$cmd\n";
}
