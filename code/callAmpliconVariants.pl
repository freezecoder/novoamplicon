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
my $tmp="tmpdir";
my $bed;
&GetOptions(
	'percent=s'=>\$percent,
    'temp=s'=>\$tmp,
	'bam=s'=>\$bam,
     'bed=s'=>\$bed,
	'samtools=s'=>\$samtools,
	'bcftools=s'=>\$bcftools,
	'cov=s'=>\$expectedcov,
	'genome=s'=>\$genome,	
	'outdir=s'=>\$outdir
);
my $version=0.1;

my $usage="#Amplicon variant caller $version\n$0  -percent 0.9  -cov 1000 regions.bed  -bam file.bam -genome genome.fasta -out <output folder>  | parallel -j+0\n";
$usage.="Requires\n\tsamtools[$samtools] and\n\tbcftools[$bcftools]\n";

my $lines=0;
print "$0 ### Input BED regions: $bed\n";
print "$0 temp  dir is $tmp, genome is $genome, bam=$bam\n";
open(IN,$bed) or die "No $bed $!";
open(SCRIPT,"+>$tmp/varcall.sh") or die "$!";
while(<IN>) {
	next if /^$|^#/;
	my @F=split;
	my($chr,$start,$end,$name,$score,$strand)=@F;
	#print "$chr\t$start\t$end\n";
	#my($fh,$fn)=tempfile(DIR=>"./bedA",SUFFIX=>".bed");
	open(BD,"+>$tmp/$name.bed") or die "!";
	#print "$fn\t$genome\t$bam\n";
	print BD  "$_";
	close BD ;
	
	my $vcf="$outdir/$name.ampl.vcf";
	my $cmd;
    $cmd=interCall("$tmp/$name.bed",$bam,$genome,$vcf,$name);
    print SCRIPT $cmd,"\n";
    $lines++ if $cmd;
    system("$cmd");

}
close IN;
close SCRIPT;

print "$0 Done\n";
exit;

if ($lines>0) {
    system("cat $script| parallel -j+0 ");
}


#system("rm -fr $tmp") if -d $tmp;

sub interCall {

	my $bed=shift;
	my $bam=shift;
	my $genome=shift;
	my $out=shift;
	my $name=shift;
	my $cmd="intersectBed  -abam $bam -b $bed  -wa -f $percent | $samtools mpileup -d 99999 -B -A  -uf $genome - 2> $tmp/mpileup.err | $bcftools call -vm -Oz -o $out.gz";
    #$cmd.="| $bcftools filter -g 10  -s lowDPQ -i 'DP>=$expectedcov | QUAL>30' -Oz -o $out.gz";
	#my $cmd2="intersectBed  -abam $bam -b $bed  -wa -f 0.8 > Avcf/$name.bam";

	return($cmd);

}
