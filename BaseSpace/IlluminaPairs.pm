#!/usr/bin/perl -w

package IlluminaPairs;

use Exporter;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename dirname/;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA         = qw(Exporter);
@EXPORT      = qw(findFiles makePairs);
@EXPORT_OK   = qw(findFiles makePairs);
%EXPORT_TAGS = ( DEFAULT => [qw(&findFiles makePairs)],
                 Both    => [qw(&findFiles &makePairs)]);



#Get all fastq gzipped files from an input directory
#Return array reference
sub findFiles {
my $indir=shift;
my $found = `find $indir -name "*gz"`;
chomp $found;
my @files=split("\n",$found);
return \@files;
}

#return object of sample pairs
#ignore the indexes with I[12]
sub makePairs {
	my $l=shift;
	my @list=@$l;
    my $count=0;
	my %h=();
	foreach $f (@list) {
		my $base=basename($f);
		my($samplename,$samplenumber,$lane,$read,$flowcellID)=split(/[_]+/,$base);
		next if $read=~/^I/;
			
		$flowcellID=~ s/\.fastq|\.gz|\.fq//g;
		#print "$f\t$samplename,$samplenumber,$lane,$read,$flowcellID\n";
		my $key="$samplename.$flowcellID";
		$h{$key}{$read}=$f;
		$count++;
	}
	printf "Illumina Object ::makePairs: $count files in %s samples\n", scalar(keys %h);
	#print Dumper(\%h);
	#printf "%s\n",join"\t",sort {$a cmp $b} keys %h;
	return \%h;
}


1;
