#Author: Abdul Malik bin Ahmad, Novocraft Technologies Sdn. Bhd., 2014
#Simple script to convert Illumina manifest file into bedfile for amplicon data; A gene list will also be generated
#Usage: perl manifestToAmplicon.pl <manifest.txt>
#
#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Data::Dumper;
use Getopt::Long;

my $outdir="outdir";

&GetOptions(
    'outdir=s'=>\$outdir
);

mkdir($outdir) unless -d $outdir;

my $usage = "\nUsage: perl convertManifestToBed.pl <IlluminaManifest.txt>\n";
print "Opening the manifest file --> ".$ARGV[0]."\n";	#Use the foreach @ARGV for multiple input, as for now, stick with one sample per run
open READMANIFEST, $ARGV[0] or die ( "\tCannot open the manifest file..\n$usage" );
my @manifestProbesHeaders;
my @manifestTargetsHeaders;
my $probesPrintCheckpoint = 0;
my $targetsPrintCheckpoint = 0;

my $manifestTargetId;
my $manifestChr;
my $manifestTargetIDChr;
my $manifestStrand;
my $manifestULSO;
my $manifestULSOcount; 
my $manifestDLSO;
my $manifestDLSOcount;
my %manifestHash;
my %manifestBedHash;

my $ampliconRegion = 0;
my $ampliconBp = 0;
my %ampliconTargetNumber;

#Extract needed column, based on Probes and Targets headers
print "\tChecking and extracting headers\n";
open WRITEPROBES,">tempProbes" or die ( "Cannot write to tempProbe" );
open WRITETARGETS,">tempTargets" or die ( "Cannot write to tempTargets" );

while ( <READMANIFEST> ) {
	if ( $probesPrintCheckpoint == 1 && $targetsPrintCheckpoint == 0 && $_ !~ /\[Targets\]/ ) { 	#Avoid printing [Targets] header
		print WRITEPROBES $_;
	}
	if ( $probesPrintCheckpoint == 1 && $targetsPrintCheckpoint == 1 ) { 			
		print WRITETARGETS $_;
	}
	if ( $_ =~ /\[Probes\]/ ) {				#Found the [Probes] head, next line is the columns name, change the checkpoint value to allow printing in above function block
		my $nextLine = <READMANIFEST>;
		$nextLine =~ s/^\s+|\s+$//g;		#Remove space(s)
		@manifestProbesHeaders = split(/\t/,$nextLine);	#Insert each column name into array
		$probesPrintCheckpoint = 1;			#The first function block will start printing the probes value
	} 
	if ( $_ =~ /\[Targets\]/ ) {
		my $nextLine = <READMANIFEST>;
		$nextLine =~ s/^\s+|\s+$//g;
		@manifestTargetsHeaders = split(/\t/,$nextLine);
		$targetsPrintCheckpoint = 1;
	}	
}
close ( READMANIFEST );
close ( WRITEPROBES );
close ( WRITETARGETS );

my $probesHeaderTargetID;
my $probesHeaderSpecies;
my $probesHeaderBuildID;
my $probesHeaderChromosome;
my $probesHeaderStartPosition;
my $probesHeaderEndPosition;
my $probesHeaderStrand;	
my $probesHeaderULSOSequence;	
my $probesHeaderDLSOSequence;

#Due to non-standardized column headers arrangement in manifest file, we're going to check which column is needed and extract that column only to a new file
print "\t\tFiltering the Probes headers column\n";
#print "\t\t\tOriginal header columns\n";
for (my $x = 0; $x < scalar( @manifestProbesHeaders ); $x++ ) {		#Go through each array which had one column name each, find the wanted header
	unless ( $manifestProbesHeaders[$x] eq '' ) {
		#print "\t\t\t".($x + 1)." for "."$manifestProbesHeaders[$x]"."\n";
	}
	if ( $manifestProbesHeaders[$x] eq "Target ID" ) {		#If the current array number match the string, assign $x + 1 for column as array is 0 based and columns start with 1
		$probesHeaderTargetID = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "Species" ) {
		$probesHeaderSpecies = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "Build ID" ) {
		$probesHeaderBuildID = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "Chromosome" ) {
		$probesHeaderChromosome = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "Start Position" ) {
		$probesHeaderStartPosition = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "End Position" ) {
		$probesHeaderEndPosition = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "Strand" || $manifestProbesHeaders[$x] eq "Submitted Target Region Strand" ) {	#Variety of headers name.
		$probesHeaderStrand = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "ULSO Sequence" ) {
		$probesHeaderULSOSequence = $x + 1;	
	}
	if ( $manifestProbesHeaders[$x] eq "DLSO Sequence" ) {
		$probesHeaderDLSOSequence = $x + 1;	
	}
}
#print "\n\t\t\tFiltered header columns
#\t\t\t$probesHeaderTargetID for Target ID
#\t\t\t$probesHeaderChromosome for Chromosome
#\t\t\t$probesHeaderStrand for Strand
#\t\t\t$probesHeaderULSOSequence for ULSO Sequence
#\t\t\t$probesHeaderDLSOSequence for DLSO Sequence\n\n";

#####################

my $targetsHeaderTargetA;
my $targetsHeaderTargetB;
my $targetsHeaderTargetNumber;
my $targetsHeaderChromosome;
my $targetsHeaderStartPosition;
my $targetsHeaderEndPosition;
my $targetsHeaderProbeStrand;
my $targetsHeaderSequence;
print "\t\tFiltering the Targets headers column\n";
#print "\t\t\tOriginal header columns\n";
for (my $x = 0; $x < scalar( @manifestTargetsHeaders ); $x++ ) {
	unless ( $manifestTargetsHeaders[$x] eq "" ) {
		#print "\t\t\t".($x + 1)." for "."$manifestTargetsHeaders[$x]"."\n";
	}
	if ( $manifestTargetsHeaders[$x] eq "TargetA" ) {
		$targetsHeaderTargetA = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "TargetB" ) {
		$targetsHeaderTargetB = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "Target Number" ) {
		$targetsHeaderTargetNumber = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "Chromosome" ) {
		$targetsHeaderChromosome = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "Start Position" ) {
		$targetsHeaderStartPosition = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "End Position" ) {
		$targetsHeaderEndPosition = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "Probe Strand" ) {
		$targetsHeaderProbeStrand = $x + 1;	
	}
	if ( $manifestTargetsHeaders[$x] eq "Sequence" ) {
		$targetsHeaderSequence = $x + 1;	
	}
}
#print "\n\t\t\tFiltered header columns
#\t\t\t$targetsHeaderTargetA for TargetA
#\t\t\t$targetsHeaderTargetB for TargetB
#\t\t\t$targetsHeaderChromosome for Chromosome
#\t\t\t$targetsHeaderStartPosition for Start Position
#\t\t\t$targetsHeaderEndPosition for End Position
#\t\t\t$targetsHeaderProbeStrand for Probe Strand
#\t\t\t$targetsHeaderSequence for Sequence\n\n";

#Write columns needed from temp files to filtered files
`cut -f$probesHeaderTargetID,$probesHeaderChromosome,$probesHeaderStrand,$probesHeaderULSOSequence,$probesHeaderDLSOSequence tempProbes > filteredProbes`;
`cut -f$targetsHeaderTargetA,$targetsHeaderTargetB,$targetsHeaderTargetNumber,$targetsHeaderChromosome,$targetsHeaderStartPosition,$targetsHeaderEndPosition,$targetsHeaderProbeStrand,$targetsHeaderSequence tempTargets > filteredTargets`;

#Open filteredProbes file
#Get main value except the start/end position as we'll get that value from [Targets] file
open READFILTEREDPROBES,"filteredProbes" or die ( "Cannot open filteredProbes file..\n" );
while ( <READFILTEREDPROBES> ) {
	#Save in the hash table for primers string (converted to length) with it's target id
	if ( $_ =~ /(.+)\t(chr..?)\t([\+|-])\t(\w+)\t(\w+)/ ) {
		#print "Target ID -> $1\nChr -> $2\nStrand -> $3\nULSO -> $4\nDLSO -> $5\n";
		$manifestTargetId = $1;
		$manifestChr = $2;
		$manifestStrand = $3;
		$manifestULSO = $4;
		$manifestDLSO = $5;
		$manifestULSOcount = length ($manifestULSO);	#This will get the length of the primer string
		$manifestDLSOcount = length ($manifestDLSO);
		
		$manifestHash{$manifestTargetId} = {
			targetId_key => $manifestTargetId,
			chr_key => "",
			startPos_key => "",
			endPos_key => "",
			strand_key => "",
			ULSOcount_key => $manifestULSOcount,	#reverse strand
			startThickPos_key => "",
			DLSOcount_key => $manifestDLSOcount,	#reverse strand
			endThickPos_key => "",
			seqLength_key => "",		
		}
	}
}
close ( READFILTEREDPROBES );
print "\tCalculating Start/End position and Thick Start/End position\n";
open READFILTEREDTARGETS,"filteredTargets" or die ( "Cannot open filteredTargets file..\n" );
while ( <READFILTEREDTARGETS> ) {
	$ampliconRegion = $ampliconRegion + 1;
	if ( $_ =~ /(.+)\t(.+)\t(\d+)\t(chr..?)\t(\d+)\t(\d+)\t([\+|-])\t(\w+)/ ) {
		my $tempTargetA = $1;
		my $tempTargetB = $2;
		my $tempTargetNumber = $3;
		my $tempChr = $4;
		my $tempStartPos = $5;
		my $tempEndPos = $6;
		my $tempTargetAChrPos = $1."_".$4."_".$5.$6;
		my $tempTargetBChrPos = $2."_".$4."_".$5.$6;
		my $tempProbeStrand = $7;
		my $tempSeq = $8;
		my $tempSeqLength = length( $tempSeq );
		
		#Extract target number value, referring to multiple location of similar amplicon sequence
		#Report later the count of how many amplicon had 1,2,3...n locations,
		if ( $ampliconTargetNumber{$tempTargetNumber}{counter_key} == 0 ) {
			$ampliconTargetNumber{$tempTargetNumber}{counter_key} = 1;
		} else {
			$ampliconTargetNumber{$tempTargetNumber}{counter_key} = $ampliconTargetNumber{$tempTargetNumber}{counter_key} + 1;
		}
		
		##Value for TargetA
		if ( $tempProbeStrand eq "+" ) {
			$manifestBedHash{$tempTargetAChrPos} = {
				targetId_key => $manifestHash{$tempTargetA}{targetId_key},
				chr_key =>$tempChr,
				startPos_key => $tempStartPos - 1,	
				endPos_key => $tempEndPos,
				strand_key => $tempProbeStrand,
				ULSOcount_key => $manifestHash{$tempTargetA}{ULSOcount_key},
				DLSOcount_key => $manifestHash{$tempTargetA}{DLSOcount_key},
				startThickPos_key => ( $tempStartPos + $manifestHash{$tempTargetA}{DLSOcount_key}) - 1,
				endThickPos_key => ( $tempEndPos - $manifestHash{$tempTargetA}{ULSOcount_key}) + 1,
				seqLength_key => $tempSeqLength,
			}
		} elsif ( $tempProbeStrand eq "-" ) {
			$manifestBedHash{$tempTargetAChrPos} = {
				targetId_key => $manifestHash{$tempTargetA}{targetId_key},
				chr_key =>$tempChr,
				startPos_key => $tempStartPos - 1,
				endPos_key => $tempEndPos,
				strand_key => $tempProbeStrand,
				ULSOcount_key => $manifestHash{$tempTargetA}{ULSOcount_key},
				startThickPos_key => ( $tempStartPos + $manifestHash{$tempTargetA}{ULSOcount_key}) - 1,
				DLSOcount_key => $manifestHash{$tempTargetA}{DLSOcount_key},
				endThickPos_key => ( $tempEndPos - $manifestHash{$tempTargetA}{DLSOcount_key}) + 1,	
				seqLength_key => $tempSeqLength,	
			}
		}
		##Value for TargetB
#		if ( $tempProbeStrand eq "+" ) {
#			$manifestBedHash{$tempTargetBChrPos} = {
#				targetId_key => $manifestHash{$tempTargetB}{targetId_key},
#				chr_key =>$tempChr,
#				startPos_key => $tempStartPos - 1,
#				endPos_key => $tempEndPos,
#				strand_key => $tempProbeStrand,
#				ULSOcount_key => $manifestHash{$tempTargetB}{ULSOcount_key},
#				DLSOcount_key => $manifestHash{$tempTargetB}{DLSOcount_key},
#				startThickPos_key => ( $tempStartPos + $manifestHash{$tempTargetB}{DLSOcount_key}) - 1,
#				endThickPos_key => ( $tempEndPos - $manifestHash{$tempTargetB}{ULSOcount_key}) + 1,
#				seqLength_key => $tempSeqLength,
#			}
#		} elsif ( $tempProbeStrand eq "-" ) {
#			$manifestBedHash{$tempTargetBChrPos} = {
#				targetId_key => $manifestHash{$tempTargetB}{targetId_key},
#				chr_key =>$tempChr,
#				startPos_key => $tempStartPos - 1,
#				endPos_key => $tempEndPos,
#				strand_key => $tempProbeStrand,
#				ULSOcount_key => $manifestHash{$tempTargetB}{ULSOcount_key},
#				startThickPos_key => ( $tempStartPos + $manifestHash{$tempTargetB}{ULSOcount_key}) - 1,
#				DLSOcount_key => $manifestHash{$tempTargetB}{DLSOcount_key},
#				endThickPos_key => ( $tempEndPos - $manifestHash{$tempTargetB}{DLSOcount_key}) + 1,	
#				seqLength_key => $tempSeqLength,	
#			}
#		}
	}
}
close (READFILTEREDTARGETS);

############
#Writing 
open WRITEDUMPER,">dumper";
print WRITEDUMPER Dumper(\%manifestBedHash);
close (WRITEDUMPER);
print "\tWriting to ampliconNovoalign.bed bedfile\n";
print "\tWriting to amplicon.bed bedfile\n";
print "\tWriting to ampliconNoPrimerRegion.bed bedfile\n";
print "\tWriting to ampliconGeneList.text file\n";
open WRITEAMPBED,">$outdir/amplicon.bed";
open WRITEAMPNOVOBED,">$outdir/ampliconNovoalign.bed";
open WRITENOPRIMERBED,">$outdir/ampliconNoPrimerRegion.bed";
open WRITEGENELIST,">$outdir/ampliconGeneList.txt";
open WRITERAWGENELIST,">$outdir/ampliconRawGeneList.txt";

foreach my $key (sort keys %manifestBedHash) {
my $geneName = "";
my $rawGeneName = "";
	if ( $key =~ /(\w+)\.?.+$/ ){
		$geneName = $1;
		$rawGeneName = $geneName;
		print WRITERAWGENELIST "$rawGeneName\n";
		if ( $geneName =~ /(.+)_\d.*$/ ) {
			my $tempGeneName = $1;
			if ( $tempGeneName =~ /(.+)_\d.*$/ ) {
				$geneName = $1;
				print WRITEGENELIST "$geneName\n";			
			} else {
				$geneName = $tempGeneName;
				print WRITEGENELIST "$geneName\n";
			}
		} else {
			print WRITEGENELIST "$geneName\n";
		}
	}
	if ( $manifestBedHash{$key}{chr_key} eq "" ) {
		print "Found empty key due to no match no probes TargetID-Chromosome\n";
	} else {
		print WRITEAMPNOVOBED $manifestBedHash{$key}{chr_key}."\t".
			$manifestBedHash{$key}{startPos_key}."\t".
			$manifestBedHash{$key}{endPos_key}."\t".
			$manifestBedHash{$key}{targetId_key}."\t0\t".
			$manifestBedHash{$key}{strand_key}."\t".		
			$manifestBedHash{$key}{startThickPos_key}."\t".
			$manifestBedHash{$key}{endThickPos_key}.
			"\n";
		
		print WRITEAMPBED $manifestBedHash{$key}{chr_key}."\t".
			$manifestBedHash{$key}{startPos_key}."\t".
			$manifestBedHash{$key}{endPos_key}."\t".
			$manifestBedHash{$key}{targetId_key}.
			"\n";
			
		print WRITENOPRIMERBED $manifestBedHash{$key}{chr_key}."\t".
			$manifestBedHash{$key}{startThickPos_key}."\t".
			$manifestBedHash{$key}{endThickPos_key}."\t".
			$manifestBedHash{$key}{targetId_key}."\t0\t".
			$manifestBedHash{$key}{strand_key}.
			"\n";
	}
}
close ( WRITEAMPNOVOBED );
close ( WRITEAMPBED );
close ( WRITENOPRIMERBED );
close ( WRITERAWGENELIST );
close ( WRITEGENELIST );
#Remove redundant gene list
`cat ampliconGeneList.txt|uniq - > uniqGeneList && mv uniqGeneList ampliconGeneList.txt`;

print "Conversion finish\n\n";
print "Amplicon Summary\n";
print "Manifest file\t: ".$ARGV[0]."\n";
print "Total Region\t: $ampliconRegion\n";
print "Target Count\t: ";
foreach ( sort keys %ampliconTargetNumber ) {
	print "$_ Region(s): ".$ampliconTargetNumber{$_}{counter_key}."\n\t\t  ";
}

print "\nFile Removal process..\n";
`rm -f tempProbes`;
`rm -f tempTargets`;
`rm -f filteredProbes`;
`rm -f filteredTargets`;
`rm -f dumper`;

print "Task end\n";
