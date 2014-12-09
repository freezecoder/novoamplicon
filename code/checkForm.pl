#!/usr/bin/perl -w
#
#Novocraft BaseSpace pipeline
#Author: Zayed Albertyn  zayed@novocraft.com
#copyright Novocraft Technologies 2013

use FindBin qw($Bin);

use Getopt::Long;
use Data::Dumper;
use File::Find;
use File::Basename qw/basename dirname/;
use Time::HiRes qw( time );
use JSON;

use lib "$Bin/..";
use BaseSpace::IlluminaPairs;

my $b;
$b=basename($0);
chomp $b;
my $version = 1.1;

#Set the pipeline script
$script="/media/ephemeral/scratch/novoamplicon/code/runApp.sh";
#Reference genome file in fasta format e.g. hg19.fasta 
my $genomefasta="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
my $isDemo=0;
#Where all external programs live e.g. samtools , freebayes,etc
my $javabin = "/home/bin/jre1.7.0_40/bin";
my $appbase="/home/bin/export";
#Form options
my $align_only=0;#default is whole pipeline
my $trimadapt=0; #default is trimming
my $tempdir="";
my $readgroup="111";
my $samplename="sampleX";
my $libname="lib123";
my $alignopts="";
my $remove_duplicates=0;
my $callvars = 0;
#####
my $MEM="2g";
my $cleanup=0;
my $variantcount= 0;
my $genomeName="hschr22";
$genomeName="bacillus:";
my $is_debug=0;
my $appsession = "/data/input/AppSession.json";
my $outdir="/data/output/appresults/1234";
my $test_indir="/vagrant/appdata/reads/";
my $test_genome="/home/bin/novopipeline/fasta/bacillus.fasta";
my $indir="/data/inputs/";
&GetOptions(
	'indir=s'=>\$indir,
	'script=s'=>\$script,
    'norun!'=>\$norun,
       'testjson!'=>\$testjson,
    'appjson=s'=>\$appsession,
    'cleanup!'=>\$cleanup,
     'mem=s'=>\$MEM,
    'rmdup!'=>\$remove_duplicates,
    'align!'=>\$align_only,
    'genome|reference=s'=>\$genomefasta,
     'outdir=s'=>\$outdir,
     'apps=s'=>\$appbase,
      'tempdir=s'=>\$tempdir,
    'samreadgroup=s'=>\$readgroup,
    'samsample=s'=>\$samplename,
    'samlibrary=s'=>\$libname,
    'aligner_opts=s'=>\$alignopts, 
    'debug!'=>\$is_debug,
     'test_indir=s'=>\$test_indir,
     'test_genome|test_reference=s'=>\$test_genome,
  );
$outdir= "foo" if ($testjson);
my $usage="Novoalign pipeline: Alignment & variant calling\n$0 --outdir <outdir>  --reference <genome fasta> read1.fastq.gz read2.fastq.gz\n";
$usage.="-apps  string. Location of dependency programs [$appbase]\n";
$usage.="-align Boolean. Only do alignmnent  [$align_only]\n";
$usage.="-samreadgroup  string. SAM RG ID [$readgroup]\n";
$usage.="-reference String. Genome FASTA location [$genomefasta ]\n";
$usage.="-appjson App session JSON [$appsession]\n";
$usage.="-testjson Test JSON inputs and exit\n";
$usage.="\nZayed Albertyn 2014\n";
$usage.="\n$0  Version:$version\n";

unless (-e $appsession) {
    print "No app JSON $appsession\n";
    print $usage;
    exit 1;
}

print  "############$0 start ################\n";
print  "############ Illumina BaseSpace: Novoalign alignment pipeline################\n";
print STDERR "$0 Version $version\n";
print  "$0 Version $version\n";

my($FormData,$manifestfile,$smname,$sampleid,$projectid,$base,$is_amplicon);
#Read data from JSON from BaseSpace form
#Here is where we can fetch form variables such as no. reads to process, reference genome name,etc
#assign indir, sampleid and projectId
checkFormInputs();

if ($testjson) {
   $norun=1; 
}
#exit 0 if $testjson;

my $start = time();

print "$b.INFO Genome used in this pipeline is $genomeName\n";

chooseGenome();
    
$base="$outdir/";

my $aref= IlluminaPairs::findFiles($indir);
#print Dumper($aref);
my $obj = IlluminaPairs::makePairs($aref);

my %pairs=();
%pairs=%$obj;
my $paired=0;

if (scalar keys %pairs ==0) {
    print STDERR "No pairs or SE reads found in $indir,exiting\n";
    exit 500;
}

#mkdir ($outdir) unless -d $outdir;
#print Dumper($obj);
my $index="genome.nix";
print "## base=$base, outdir=$outdir\n";

my $tmp="$outdir/tmpaln";
$base="$outdir/novoalign";

#Set some outputs
my $align1 = "$base/$smname.bam";
my $bedgraph = "$base/$smname.bedgraph";
my $vcf = "$base/$smname.vcf";

$samplename=$smname;
$samplename=~s/\s+|\?|\$|\^//g;
$libname=$samplename;




my $addopt=" ";
addAlignerOpts();

system("mkdir -p $tmp");
system("mkdir -p $base");
novoindex($genomefasta,$index) unless -e $index;

#my $novoversion=`novoalign version`;

##############################################
#####  Alignment stage
##############################################
#
print "
 ##############################################
 ######  Alignment stage
 ###############################################
";
if (scalar keys %pairs ==1) {
    #only 1 set of files
     NovoSinglePairedReads();
}else {
    #else we combine multiple sample files
    NovoCombineReads();
}





exit;
############################################
#Run some utilities to summarize alignments
############################################
unless ($norun) {
    system("bash /home/bin/novopipeline/scripts/novoversion.sh > $base/summary.txt");
    system("bamtools stats -in $align1  |sed 's/:/,/' >> $base/summary.txt");
    print "$b.INFO Calculating reference coverage stats\n";
    system("mkdir -p $base/rplots; bash /home/bin/novopipeline/scripts/genomeCoverage.sh $align1 > $base/genomeCoverage.txt && $appbase/genomeCovplot.R  $base/genomeCoverage.txt $base/rplots/genomecov.png");
    system("bash /home/bin/novopipeline/scripts/covsummary.sh  $base/genomeCoverage.txt > $base/coveragestats.txt ") if -e "$base/genomeCoverage.txt" ;
    system("bedtools genomecov -ibam $align1 -bg > $bedgraph");
    system("gzip -f $bedgraph") if -e $bedgraph;

}
#clean up genome index;
unlink($index) if -e $index;
########################################
##### Variant Calling ##################
########################################
if ($callvars) {
    print "
     ##############################################
     ######  Variant Calling stage
     ###############################################
    ";
    print "# $b.INFO\tVariant Calling starting\n";
    my $cmd="perl $script call -reference $genomefasta $align1 $vcf";
    $cmd.=" -rmdup ";
    print "$cmd\n:";
    system($cmd) unless $norun;
    if (-e $vcf) {
        print "Counting variants\n";
        system("bash /home/bin/novopipeline/scripts/variantcount.sh $vcf >> $base/summary.txt");
        print STDERR "Compressing and indexing VCF";
        system("bgzip -f $vcf");
        system("tabix -p vcf $vcf.gz");
    }
    system(" mkdir -p $base/logs && mv $base/*stdout $base/logs/");
}


##### Final actions ##################################
print "$b.INFO:\tResults written to $base\n";
print "$b.INFO:\tCleaning up files\n";
#system("cp -fr $tmp/*.txt $tmp/*.rplots $base/") unless $norun;
system("mkdir -p $base/logs && mv $tmp/*inserts.txt   $tmp/*err.txt $tmp/*log.txt $tmp/*.stdout.txt $base/logs/") unless $norun;
system("mkdir -p $base/qcal && cp $tmp/*qcal.txt $base/qcal ");
system("mv $base/*log.txt $base/*err.txt $base/logs");
print "Expunging temp directory $tmp\n";
system("rm -fr $tmp") unless $norun;


cleanSummary() unless $norun;
system("mv $base $outdir/$smname ");
print "##### Novoalign App ended #####\n";
exit 0;


#### END of Pipeline





#### Functions/subroutines
#Align reads function
sub alignReads {
    my $cmd=shift;
    my $reads=shift;

    unless ($reads) {
        print STDERR "No reads found in this sample, cant run $cmd\n";
        print STDERR "App exiting\n";
        exit 500;
    }
     my $out=shift;
    my $runcmd = "bash $script $reads $out ";
    print "\t$b.INFO Running batch: $runcmd\n";  
    system($runcmd) unless $norun;  
}


sub addAlignerOpts {
    if ($readNum && $trimadapt ) {
        $addopt=" -aligner_opt -#$readNum  ";
    }elsif ($readNum && !$trimadapt) {
        $addopt="  -aligner_opt  -#$readNum  ";
    }elsif ($trimadapt && !$readNum ) {
        #$addopt=" -aligner_opt \" -a \" ";
    }
}

sub cleanSummary {
    if (-e "$base/Picard_duplicate_metrics.txt"){
        system("bash /home/bin/novopipeline/scripts/rmdup_summary.sh $base/Picard_duplicate_metrics.txt >>  $base/summary.txt");
    }
    open(IN,"$base/summary.txt") or die "$!";
    open(OUT,"+>$base/summary2.txt") or die "$!";
    while(<IN>){
        next if /^$/;
        next if /\*/;
        next if /Stats\s+for\s+BAM\s+/;
        print OUT $_;
    }
    close IN;
    close OUT;
    system("mv $base/summary2.txt $base/summary.txt");
}

#call the pipeline here
#for only 1pair of reads
sub NovoSinglePairedReads {
    print "$b.INFO Only 1 set of read files in this sample $sampleid\n";   
     my @ray=keys %pairs;
     my $inreads="";
        $sid=$ray[0];
        $read1= $pairs{$sid}{R1};
        $read2= $pairs{$sid}{R2};   
         if ($read1 && $read2) {
                   $inreads="$read1 $read2";
        }else {
                    $inreads=" $read1 ";
        }
	my $manifestpath=`find /data/ -name $manifestfile `;
	chomp $manifestpath;
        $cmd="bash $script $outdir/$smname  $inreads $genomefasta $manifestpath";
        print "AlignCMD: $cmd\n";
        system($cmd) unless $norun;
        #system("mv $tmp/rplots $outdir/$smname.rplots  ");
}



sub NovoCombineReads {
     print "$b.INFO\t > 1  read pair found, aligning separately\n";
    #Align FASTQs read1/read2 for each lane
    my $cmd="perl $script align  -samsample $smname -samlib $smname $addopt  -novoindex $index ";
    my $bamfiles="";
    my $cnt=0;
    foreach $sid (keys %pairs) {
            $cnt++;
            print "-------------- $sid Lane $cnt Aligning -----------------\n";
            my ($read1,$read2);
            $read1= $pairs{$sid}{R1};
            $read2= $pairs{$sid}{R2};
            if ($read1 && $read2) {
                     $paired=1; 
                     alignReads($cmd,"$read1 $read2","$tmp/$sid/out.bam");
            }else {
                    $paired=0;
                     alignReads($cmd,$read1,"$tmp/$sid/out.bam");
            }
            $bamfiles.=" $tmp/$sid/*.aln.bam ";
            if (-d "$tmp/$sid") {
                print "Copying files from $tmp/$sid  to $tmp/ \n";
                system("cp $tmp/$sid/*novoalign.log.txt $tmp/");
                system("cp $tmp/$sid/*out.txt $tmp/");
                system("cp $tmp/$sid/*err.txt $tmp/");
                system("cp $tmp/$sid/*qcal.txt $tmp/");
                system("cp $tmp/$sid/*inserts.txt $tmp/");
                system("mv $tmp/$sid/rplots $outdir/$sid.rplots  ");
                system("mv $tmp/$sid/*aln.bam $tmp/$sid.aln.bam; ");

               if ($cnt==1) {
                    system("cp -fr $outdir/$sid.rplots $base/rplots");
               }
               #cleanup folder
               system("rm -fr $tmp/$sid");
            }
    }
 }





sub up_to_date {
    my $fileA=shift;
    my $fileB=shift;
    my $timeA = ( -M $fileA );
    my $timeB = ( -M $fileB );
    #$current_time = time;

    #if fileA is older than fileB
    unless ($timeB) {
        return 0;
    } 
    if ($timeA > $timeB) {
        return 1;
     }else {
        return 0;
     }
}

    
sub novoindex {
    my $fasta=shift;
    my $index=shift;
    my $CMD="novoindex $index $fasta";
    print "INFO RUN: Indexing Reference $fasta\n$CMD\n";
    system($CMD) unless $norun;
}



sub parse_JSON {
        my $filename =shift;
        my $json_text = do {
           open(my $json_fh, "<:encoding(UTF-8)", $filename)
              or die("Can't open \$filename\": $!\n");
           local $/;
           <$json_fh>
        };

        my $json = JSON->new;
        my $data = $json->decode($json_text);


        my $readNum="";
        my $is_softclip=0;
        my $genomeId="Human";
        my $alignonly=0;
        my $trimadapt=0;
        my $sampleid="";
        my $projectid="";
        my $is_amplicon=0;
	my $manifestfile="";
	my $readdepth="";
       
        #  print Dumper($data) if $testjson; 
        
        my $prop = $data->{Properties}{Items};
        foreach $obj (@{$prop}) {
	
            #print Dumper($obj);
            #  print "######################################################################\n";
                foreach $key (keys %{$obj} ) {
	                if ($obj->{$key} eq "Input.readNum") {
                                $readNum= $obj->{Content};
                                if ($readNum =~/unlimited/ ) {
                                        $readNum="";
                                }
                        }
                        #Assign sample ID e.g. NA12878
                        if ($obj->{$key} eq "sample") {
                          $smname = $obj->{Content}{SampleId} if  $obj->{Content}{SampleId};
                              
                        }    
                        #Assign genome ID
                        if ($obj->{$key} eq "Input.genome-id") {
                                $genomeId= $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.softclip") {
                                $is_softclip= $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.alignonly") {
                                $alignonly = $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.callvars") {
                                $callvars = $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.trimadapt") {
                                $trimadapt = $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.sample-id") {
                                $sampleid = $obj->{Content}{Id};
                                print $key," is something \n";
                        }
                        if ($obj->{$key} eq "Input.project-id") {
                                $projectid = $obj->{Content}{Id};
                        }
                        #Is this an amplicon pipeline
                        if ($obj->{$key} eq "Input.is_amplicon") {
                                $is_amplicon = $obj->{Content};
                        }
                        if ($obj->{$key} eq "Input.manifestfile") {
				$manifestfile= $obj->{Items}[0]{Path};
                        }
                        if ($obj->{$key} eq "Input.readDepth") {
                                $readdepth = $obj->{Content};
                        }


                }

        }
	print "--------------------------------------------------------------";
        print "\nForm: sampleid is $sampleid\n";
        print "Form: projectid is $projectid\n";
        print "Form: readnum is $readNum\n";
        print "Form: softclip is $is_softclip\n";
        print "Form: genome  is $genomeId\n";
        print "Form: alignonly  is $alignonly\n";
        print "Form: trimadapt is $trimadapt\n";
        print "Form: call variants is $callvars\n";
        print "Form: sample name is $smname\n";
        print "Form: Amplicon button set to $is_amplicon\n";
	print "Form: Read depth $readdepth\n";
        print "Form: Manifest set to $manifestfile\n";
        $r{readNum}=$readNum;
        $r{genomeId}=$genomeId;
        $r{is_amplicon}=$is_amplicon;
        $r{is_softclip}=$is_softclip;
        $r{align_only}=$alignonly;
        $r{projectid}=$projectid;
        $r{sampleid}=$sampleid;
        $r{samplename}=$smname;
        $r{trimadapt}=$trimadapt;
        $r{manifestfile}=$manifestfile;
        $r{readdepth}=$readdepth;
          exit 0 if $is_debug;
        return \%r;
}


sub chooseGenome  {
    if ($genomeName =~/chr22/) {
        print "$b.INFO Downloading reference  hg19 chr22 genome sequence\n";
        #system("wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz -O chr22seq.fa.gz");
        #system("gunzip -f chr22seq.fa.gz");
        $genomefasta="chr22seq.fa";
    } elsif ($genomeName=~/coli/) {
            $genomefasta="/home/bin/novopipeline/fasta/ecoli.fasta";
     }elsif ($genomeName =~ /bacil/) {
            $genomefasta="/home/bin/novopipeline/fasta/bacillus.fasta";
     }elsif ($genomeName =~ /phix|phiX/) {
            $genomefasta="/home/bin/novopipeline/fasta/phiX.fasta";
     }elsif ($genomeName =~ /hiv1/) {
            $genomefasta="/home/bin/novopipeline/fasta/hiv1.fasta";
     }elsif ($genomeName =~ /hiv2/) {
            $genomefasta="/home/bin/novopipeline/fasta/hiv2.fasta";
     }else {
        print "$b.INFO: Genome is full hg19 ";
         $genomeName="hg19";
       # $genomefasta="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
    }

    if (! -e $genomefasta) {
        warn "Cannot locate genome reference file: $genomefasta\n";
        #exit 555;
    }
}

sub checkFormInputs {
    if (-e $appsession) {
        print STDERR "$b MSG: Reading in Form data from $appsession\n";
        $FormData = parse_JSON($appsession);
        if (!$FormData->{is_softclip}) {
            $alignopts.=" -o FullNW";
        }
        if ($FormData->{align_only}) {
            $align_only=1;
        }
        if ($FormData->{trimadapt}) {
            $trimadapt=1;
        }
        if ($FormData->{readNum}  && $FormData->{readNum} !~ /\;|\~|\%|\*|\)|\(|\@|\!/ ) {
            $alignopts.=" -#".$FormData->{readNum}." ";
         }
        if ($FormData->{genomeId}) {
                $genomeName= $FormData->{genomeId};
         }
         if ($FormData->{projectid}){
                $projectid=$FormData->{projectid}; 
           }
         if ($FormData->{sampleid}){
                $sampleid=$FormData->{sampleid}; 
           }
         if ($FormData->{samplename}){
                $smname=$FormData->{samplename}; 
           }
         if ($FormData->{is_amplicon}){
                $is_amplicon=$FormData->{is_amplicon}; 
           }
         if ($FormData->{manifestfile}){
                $manifestfile=$FormData->{manifestfile}; 
           }
        $outdir= "/data/output/appresults/$projectid/Novoalign.$sampleid"."_analysis";
        $outdir= "/data/output/appresults/$projectid";
        $indir = "/data/input/samples/$sampleid";
        print "$b.INFO: reads folder: $indir\n";
        print "$b.INFO: Output folder: $outdir\n";

     }else {
        print STDERR "$b: MSG No Appsession.json found\n";
     }
}


sub moreChecks {    
    if (scalar @ARGV==0 ){
        print $usage;
        exit 255;
    }

    unless ($outdir) {
        print STDERR  $usage;
        print STDERR "No outdir given, error";
        exit 255;
    }

    unless ( $genomefasta) {
        print STDERR $usage;
        print STDERR "No Genome reference given, error\n";
        exit 255;
    }
}
