
use Data::Dumper;
my %h=();

my $manifest=shift;
open(IN,$manifest) or die "$!";

my @hd;
my $nn;
while(<IN>) {
    s/\cM//g;
        chomp;
        my @F=split(/\t/,$_);
     if (/^\[Prob/ .. /^\[Targets/)  {
        $nn=$. if /^\[Targets/;
        next if /^\[/;
        if (/^Target/) {
                 @hd=@F;
        }else {
              head2dat(\@hd,\@F);
        }
      } elsif ($nn && $nn+1 .. eof)  {
        if (/^TargetA/) {
                print "Now T\n";
                 @hd=@F;
        }else {
             head2dat(\@hd,\@F);
        }
      } elsif ($nn && $nn+1 .. eof)  {
            
     
     }

}


sub head2dat {
    my $hd=shift;
    my $d=shift;
    my @data=@$d;
    my @head=@$hd;
    my %h=();
    my $i=0;
    foreach $v (@data) {
          my $k= $head[$i];
          $h{$k}=$v;
          $i++;
     }
     print Dumper(\%h);
}
