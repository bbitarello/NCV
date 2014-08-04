#!/usr/bin/perl -w

#clean slim output files

use strict;

#declare variables:

my $input;
my $yes;
my $no;
my $p1;
my $p2;
my $p3;
my $p23;
my $s1;
my $s2;
my $s3;
my $startp2;
my $endp2;
my $startp3;
my $endp3;
my $folder;
my $out;
my @p1;
my @p2;
my @p3;
my $file;
my $temp;
my $i;
my @fields;
my $fields;
my $outp2;
my $outp3;

#must improve this to apply to all subdirectories;;


#$s1=200;
#$s2=28948;
#$s3=67004;

$file="1_2pops_withOUTGROUP.log";
#$temp="temp.txt";
$out="outp1.txt";
#$outp2="outp2.txt";
#$outp3="outp3.txt";

open (FILE, "<$file") or die "cannot open $file!\n";
#open (TEMP, ">$temp") or die "cannot open $temp!\n";
open (OUT1, ">$out") or die "cannot open $out!\n";
#open (OUT2, ">$outp2")or die "cannot open $outp2\n";
#open (OUT3, ">$outp3") or die "cannot open $outp3\n";
#undef $/; #read entire file at once.


while (<FILE>){

 #       (@fields)=split(/Genomes:\n/,$_);

	

#	print TEMP $fields[1];
#	print TEMP "\n";
#
	#now I must take the first 200 lines and -place in $p1, the following X and place in $p2 and the ramaining in $p3
#}

#close(FILE);
#close(TEMP);

#	open (YES, "<$temp") or die "oh, wait...bummer\n";

	#push lines to @p1, @p2 and @p3.
#	$/="\n";

#	while (<YES>){

	if($_=~ /^p1/){
	push(@p1, $_);
	}
	elsif($_ =~ /^p2/){
	push(@p2, $_);
	}
	elsif($_=~/^p3/){
	push(@p3,$_);
	}	
}

print OUT1 @p1;
print OUT1 @p2;
print OUT1 @p3;
#print OUT2 @p2;
#print OUT3 @p3;

#print OUT @p2;
#print OUT @p3;
close(OUT1);
close(FILE);

exit;

