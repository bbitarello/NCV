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
my @lines;
my $line;

#must improve this to apply to all subdirectories;;

#$s1=200;
#$s2=28948;
#$s3=67004;

#$startp2=$s1+1;
#$endp2=$startp2+$s2;
#$startp3=$endp2+1;
#$endp3=$startp3+$s3;


for ($i=1; $i<=10; $i++){

	$file= "/mnt/sequencedb/PopGen/barbara/simulations/model2pop/s1/".$i."/1_2pops_withOUTGROUP.out";
#	$temp=$file."_".$i."txt";
#	$out=$file."new";

	open (FILE, "<$file") or die "cannot open $file!\n";
#	open (TEMP, ">$temp") or die "cannot open $temp!\n";
	undef $/; #read entire file at once.
	while (<FILE>){
        	(@fields)=split(/Genomes:\n/,$_);
	}
#	@lines=split(/\n/, @fields[1]);
	#print "@fields[1]\n//";
	$/="\n";
	while(<$fields[1]>){
	if ($line =~ /^p1/){
	push(@p1,$line);
	}
	elsif($line=~/^p2/){
	push(@p2,$line);
	}
	else{
	push(@p3,$line);
	}
	}
	
}
print "@p2\n";
close(FILE);
#sample 1, 60 and 60 from p1, p2 and p3 respectively.
exit;

