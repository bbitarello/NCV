#!/bin/perl


#####################################################################################
#	Make positions file for splitting jobs for NCV
#	Barbara Bitarello
#	Last modified: 16.01.2014
#################################################################


$input=$ARGV[0];

$first=$ARGV[1];

$last=$ARGV[2]

$skip=$ARGV[3];



$out='out_positions.txt';


open(INPUT, $input) or die 'cannot open $input\n';
open(OUT,$out) or die 'cannot open $out\n';


print OUT $

	
