#!/usr/bin/perl -w

use strict;

# Number of lines to pick and file to pick from
# Error checking omitted!
my ($pick, $file) = @ARGV;

open(my $fh, '<', $file)
    or die "Can't read file '$file' [$!]\n";

# count lines in file
my ($lines, $buffer);
while (sysread $fh, $buffer, 4096) {
    $lines += ($buffer =~ tr/\n//);
}

# limit number of lines to pick to number of lines in file
$pick = $lines if $pick > $lines;

# build list of N lines to pick, use a hash to prevent picking the
# same line multiple times
my %picked;
for (1 .. $pick) {
    my $n = int(rand($lines)) + 1;
    redo if $picked{$n}++
}

# loop over file extracting selected lines
seek($fh, 0, 0);
while (<$fh>) {
    print if $picked{$.};
}
close $fh;
