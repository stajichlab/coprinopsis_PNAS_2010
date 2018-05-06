#!/usr/bin/perl -w
use strict;
my $l = <>;
$l=~s/dS/dS_inv/;
print $l;

while(<>) {
    my @row = split;
    next if $row[2] == 0;
    $row[2] = 1/$row[2];
    print join("\t", @row), "\n";
}
