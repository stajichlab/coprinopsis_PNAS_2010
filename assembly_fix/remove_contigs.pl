#!/usr/bin/perl -w
use strict;

my $skip = shift;
open(my $fh => $skip) || die;
my %contigs;
while(<$fh>) {
    my ($ctg) = split;
    $contigs{$ctg}++;
    $ctg =~ s/ccin_Contig/U/;
    $contigs{$ctg}++;
}
while(<>) {
    my @row = split(/\t/,$_);
    if(exists $contigs{$row[0]} ) {
	next;
    } 
    print;
}
