#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $skip = shift;
open(my $fh => $skip) || die;
my %contigs;
while(<$fh>) {
    my ($ctg) = split;
    $contigs{$ctg}++;
    $ctg =~ s/ccin_Contig/U/;
    $contigs{$ctg}++;
}
my $in = Bio::SeqIO->new(-format => 'fasta',-fh => \*ARGV);
my $out = Bio::SeqIO->new(-format => 'fasta');
while(my $seq = $in->next_seq ) {
    if( exists $contigs{$seq->id} ) {
	warn("skipping ",$seq->id,"\n");
	next;
    } 
    $out->write_seq($seq);
}
