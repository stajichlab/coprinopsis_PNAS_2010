#!/usr/bin/perl -w
use strict;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::PrimarySeq;
my $db = Bio::DB::Fasta->new(shift);
my $out = Bio::SeqIO->new(-format => 'fasta');
while(<>) {
 my @line = split;
 my $segment = $db->seq($line[0], $line[3] => $line[4]); 
 my $id = $line[-1];
 $id =~ s/ID=//;
 $out->write_seq(Bio::PrimarySeq->new(-id => $id, -desc => sprintf("%s:%d..%d",$line[0],$line[3],$line[4]),
				      -seq=> $segment));

}
