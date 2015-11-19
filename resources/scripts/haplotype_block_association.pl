#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::DB::SeqFeature;
my $debug = 0;

GetOptions('v|verbose|debug!' => \$debug,
	   );

my $dir = shift;
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
					 -autoindex=>1,
					 -dir     => $dir);
chomp(my $header = <>);
my @rows = split(/\t/,$header);
my $i =0;
my %hdr = map { $_ => $i++ } @rows;
print join("\t", qw(CHROM START END NAME LABEL STATUS LENGTH
		    GENE_DENSITY REPEAT_DENSITY)),"\n";
while(<>) {
    chomp;
    my ($chrom,$start,$end,$name,$label,$hot) = split(/\t/,$_);
    my $segment = $db->segment($chrom,$start => $end);
    my @genes = $segment->features(-type => 'gene');
    my $gene_count = scalar @genes;
    my $gene_density = sprintf( "%.2f", $gene_count / 
				($segment->length / 1000));
    my @repeats = $segment->features(-type => "match_part:repeatmasker");
    my $repeat_density = sprintf( "%.2f", scalar @repeats / 
				  ($segment->length / 1000));
    print join("\t", $chrom,$start,$end,$name,$label,$hot,
	       int($segment->length / 1000),
	       $gene_density,$repeat_density),"\n";
}
