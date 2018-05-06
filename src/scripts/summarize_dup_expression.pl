#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(max min sum);
my %shorter = ( 'HeadtoHead' => 'HH',
		'HeadtoTail' => 'HT',
		'TailtoTail' => 'TT',
		'NotAdjacent'=> 'NAD');
my @header = split(/\s+/,<>);
my $i = 0;
my %header = map { $_ => $i++ } @header;

my (%dcount,%ar);
while(<>) {
    my @row = split;
    my $type = $row[ $header{geneconfig} ];
    my $id = join(",",sort { $a <=> $b}
		  ($row[$header{'model1_cat_exp'}],
		   $row[$header{'model2_cat_exp'}]));
    $ar{$id}++;
    $dcount{$type}->{$id}++;    
    $dcount{$type}->{ALL}++;
}
my @combos = (sort keys %ar, 'ALL');
print join("\t", qw(TYPE), @combos),"\n";
my %r;
for my $type ( keys %dcount ) {
    print join("\t", $type, 
	       map { $r{$_} += $dcount{$type}->{$_} || 0;
		     $dcount{$type}->{$_} || 0} @combos),"\n";
}
print join("\t", 'ALL', (map { $r{$_} } @combos)),"\n";

my @range = (-1,0,1);
for my $type ( keys %dcount ) {
    print "\n";
    my $t = $shorter{$type};
    print join("\t",$t, @range),"\n";
    for my $r ( @range ) {
	print join("\t", $r, 		   
		   map { $dcount{$type}->{"$r,$_"} || ''}
		   @range),"\n";
    }
}
