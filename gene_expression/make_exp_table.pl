#!/usr/bin/perl -w
use strict;

my $oligomap = 'oligo_map/oligo_map_33009_01-28-10.csv';


open(my $fh => $oligomap) || die $!;
my (%oligo2gene,%lookup);

my $oligohdr = <$fh>;
chomp($oligohdr);
my $j = 0;
my %oligohdr = map { $_ => $j++ } split(',',$oligohdr);

while(<$fh>) {
    chomp($_);
    my @row = split(',',$_);
    next unless( defined $row[ $oligohdr{'IU name'} ] );
    $oligo2gene{$row[0]} = $row[ $oligohdr{'IU name'} ];
}
close($fh);

my @id_cols = map { uc($_)} ('UNC NAME', 'Broad v2');#'Broad v2','Broad v2 Location');

my %table;
my @experiments;
for my $file ( @ARGV) {
    open($fh => $file) || die $!;
    my $stem;
    if ($file =~ /CC5_A2B2g_(\w+)_/ ) {
	$stem = $1;
    } else {
	warn("Cannot find stem for $file\n");
    }
    push @experiments, $stem;
    my $header = <$fh>;
    chomp($header);
    my @header = split(',',$header);
    my $value_col;
    my $i = 0;
    my %header2col;
    for my $c ( @header ) {
	
	if( $c =~ /M\(/) {
	    $value_col = $i;
	}
	$header2col{uc $c} = $i++;
    }

    for my $idcol ( @id_cols) {
	if( ! exists $header2col{$idcol} ) {
	    die("cannot find '$idcol' in header (",join(",",@header),")\n");
	}
    }
    while(<$fh>) {
	chomp;
	my @row = split(',',$_);
	$lookup{$row[0]} = [ map { $row[$header2col{$_} || 'NA' ] } @id_cols ];
	$table{$row[0]}->{$stem} = $row[ $value_col ];
    }
}

#print join("\t", qw(CCIN_ID OLIGO_NAME BROAD_ID BROAD_LOCATION),
#	   @experiments),"\n" ;
print join("\t", qw(CCIN_ID OLIGO_NAME BROAD_GENE IU_GENE),@experiments),"\n" ;

for my $id ( sort { $a <=> $b } keys %table ) {
    print join("\t", $id, @{$lookup{$id}}, 
	       $oligo2gene{$lookup{$id}->[0]},
	       map { defined $table{$id}->{$_} ? 
			 $table{$id}->{$_} : 'NA' } @experiments),"\n";
}
