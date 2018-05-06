#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $exp = 'ccin_clean2.exp';
my $dup = 'dup_pair.tab';
my $all = 'genome_pair_config.dat';

GetOptions(
	   'e|exp:s' => \$exp,
	   'd|dup:s' => \$dup,
	   'a|all:s' => \$all,
	   );

my %all_adj;
my %pairs;

open(my $fh => $dup) || die $!;
my $header = <$fh>;
my %duplicates;
while(<$fh>) {
    my @row = split;
    my ($pair1,$pair2,$config,$chrom,@rest) = split;
    my ($ds,$dn,$dsdn) = splice(@rest,-3,3);
    my $dnds = $dsdn ne 'NA' ? 1/$dsdn : '-1';
    $ds = -1 if $ds eq 'NA';
    $dn = -1 if $dn eq 'NA';

    $duplicates{$pair1}->{$pair2} = [$ds,$dn,$dnds];
    $duplicates{$pair1}->{$pair1} = [$ds,$dn,$dnds];
}

open($fh => $exp) || die $!;
my ($id,$oligo,$broad,$iu,@experiments) = split(/\s+/, <$fh>);
my %fhs;
for my $exp ( @experiments) { 
    open(my $ofh => ">$exp\_dups.exp")|| die $!;
    print $ofh join("\t", qw(LEFT RIGHT CONFIG DUPLICATED DN DS DNDS 
			     EXP_LEFT EXP_RIGHT)),"\n";
    $fhs{$exp} = $ofh;
}

my @results;
my %exp_results;
while(<$fh>) {
    ($id,$oligo,$broad,$iu,@results) = split;
    $exp_results{$broad} = [@results];
}


open($fh => $all) || die $!;
while(<$fh>) {
    my ($left,$right,$config) = split;
    my ($dn,$ds,$dnds) = (-1,-1,-1);
    my $duplicated = 'NO';
    if( $duplicates{$left}->{$right} ) {
	($dn,$ds,$dnds) = @{$duplicates{$left}->{$right}};
	$duplicated = 'YES';
    }
    my $i = 0;
    for my $e ( @experiments ) {
	my $ofh = $fhs{$e};
	if( defined $exp_results{$left}->[$i] &&
	    $exp_results{$left}->[$i] ne 'NA' &&
	    defined $exp_results{$right}->[$i] &&
	    $exp_results{$right}->[$i] ne 'NA') {
	    
	    print $ofh join("\t", $left,$right,$config, 
			    $duplicated,$dn,$ds,$dnds,
			    $exp_results{$left}->[$i],
			    $exp_results{$right}->[$i]),"\n";
	}
	$i++;
    }
}

