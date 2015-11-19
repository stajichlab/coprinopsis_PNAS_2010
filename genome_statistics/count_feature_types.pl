#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;
use constant { 
    start  => 0,
    end    => 1,
    strand => 2,
};

my $genome_size = 36_294_355;

my $introns_per_gene;
my %types;
my %genes;
my $gene;
while(<>) {
    next if /^\#/;
    my ($chrom,$src,$type,$start,$end,undef,$strand) = split(/\t/,$_);
    if( $type eq 'mRNA' ) {
	if( /Name=(\S+)/ ) {
	    $gene = $1;
	} else {
	    warn("cannot find gene name for $_");
	}
    }
    my $length = ($end - $start) +1;
#     if( $length == 1 ) { 
# 	warn $_;
#     }
    if( $type eq 'exon' ) {
	push @{$genes{$gene}}, [$start,$end,$strand eq '+' ? 1 : -1];
    }
    if( $type eq 'mRNA' ) {
	next;
    }
    push @{$types{$type}}, $length;
}


for my $gene ( keys %genes ) {
    my $intron_count = 0;
    my $transcript_length = 0;
    if( @{$genes{$gene}} == 1 ) {
	push @{$types{'single_exon'}}, 0;
	my ($exon) =  @{$genes{$gene}};
	$transcript_length = $exon->[end] - $exon->[start] +1;
    } else {
	my $lastexon;
	for my $exon ( sort { $a->[start] * $a->[strand] <=> 
				  $b->[start] * $b->[strand] }
		       @{$genes{$gene} } ) {
	    $transcript_length += ($exon->[end] - $exon->[start] +1);
	    if( $lastexon ) {
		$intron_count++;
		my $intron_len = 0;
		if( $exon->[strand] == 1 ) {
		    $intron_len = ( $exon->[start]-1 - $lastexon->[end]+1)+1;
#		    warn(sprintf("len=%d rev intron %d..%d for %d..%d - %d..%d\n",
#				 $intron_len,
#				 $lastexon->[end]+1,$exon->[start]-1,
#				 $lastexon->[start],$lastexon->[end],
#				 $exon->[start],$exon->[end],
#				 ));
		} else {
		    $intron_len = ( $lastexon->[start]-1 -
				    $exon->[end]+1) + 1;
#		    warn(sprintf("len=%d rev intron %d..%d for %d..%d - %d..%d\n",
#				 $intron_len,
#				 $exon->[end]+1, $lastexon->[start]-1,
#				 $exon->[start],$exon->[end],
#				 $lastexon->[start],$lastexon->[end],
#				 ));
#		    exit;
		}
		push @{$types{'intron'}}, $intron_len;
	    }
	    $lastexon = $exon;
	}
    }
    push @{$types{'intron_count'}}, $intron_count;
    push @{$types{'mRNA'}}, $transcript_length;	    
}

print join("\t", qw(TYPE MAX MIN N MEAN MEDIAN SUM SUM_MB)),"\n";
for my $type ( sort keys %types ) {
    my $stats = Statistics::Descriptive::Full->new();
    $stats->add_data(@{$types{$type}});
    print join("\t", $type,
	       $stats->max,
	       $stats->min,
	       $stats->count,
	       sprintf("%.2f",$stats->mean),
	       $stats->median,
	       $stats->sum,
	       sprintf("%.2f",$stats->sum / 1_000_000),
	       sprintf("%.2f",(100 * ($stats->sum /$genome_size))),
	       ),"\n";
}
