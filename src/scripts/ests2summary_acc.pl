#!/usr/bin/perl -w
use strict;

# turn the EST library file (fasta) with these headers
# >gi|186828799|gb|FG068291.1|FG068291 CCOS1H9_T7 Coprinus cinereus osmotic shock mycelial cDNAs Coprinopsis
# >gi|61244614|gb|DN593839.1|DN593839 CCMIN9E12 Coprinus cinereus minimal medium exponentially 
# into a table of per-library counts and accession numbers (collapsed)
# the file used for the paper is in data/coprinus_ests.headers.fa.bz2

=head1 NAME

 est2summary_acc.pl - convert EST file to summary table by-library

=head1 SYNOPSIS 

 bzcat data/coprinus_ests.headers.fa.bz2 | perl scripts/est2summary_acc.pl > summary.table

=head1 DESCRIPTION

For the genome paper, summarizing the subimitted ESTs by the their library and accession numbers (condensed).

=head1 AUTHOR - Jason Stajich

Jason Stajich - jason_stajich[AT]berkeley.edu

=cut

my %groups;
while (<>) {
  if(s/^>gi\|(\d+)\|gb\|[\d\w\.]+\|(\S+)\s+(\S+)\s+//) {
    my ($gi,$acc,$clone) = ($1,$2,$3);
    my $comment = $_;
    my ($acc_pref,$acc_num);
    if ( $acc =~ /^([A-Z]+)(\d+)/ ) {
      ($acc_pref,$acc_num) = ($1,$2);
    } else {
      die;
    }
    if ($clone =~ /^(CC[A-Z]+)(\d+[A-H]\d+)(_T7)?/) {
      my ($group,$well) = ($1,$2);
      push @{$groups{$group}}, [$acc_pref,$acc_num,$gi,$well];
    } elsif ( $clone =~ /^CK/ ) {
      push @{$groups{'CK'}}, [$acc_pref,$acc_num,$gi,$clone];
    } elsif ( $clone =~ /^[a-z]\d[a-z]\d+cc/ ) {
      push @{$groups{'K+6'}}, [$acc_pref,$acc_num,$gi,$clone];
    } else {
      warn("pattern couldn't match $clone\n");
    }
  }
}

my $sum;
for my $group (keys %groups ) {
  my %pref;
  for my $p ( @{$groups{$group}} ) {
    push @{$pref{$p->[0]}}, $p->[1];
  }
  my @keep;
  while ( my ($p,$nums) = each %pref ) {
    push @keep, join(",",map { $_ = $p.$_;
			     s/-/-$p/g; 
			     $_;} 
		     &collapse_nums(sort { $a <=> $b }@{$nums}));
  }
  printf "%s\t%d\t%s\n", $group, scalar @{$groups{$group}},join(",", @keep);
  $sum += scalar @{$groups{$group}};
}

print "$sum total ESTs\n";


sub collapse_nums {
  # This is probably not the slickest connectivity algorithm, but will do for now.
  my @a = @_;
  my ($from, $to, $i, @ca, $consec);
    
  $consec = 0;
  for ($i=0; $i < @a; $i++) {
    not $from and do{ $from = $a[$i]; next; };
    # pass repeated positions (gap inserts)
    next if $a[$i] == $a[$i-1];
    if ($a[$i] == $a[$i-1]+1) {
      $to = $a[$i];
      $consec++;
    } else {
      if ($consec == 1) {
	$from .= ",$to";
      } else {
	$from .= $consec>1 ? "\-$to" : "";
      }
      push @ca, split(',', $from);
      $from =  $a[$i];
      $consec = 0;
      $to = undef;
    }
  }
  if (defined $to) {
    if ($consec == 1) {
      $from .= ",$to";
    } else {
      $from .= $consec>1 ? "\-$to" : "";
    }
  }
  push @ca, split(',', $from) if $from;

  @ca;
}
