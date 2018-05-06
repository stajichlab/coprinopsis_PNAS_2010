#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::SeqFeature::Store;
use Env qw(HOME);

my %strand2cfg = ( '-,-' => 'HeadtoTail',
		   '-,+' => 'TailtoTail',
		   '+,+' => 'HeadtoTail',
		   '+,-' => 'HeadtoHead',
		   );

my $type = 'gene:CC2_FINAL_CALLGENES_1';
my $target_sp = 'Ccin';
my $pairfile = 'Ccin_dsdn.tab';

my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $output;
my $exp = 'ccin_clean2.exp';
my $pairsfile = 'Ccin_dsdn.tab';

GetOptions(
	   'v|verbose!'  => \$debug,
	   'u|user:s'    => \$user,
	   'p|pass:s'    => \$pass,
	   'host:s'      => \$host,
	   'db|dbname:s' => \$dbname,
	   'pair:s'      => \$pairfile,
	   'e|exp:s'     => \$exp,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

my @names;
my $iter = $dbh->get_seq_stream(-type => 'scaffold:chromosome');
while( my $chrom = $iter->next_seq ) {
    push @names, $chrom->seq_id;
}


my %adj;
for my $seqid ( @names ) {
    my $segment = $dbh->segment($seqid);
    my @genes = map { $_->[0] } 
    sort { $a->[1] <=> $b->[1] } 
    map { [$_, $_->start] } $segment->features(-type => $type);
    
    my $last_gene;
    for my $gene ( @genes ) {
	my $gene_name = $gene->name;

	$adj{$gene_name}->{gene}   = $gene;
	$adj{$gene_name}->{seqid}  = $seqid;	
	$adj{$gene_name}->{start}  = $gene->start;
	$adj{$gene_name}->{end}    = $gene->end;
	$adj{$gene_name}->{strand} = $gene->strand < 0 ? '-' : '+';

	if( $last_gene ) {
	    $adj{$gene_name}->{left} = $last_gene;	    
	    $adj{$last_gene}->{right} = $gene_name;
	    my $cfg = $strand2cfg{$adj{$gene_name}->{strand}.",".
				      $adj{$last_gene}->{strand}};
	    #print $ofh join("\t", $last_gene, $gene_name, $cfg),"\n";
	} else {
	    $adj{$gene_name}->{left} = undef;
	}	
	$last_gene = $gene_name;
    }
    if( $last_gene ) {
	$adj{$last_gene}->{right} = undef;
    }
}

open(my $fh => $exp) || die $!;
my ($id,$oligo,$broad,$iu,@experiments) = split(/\s+/, <$fh>);
my %fhs;

my @results;
my %exp_results;
while(<$fh>) {
    ($id,$oligo,$broad,$iu,@results) = split;
    $exp_results{$broad} = [@results];
}

open($fh => $pairfile) || die $!;
my @header = split(/\s+/,<$fh>);
my %header;
{ 
    my $i = 0;
    %header = map { $_ => $i++ } @header;
}

my %seen; # unique counter
splice(@header,2,0,qw(geneconfig  
		      model1_expression model1_cat_exp
		      model2_expression model2_cat_exp
		      model1_chrom model1_start model1_stop model1_strand 
		      model2_chrom model2_start model2_end model2_strand));
$header[-1] = 'dnds';
for my $e ( @experiments) { 
    open(my $ofh => ">$e\_onlydups.exp")|| die $!;
    print $ofh join("\t", @header),"\n";
    $fhs{$e} = $ofh;
}

while(<$fh>) {
    my $line = $_;
    my $ii =0;
    my ($one,$two,$ds,$dn,$dsdn) = split;
    for ( $one, $two ) {
	s/T\d+$//;
	s/^Ccin://;
    }
    my $dnds;
    if( $dsdn =~ /^\d+(\.\d+)?/ && $dsdn != 0) {
	$dnds = 1/ $dsdn;
    } else {
	$dnds = -1;
    }
    my $configuration = 'NotAdjacent';
    if( (defined $adj{$one}->{right} && $adj{$one}->{right} eq $two) ||
	(defined $adj{$one}->{left} && $adj{$one}->{left} eq $two )) {
	if( $seen{"$one,$two"}++) {
	    #   warn("seen $one,$two\n");
	    next;
	}
	my ($left,$right) = sort { $adj{$a}->{start} <=> 
				   $adj{$b}->{start} } ($one,$two);
	$configuration = $strand2cfg{ $adj{$left}->{strand}.",".
					  $adj{$right}->{strand} };
    }
    my $i = 0;
    for my $e ( @experiments ) {
	my $ofh = $fhs{$e};
#	if( defined $exp_results{$one}->[$i] &&
#	    $exp_results{$one}->[$i] ne 'NA' &&
#	    defined $exp_results{$two}->[$i] &&
#	    $exp_results{$two}->[$i] ne 'NA') {
	next if( ! exists $adj{$one} || ! exists $adj{$one}->{seqid} ||
		 ! exists $adj{$two} || ! exists $adj{$two}->{seqid} );

	my $one_exp = $exp_results{$one}->[$i] || 0;
	my $two_exp = $exp_results{$two}->[$i] || 0;
	$one_exp = 0 if $one_exp eq 'NA';
	$two_exp = 0 if $two_exp eq 'NA';
	print $ofh join
	    ("\t",$one,$two,$configuration,
	     $one_exp,
#	     ( $one_exp < -0.5 ? -1 : $one_exp > 0.5 ? 1 : 0), 
	     ($one_exp < 0 ? -1 : $one_exp > 0 ? 1 : 0),
	     $two_exp,
	     ($two_exp < 0 ? -1 : $two_exp > 0 ? 1 : 0),
#	     ( $two_exp < -0.5 ? -1 : $two_exp > 0.5 ? 1 : 0), 
	     (map { $adj{$one}->{$_} } qw(seqid start end strand)), 
	     (map { $adj{$two}->{$_} } qw(seqid start end strand)),
	     $ds, $dn, $dnds), "\n";    
	$i++;
    }
    #   }
}

sub read_cnf {
    my ($user,$pass) = @_;
    if( -f "$HOME/.my.cnf") {
        open(IN,"$HOME/.my.cnf");
        while(<IN>) {
            if(/user(name)?\s*=\s*(\S+)/ ) {
                $user = $2;
            } elsif( /pass(word)\s*=\s*(\S+)/ ) {
                $pass = $2;
            }
        }
        close(IN);
    }
    return ($user,$pass);
}
