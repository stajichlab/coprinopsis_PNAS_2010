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
GetOptions(
	   'v|verbose!'  => \$debug,
	   'u|user:s'    => \$user,
	   'p|pass:s'    => \$pass,
	   'host:s'      => \$host,
	   'db|dbname:s' => \$dbname,
	   'pair:s'  => \$pairfile,
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

# compute the adjacency
my %adj;

my @names;
my $iter = $dbh->get_seq_stream(-type => 'scaffold:chromosome');
while( my $chrom = $iter->next_seq ) {
    push @names, $chrom->seq_id;
}


my %genome_dist;
open(my $ofh => ">genome_pair_config.dat") || die $!;
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
	    print $ofh join("\t", $last_gene, $gene_name, $cfg),"\n";
	} else {
	    $adj{$gene_name}->{left} = undef;
	}	

#	$genome_dist{$strand2cfg{$adj{$gene_name}->{strand}.",".
#				     $adj{$last_name}->{strand}}}++;
	$last_gene = $gene_name;
    }
    if( $last_gene ) {
	$adj{$last_gene}->{right} = undef;
    }
}
close($ofh);
open(my $fh => $pairfile) || die $!;
my @header = split(/\s+/,<$fh>);
my %header;
{ 
    my $i = 0;
    %header = map { $_ => $i++ } @header;
}

my %seen; # unique counter
splice(@header,2,0,qw(geneconfig chrom 
		      model1_start model1_stop model1_strand 
		      model2_start model2_end model2_strand));
print join("\t", @header),"\n";
while(<$fh>) {
    my $line = $_;
    my $ii =0;
    my ($one,$two,$ds,$dn,$dsdn) = split;
    for ( $one, $two ) {
	s/T\d+$//;
	s/^Ccin://;
    }
    if( (defined $adj{$one}->{right} && $adj{$one}->{right} eq $two) ||
	(defined $adj{$one}->{left} && $adj{$one}->{left} eq $two )) {
	if( $seen{"$one,$two"}++) {
	 #   warn("seen $one,$two\n");
	    next;
	}
	my ($left,$right) = sort { $adj{$a}->{start} <=> $adj{$b}->{start} } ($one,$two);
	my $configuration = $strand2cfg{ $adj{$left}->{strand}.",".$adj{$right}->{strand} };
	
	print join("\t",$left,$right,$configuration,
		   (map { $adj{$left}->{$_} } qw(seqid start end strand)), 
		   (map { $adj{$right}->{$_} } qw(start end strand)),
		   $ds, $dn, $dsdn), "\n";
    } 

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
