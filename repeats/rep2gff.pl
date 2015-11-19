use strict;
use warnings;

my $ref;
while(<>) {
 if( /reference = (\S+)/ ) {
   $ref = $1;
 } else {
	my @line = split;
	my ($start,$end) = split(/\-/,$line[-1]);
	print join("\t", $ref, 'RepeatMasker','Repeat',$start,$end,'.','+','.',sprintf("ID=%s",$line[1])),"\n";
 }
}
