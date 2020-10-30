use Bio::DB::Fasta;
use strict;

my $list_file = shift;
my $fasta_file = shift;

my $db = Bio::DB::Fasta->new($fasta_file);

open(FH, $list_file);
while (<FH>) {
	chomp;
	my ($seqid, $start, $end) = split('\t', $_);
	my $seqstr = $db->seq($seqid, int($start) => int($end));
	print $seqstr;
}
