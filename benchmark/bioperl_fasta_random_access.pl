use Bio::DB::Fasta;
use strict;

my $list_file = shift;
my $fasta_file = shift;

my $db = Bio::DB::Fasta->new($fasta_file);

open(FH, $list_file);
while (<FH>) {
	chomp;
	my $seqstr = $db->seq($_);
	print $seqstr;
}
