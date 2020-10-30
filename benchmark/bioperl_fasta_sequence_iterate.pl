use strict;
use Bio::SeqIO;

my $fasta_file = shift;
my $fasta_db = Bio::SeqIO->new(-file=>$fasta_file, -format=>'fasta');

while (my $seqobj = $fasta_db->next_seq) {
	$_ = $seqobj->seq;
}
