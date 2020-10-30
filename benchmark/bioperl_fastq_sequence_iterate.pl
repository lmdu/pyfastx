use strict;
use Bio::SeqIO;

my $fastq_file = shift;
my $fastq_db = Bio::SeqIO->new(-file=>$fastq_file, -format=>'fastq');

while (my $data = $fastq_db->next_dataset) {
	print $data->{-seq};
}
