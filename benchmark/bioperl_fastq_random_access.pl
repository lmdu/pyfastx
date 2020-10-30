use strict;
use Bio::Index::Fastq;

my $list_file = shift;
my $fastq_file = shift;
my $index_file = $fastq_file . ".index";

my $inx = Bio::Index::Fastq->new('-filename' => $index_file);

open(FH, $list_file);

while (<FH>) {
	chomp;
	my $seq = $inx->get_Seq_by_id($_);
	print $seq->seq;
}
