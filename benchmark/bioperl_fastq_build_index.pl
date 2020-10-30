use strict;
use Bio::Index::Fastq;

my $fastq_file = shift;
my $index_file = $fastq_file . ".index";

my $inx = Bio::Index::Fastq->new('-filename' => $index_file, '-write_flag' => 1);
$inx->make_index($fastq_file);
