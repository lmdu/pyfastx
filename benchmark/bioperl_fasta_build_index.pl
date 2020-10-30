use strict;
use Bio::DB::Fasta;

my $fasta_file = shift;

Bio::DB::Fasta->new($fasta_file);
