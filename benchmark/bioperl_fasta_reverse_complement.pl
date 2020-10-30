use strict;
use Bio::DB::Fasta;

my $longname = shift;
my $fafile = shift;

my $db = Bio::DB::Fasta->new($fafile);
my $seqobj = $db->get_Seq_by_id($longname)->revcom;
print $seqobj->seq;
