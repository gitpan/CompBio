use strict;
use Test::More tests => 6;
use warnings;

# use a BEGIN block so we print our plan before MyModule is loaded
BEGIN {
	use_ok('CompBio::DB');
    $| = 1;
} # BEGIN
my $DEBUG = 0;

my $cbdb1 = CompBio::DB->new;
ok(ref $cbdb1 eq "CompBio::DB","Default datbase handle creation");

my $cbdb2 = CompBio::DB->new({user => "nobody"});
ok(ref $cbdb2 eq "CompBio::DB","Database handle using hashref");

my $cbdb3 = CompBio::DB->new((user => "nobody"));
ok(ref $cbdb3 eq "CompBio::DB","Database handle using hash");

my $AR_seqs = $cbdb1->get_aa_seq("SSUB_ECOLI");
ok($$AR_seqs[0][0] eq "SSUB_ECOLI"
	&& $$AR_seqs[0][1] =~ /TPLLLNAVSKHYAENIVLNQLDLHIPAGQFVAVVGRSGGGKSTLLRL/,
	"Single aa sequence selected");

my $AR_ids = [qw(UUP1_HAEIN PHNL_ECOLI LOLD_NEIMA BIO5_YEAST TYRP_ECOLI)];
$AR_seqs = $cbdb2->get_aa_seq($AR_ids);
ok(@$AR_ids == 5,"Five sequences selected");

my $AR_defs = $cbdb3->get_annotation("SSUB_ECOLI");
printout($AR_defs);
print "====\n";
my $AR_cds = $cbdb2->get_cds("b2891",{database => "genome"});
printout($AR_cds);
print "====\n";
$AR_cds = $cbdb1->get_cds([qw(26738 25277 26696)],{database => "genome",
	ISKEY => 1});
printout($AR_cds);
print "====\n";


sub printout {
	my $AR = shift;
	foreach my $AR_row (@$AR) {
		next unless @$AR_row;
		print join("\t",@$AR_row) , "\n";
	} # return list
} # printout
