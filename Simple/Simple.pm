package CompBio::Simple;

require 5.005_62;
use strict;
use warnings;
use Carp qw(cluck confess croak carp);
use CompBio;
use CompBio::Simple::ArrayFile;

use vars qw($AUTOLOAD);

our $VERSION = '0.44';
our $DEBUG = 0;

our %Auto_method = (
    tbl_to_fa => 1,
    tbl_to_ig => 1,
    fa_to_tbl => 1,
    ig_to_tbl => 1    
);

=head1 NAME

CompBio::Simple.pm  - Simple OO interface to some basic methods useful in
bioinformatics.

=head1 SYNOPSIS

use CompBio::Simple;

C<my $cbs = CompBio::Simple->new;>

C<my $fa_seq = $cbs->tbl_to_fa($tblseq,%params);>

$tblseq could be either scalar, array by reference, 2d array by reference, or a
hash by reference; return data type is same as submitted (see RETURN_TYPE option
under general method description). Scalar may be either sequence
records (one record per line), or a filename. Array's should contain an entire
sequence record (loci\tsequence) in each indexed position. 2D arrays should be
ids in the first column, sequences in the second, such that $array[2][1] would
be the aa sequence for the third record. Hash's should have the loci as the key,
and only the sequence as the value. Newlines at the ends of sequences in array and
hash types are unneccisarry and will be stripped off 

=head1 DESCRIPTION

Originally developed at the BioMolecular Engineering Research Center
(http://bmerc-www.bu.edu), this module is intended to take a number
of small commonly used methods, and make them into a single package.

The early versions of this module assumed installation on our local system.
Although I have tried to correct this in the current version, you may find
this package requires a litle twidling to get working. I'll try to leave
comments where I think it is most likely, but hopefully use of a relational
database and local setting changes in the base module CompBio.pm will have
taken care of it. If not _please_ email me at seanq@darwin.bu.edu with the
details.

Thanks!

=head1 Methods

A couple of important general notes. All methods will describe the key required
arguments and will also indicate which arg is the optional %params hash. %params
can be passed as a hash or reference. Also, most arguments that handle sequences
or ids accept input as scalar (as is or by reference), array (by reference),
2D array (by reference), or hash (by reference). scalar may be a filename. 
If I miss documenting it, try it just in case. Unless a parameter
option is passed to signify otherwise, I assume 'doing the right thing' is to
return the same type of data construct as recieved. The prefered method for
submitting sequence data is as a reference to an array, each index containing
one sequence record (and that preferably in table format).

All methods in Simple.pm accept a %params hash as the final argument. Options
available for a specific method are described in that methods section. Options
that can be defined for any (well, more than one at least) method are:

=over 4

=item RETURN_TYPE <type>

Value from: (SCALAR,REFSCALAR,ARRAY,2D,HASH)

This option overrides the default behavior of returning data in the same format
submitted and return data in the manner specified.

=item OUTFILE <filename>

Opens up <filename> for writing. Results of operation are written to file and 0
is returned.

=item DEBUG <value>

A numeric value supplied with this option sets the debug level for the specific
method call only. Does not affect the 'global' debug level set when creating a
new CompBio::Simple object.

=back

Also note that the params hash is shared by all methods called by the method
invoked, so params for those internal methods could also be added. The most
notable case of where this might be used is that check_type is actually called
by almost every other method (Simple tries to never assume or restrict what
sequence format is going to be used). So check_type's CONFIDENCE option could
be included to affect it's behavior.

=head2 new

Construct an object for invoking methods in CompBio::Simple.

Options:

DEBUG: Sets default debug level for all method calls using this object. Can
be overridden for a specific operation by defining a DEBUG option with the
method call. Higher integer values indicate more output. 1 and 2 are good
enough generally, 3 and higher can produce overwhelming amounts of output, 4
or higher will also cause the program to die where it would otherwise warn.

=cut
sub new {
    my ($proto,%parameters) = @_;
    my $class = ref($proto) || $proto;
    my $self = {};

    #handle params as nec. such as setting debug or changing env. variables
    $DEBUG = $parameters{'DEBUG'} || 0;
    (require warnings && import warnings) if $DEBUG; #??
    (require diagnostics && import diagnostics) if $DEBUG >= 2; #??
    $self->{'cbc'} = CompBio->new(%parameters);
    $self->{'_created'} = 1;

    return bless ($self,$class);
} # new

# maybe should explicitly perldoc module instead?
=head2 _help

Quits current application and uses perldoc to display the POD within it.

=cut
sub _help() {
    exec 'perldoc $0';
} # perldoc cheat

sub check_type {
    my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
        
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($AR_seqs),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    $params{CONFIDENCE} ||= 3;
    
    unless ($params{CONFIDENCE} > 0) {
        _error("check_type can't check type 0 or fewer times ($params{CONFIDENCE})");
        $params{CONFIDENCE} = 3;
    } # bad call
    
    if (@$AR_seqs < 3) {
       return ($$self{cbc}->check_type($AR_seqs,%params));
    } # no sense trying a bunch of times for just a sequence or two
    elsif (@$AR_seqs <= $params{CONFIDENCE}) {
        $params{CONFIDENCE} = int(@$AR_seqs / 1.5);
    } # No sense checking more times than there are items
    
    my %types = ();
    foreach (1 .. $params{CONFIDENCE}) {
        my $AR_test_subset = [(@$AR_seqs[int(rand(@$AR_seqs)),int(rand(@$AR_seqs))])];
        $types{$$self{cbc}->check_type($AR_test_subset,%params)}++;
    } # check type on set given # of times
    
    my @keys = keys %types;
    if (@keys > 1) {
        my $keys = join("\t",@keys);
        _error("Multiple types found in set using $params{CONFIDENCE} checks:\n$keys",1);
    } # different types found in set, die
    
    return _munge_seq_return($keys[0],$return_type,%params);
} # check_type

# man this error checking shit is a pain! :P
# get the rest of the methods docs from CompBio.pm
sub six_frame {
    my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
        
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($AR_seqs),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    my $AR_ret = _munge_array_to_scalar($self,$AR_seqs,%params);
    return _munge_seq_return($AR_ret,$return_type,%params);
} # six_frame

sub complement {
} # complement

# why aren't I checking format type and converting as necessary?!?
# simple should DTRT here! Think about data storage situations where
# converting to array will ruin type formating - and data. And where
# check_type will fail on sequence format incorrect. This should end
# being part of sanity shecking as well. Make sure I structure to fall
# as fast as possible.

sub dna_to_aa {
    my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
    
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($AR_seqs),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    my $AR_ret = _munge_array_to_scalar($self,$AR_seqs,%params);
    return _munge_seq_return($AR_ret,$return_type,%params);
} # dna_to_aa

# this is still way to brute force
sub _munge_array_to_scalar {
    my $self = shift;
    my $AR_seqs = shift;
    my %params = @_ ? @_ : ();
    my $caller = (caller(1))[3]; # fully qualified subroutine
    print "\t++  Got $caller from caller\n" if $DEBUG >= 3;
    $caller =~ s/.+:://;
    
    my $guess = $$self{cbc}->check_type($AR_seqs,%params);
    my $AR_ret = [];
    print "check_type returned $guess\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { print "seqcount in AR_seqs = ",scalar(@$AR_seqs),"\n" }
    
    if ($guess eq "RAW") {
        foreach my $dna (@$AR_seqs) {
            print "$dna\n" if ($DEBUG >= 2);
            my $result = $$self{cbc}->$caller(\$dna,%params);
            if (ref $result) { push(@$AR_ret,$$result) }
            else { push(@$AR_ret,$result) }
        } # easier to loop over one than add more tests or disallow mult. raw seqs
    } # no id's or fields to parse
    elsif ($guess eq "CDNA") {
        foreach my $seq (@$AR_seqs) {
            print "seq from AR_seqs:\n$seq<-\n" if ($DEBUG >= 2);
            my ($id,$dna) = split(/\s+/,$seq);
            print "Translating $dna\n" if ($DEBUG >= 2);
            my $result = $$self{cbc}->$caller(\$dna,%params);
            my $ret = "$id\t";
            if (ref $result) { $ret .= $$result }
            else { $ret .= $result }
            print "recieved $ret from $caller\n" if ($DEBUG >= 2);
            push(@$AR_ret,$ret);
        } # foreach
    } # translate cdna's
    # can't I compress these next two?
    elsif ($guess eq "FA") {
        $AR_seqs = $$self{cbc}->fa_to_tbl($AR_seqs,%params);
        print join("\n",@$AR_seqs),"\nReturned from fa_to_tbl\n" if ($DEBUG >= 2);
        foreach my $seq (@$AR_seqs) {
            (my $id,my $dna) = split(" ",$seq);
            print "Passing sequence from $id to $caller\n"  if ($DEBUG >= 1);
            my $result = ${$$self{cbc}->$caller(\$dna,%params)};
            my $ret = "$id\t";
            if (ref $result) { $ret .= $$result }
            else { $ret .= $result }
            print "Recieved $ret from $caller.\n" if ($DEBUG >= 2);
            push(@$AR_ret,$ret);
        } # foreach
        $AR_ret = $$self{cbc}->tbl_to_fa($AR_ret,%params);
    } # convert fa and then to aa
    elsif ($guess eq "IG") {
        $AR_seqs = $$self{cbc}->ig_to_tbl($AR_seqs,%params);
        foreach my $seq (@$AR_seqs) {
            my ($id,$dna) = split($seq);
            my $ret = $id . ${$$self{cbc}->$caller(\$dna,%params)};
            push(@$AR_ret,$ret);
        } # foreach
        $AR_ret = $$self{cbc}->tbl_to_ig($AR_ret,%params);
    } # convert fa and then to aa
    else { _error("Can't determine how to try to handle $guess in this context\n") }

    return $AR_ret;
} # really need a better name for this!

# allow for optional param FORMAT to declare format type for munging
# return_type should be declarable in params?
sub _munge_seq_input {
    my $self = shift;
    my $seq = (ref($self) && ref($self) !~ /ARRAY|HASH|SCALAR/) ? shift : $self;
    _help() if (! ref($seq) && $seq eq "help");
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($seq),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    my $return_type = "";
    my $AR_seqs = [];
    
    my $ref = ref($seq);
    if ($ref eq "ARRAY") {
        if (ref($$seq[0]) eq "ARRAY") {
            _error("Can't determine how to handle $seq\n") unless @{$$seq[0]} == 2;
            foreach (@$seq) { push(@$AR_seqs,join("\t",@$_)) }
            $return_type = "2D";
        } # collapse 2d array
        else {
            $AR_seqs = $seq;
            $return_type = "ARRAY";
        } # else
    } # array ref passed in
    elsif ($ref eq "HASH") {
        foreach my $k (keys %$seq) { push(@$AR_seqs,join("\t",($k,$$seq{$k}))) }
        $return_type = "HASH";
    } # if hash submitted
    elsif (! $ref || $ref eq "SCALAR") {
        my $new_ref = "";
        my @tmparray = ();
        # this was needed to avoid turning off warn for whole section under filetest
        # as a value with a newline blew up the file stat, but chomp ruins refs
        my $filetest = $seq;
        chomp($filetest);
        
        if (-s $filetest) {
            unless (open(DAT,$filetest)) {
                _error("Can't read table file $filetest: $!");
                return;
            } # can't open DAT
            # 6 lines @ 80 char/line or enough of a tbl, raw or cdna seq to be fairly sure
            read(DAT,$tmparray[0],480);
            my $guess = $$self{cbc}->check_type(\@tmparray,%params);
            if ($guess eq "UNKNOWN") { _error("Can't determine sequence format from $tmparray[0]") }
            else {
                tie @$AR_seqs,"CompBio::Simple::ArrayFile",$filetest,$guess;
            } # tie array to file
            $return_type = "ARRAY";
        } # we have a file to read
        # maybe I should tie array differently and not copy whole 
        # data set split into array?
        else {
            $new_ref = $ref ? $seq : \$seq;
            $return_type = $ref ? "REFSCALAR" : "SCALAR";
            $tmparray[0] = join("\n",split(/\n/,$$new_ref,6));
            my $guess = $$self{cbc}->check_type(\@tmparray,%params);
            $params{FORMAT} = $guess;
            if ($guess =~ /TBL|RAW|CDNA/) {
                $AR_seqs = [(split(/\n/,$$new_ref))];
            } # doesn't matter which, they are all single 'line' records
            elsif ($guess eq "FA") {
                $AR_seqs = [(split(/\n?(?=>)/,$$new_ref))];
                $$AR_seqs[-1] =~ s/\s+$//; # for some reason chomp($$AR_seqs[-1]) didn't work?
            } # split appropriate to fasta
            elsif ($guess eq "IG") {
                $AR_seqs = [(split(/\n?(?<=1\n)(?=;)/,$$new_ref))];
                $$AR_seqs[-1] =~ s/\s+$//; # for some reason chomp($$AR_seqs[-1]) didn't work?
            } # split appropriate to ig
            else { _error("Got $guess as format type for sequence(s)\n") }
        } # else examine the scalar and populate ar_seqs
    } # not 
    elsif ($ref) { _error("Got an unusable reference ->$ref\n") }
    else {
        _error("Can't determine how to handle $seq\n");
        return;
    } # else

    if ($params{'OUTFILE'}) {
        my $cat = ">";
        if (-s $params{'OUTFILE'}) {
            print "Requested output file $params{'OUTFILE'} exists, overwrite(1), append(2) or quit(3)? [2] ";
            chomp(my $input = <STDIN> || 2);
            if ($input == 2) {
                $cat = ">>";
            } # append
            elsif ($input == 1) { } # just fall out
            elsif ($input == 3) {
                _error("User canceled operation.\n");
                return;
            } # quit operation
            else {
                _error("Don't know what to do given $input\n");
                return;
            } # huh?
        } # if file already exists
        
        unless (open(OUT,"$cat $params{OUTFILE}")) {
            _error("Can't open OUT $params{OUTFILE}: $!");
            return;
        } # unless OUT opens
        $return_type = "FILE";
        $params{'FH'} = *OUT;
    } # output to file requested
    
    $return_type = $params{RETURN_TYPE} || $return_type;
    
    return($self,$AR_seqs,$return_type,%params);
} # _munge_seq_input

sub _munge_seq_return {
    my $AR_seqs = shift;
    my $return_type = shift;
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my $items = ref($AR_seqs) ? scalar(@$AR_seqs) : 1;

    print "_munge_seq_return recieved $items entries for return type $return_type\n" if $DEBUG;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }  
    
    return $AR_seqs unless (ref($AR_seqs) && $return_type ne "ARRAY");
    
    if ($return_type eq "SCALAR") {
        my $ret = join("\n",@$AR_seqs) . "\n";
        return $ret;
    } # if
    elsif ($return_type eq "FILE") {
        my $fh = $params{'FH'}; # print doesn't like this being used directly
        print $fh join("\n",@$AR_seqs) , "\n";
        close $fh or _error("$params{'OUTFILE'} failed to close properly: $!");
        return 0;
    } # elsif
    elsif ($return_type eq "HASH") {
        my %h = ();
        foreach (@$AR_seqs) {
            my ($id,$seq) = split;
            $h{$id} = $seq;
        } # foreach
        return \%h;
    } # elsif
    elsif ($return_type eq "2D") {
        foreach (@$AR_seqs) {
            $_ = [(split)];
        } # foreach
        return $AR_seqs;
    } # elsif
    elsif ($return_type eq "REFSCALAR") {
        my $ret = join("\n",@$AR_seqs) . "\n";
        return \$ret;
    } # elsif
    
    return $AR_seqs;
} # _munge_seq_return

sub AUTOLOAD {
    my $program = $AUTOLOAD;
    $program =~ s/.*:://;
    print "AUTOLOAD being used to call $program\n" if $DEBUG >= 3;

    if ($Auto_method{$program}) {
        # first, because this is Simple, we'll make it handle seq data input
        my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
        my $AR_ret = $$self{cbc}->$program($AR_seqs,%params);
        # And then we'll munge the output back to format provided
        return _munge_seq_return($AR_ret,$return_type,%params);
    } # call is for a a basic core method that takes a set of sequences in an array
} # AUTOLOAD

sub _error {
    my $msg = shift || "";
    my $die = shift || "";
    if (($die && $DEBUG) || $DEBUG >= 5) {
        confess $msg;
    } # if
    elsif ($die || $DEBUG >= 4) {
        croak $msg;
    } # elsif
    elsif ($DEBUG >= 2) {
        cluck $msg;
    } # elsif
    else { carp $msg }
} # _error - throw $DEBUG dependant exception

1;

__END__

=head2 check_type

Checks a given sequence or set of sequences for it's type. Currently groks 
fasta(.fa), table(.tbl), raw genome(.raw), intelligenics[?](.ig) and coding
dna sequence(.cdna) types. Each index of the referenced array should
be an entire sequence record. * It is however B<not> recomended that you load
up an entire raw genome into memory to do this test - see L<perlfunc read> *

C<$type = check_type(\@seqs,%parameters);>
        
Possible return types are CDNA, TBL, FA, IG, RAW and UNKNOWN.

OPTION:

CONFIDENCE: Set this parameter to an positive integer representing the number of
checks against your sequence set you would like performed. The default value is 3,
however if only one or two sequences are in your set only one check will be performed.
If CONFIDENCE is set to a value equal to or greater than the number of sequences in
your set, it will be reset to a value approx. 66% of the number of sequences in your
set (a I<very> high level of confidence).

Check type will then make CONFIDENCE # of tests on youe sequence set, randomly
selecting 2 members into a subset and checking the format. If more than one format
type guesses check_type will throw an exception reporting all the types found.

=head2 tbl_to_fa

Converts sequence records in table (tab delimited, usually .tbl file extension)
format to fasta format.

C<$fa_seq = $cbs->tbl_to_fa($tblseq,%params);>

Extra data fields in the table format will be added to the annotation line (>)
in the fasta output.

=head2 tbl_to_ig

Converts sequence records in table (tab delimited) format to .ig format.

C<$aref_igseqs = $cbc->tbl_to_ig(\@tbl_seqs,%params);>

Extra data fields in the table format will be placed on a single comment
line (;) in the ig output.

=head2 fa_to_tbl

Converts sequence records in fasta format to table format.

C<$aref_faseqs = $cbc->fa_to_tbl(\@fa_seq);>

Extra data in the fasta format will be placed in an extra tab delimited
field after the sequence record.

Options:

CLEAN: Reduces signifier to the first strech of non white space characters
found at least 4 characters long with none of the characters (| \ ? / ! *)

=head2 ig_to_tbl

Accepts ig format sequences in a referenced array, one complete sequence
record per index. This method returns the sequence(s) in table(.tbl) format
contained in a referenced array.

C<$aref_igseqs = $cbc->ig_to_tbl(\@fa_seq);>

Extra comment lines in the ig format will be placed as extra tab delimited
values after the sequence record.

=head2 dna_to_aa

Convert dna sequences, containing no whitespace, to the amino acid residues
as coded by the standard 'universal' genetic code. dna sequence may contain
standard special characters (i.e. R,S,B,N, ect.).

C<$aa = dna_to_aa(\$dna_seq,%params);>

Options:

C: Set to a true value to indicate dna should be converted to it's complement
before translation.

ALTCODE: A reference to a hash containing alternate coding keys where the value 
is the new aa to code for. Stop codons are represented by ".".

SEQFIX: Set to true to alter first position, making V or L an M, and
removing stop in last position.

=head2 complement

Converts dna to it's complementary strand. DNA sequence is submitted as scalar
by reference. There is no return as sequence is modified. To maintain original
sequence, send a reference to a copy.

C<complement(\$dna);>

=head2 six_frame

C<$result = six_frame($raw_file,$id,$seq_len,$out_file);>
	
Converts a submitted dna sequence into all 6 frame translations to aa.
Output id's have start and stop positions encoded; if first value is larger,
translated from anti-sense.

Options:

ALTCODE: A reference to a hash containing alternate coding keys where the value 
is the new aa to code for. Stop codons are represented by ".". Multiple
allowable values for the final dna in the codon may be provided in the
format AT[TCA].

ID: Prefix for translated segment identifiers, default is SixFrame.

SEQLEN: Minimum length of aa sequence to return, default is 10.

OUTPUT: A filename to pipe results to. This is recomended for large dna
sequences such as large contigs or whole chromosomes, otherwise results are
stored in memory until process complete ** this includes the use of the general
method option of OUTFILE ** . If the value of OUTFILE is STDOUT
results will be sent directly to standard out.


=head2 AA_HASH

Creates a hash table with amino acid lookups by codon. Includes all cases where
even an alternate na code (such as M for A or C) would return an unambiguous aa.
Also consistent with the complement method in this package, ie, lower cases in
some contexts, for ease of use with six_frame.

 C<%aa_hash = aa_hash;>

=head1 EXPORT

None by default.

=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC -n CompBio::Simple

=item 0.42

Got initial set of methods linking to CompBio in. Developed the _munge stuff.

=item 0.43

Brought tests up to same point as in CompBio. Added six_frame. Moved the data handeling
for dna_to_aa to own _munge type internal method. Added RETURN_TYPE parameter to
_munge input so user can pick. Created the CompBio::Simple::ArrayFile module to TIE
a sequence file to an array. Added six_frame using same array processing _munge as
dna_to_aa.

=item 0.44

Made seriouse attempt to get the docs caught up.
Added random sampling and CONFIDENCE to check_type.

=back

=head1 TO DO

_munge_seq_imput needs to also detect when it just has a list of loci and hand
request to DB to attempy to fulfil, and/or accept a DB param.

Look at using a tie type method similar to ArrayFile.pm for fetching data from
a database - simply throwing an exception if no db module provided (eval use?).

rewrite (<sigh>) _munges to detect submited format (fasta, etc) and return same
format whenever possible by default, and adding RETURN_FORMAT parameter, negating
the need for user to call converter after dna_to_aa for example. *This is probably
over the top since most methods _are_ the converters, but why make myself write
if/elsif block in all the utils that offer different format returns?*

=head1 COPYRIGHT

Developed at the BioMolecular Engineering Research Center at Boston
University under the NHLBI's Programs for Genomic Applications grant.

Copyright Sean Quinlan, Trustees of Boston University 2000-2001.

All rights reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Sean Quinlan, seanq@darwin.bu.edu

Please email me with any changes you make or suggestions for changes/additions.
Latest version is available under ftp://mcclintock.bu.edu/BMERC/perl/.

Thank you!

=head1 SEE ALSO

L<perl(1)>, L<CompBio(3)>.

=cut
