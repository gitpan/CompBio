package CompBio::Simple;

require 5.005_62;
use strict;
use warnings;
use Carp qw(cluck confess croak carp);
use CompBio;
use CompBio::Simple::ArrayFile;

our @ISA = qw(Exporter);
use vars qw($AUTOLOAD);

our $VERSION = '0.43';
our $DEBUG = 0;

our %Auto_method = (
    check_type => 1,
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

my $cbs = CompBio::Simple->new;

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
or ids accept input as scalar, arrays (by reference) or hashes (by reference).
If I miss documenting it, try it just in case. Unless a parameter option is
available to signify otherwise, I assume 'doing the right thing' is to return the
same type of data construct as recieved.

Sorry, not all of the methods available trough Simple.pm have been documented.
See CompBio.pm for details on the basic usage and purpose of the methods. I'm
planning on finishing up the docs to the current state for 0.46.

=head2 new

Construct an object for invoking methods in CompBio::Simple.

=cut
sub new($%) {
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

=head2 help

Quits current application and uses perldoc to display the POD.

=cut
sub _help() {
    exec 'perldoc $0';
} # perldoc cheat

=head2 tbl_to_fa

Converts a sequence record in table (tab delimited, usually .tbl file extension)
format to fasta format.

	C<$fa_seq = $cbs->tbl_to_fa($tblseq,%params);>

$seqdat can be either scalar, array by reference, 2d array by reference, or a
hash by reference; return data type is same as submitted. Scalar may be either sequence
records (one record per line), or a filename. Array's should contain an entire
sequence record (loci\tsequence) in each indexed position. 2D arrays should be
ids in the first column, sequences in the second, such that $array[2][1] would
be the aa sequence for the third record. Hash's should have the loci as the key,
and only the sequence as the value.

=cut
# man this error checking shit is a pain! :P
# get the rest of the methods docs from CompBio.pm
sub six_frame {
    my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
        
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($AR_seqs),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    my $AR_ret = _munge_array_to_scalar('six_frame',$self,$AR_seqs,%params);
    return _munge_seq_return($AR_ret,$return_type,%params);
} # six_frame

sub complement {
} # complement

# why aren't I checking format type and converting as necisary?!?
# simple should DTRT here! Think about data storage situations where
# converting to array will ruin type formating - and data. And where
# check_type will fail on sequence format incorrect. This should end
# being part of sanity shecking as well. Make sure I structure to fall
# as fast as possible.

sub dna_to_protein {
    my ($self,$AR_seqs,$return_type,%params) = _munge_seq_input(@_);
    
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    if ($DEBUG >= 3) { foreach my $k (keys %$self) { print "$k\t$$self{$k}\n" } }
    print ref($AR_seqs),"<- ref?\n" if $DEBUG >= 2;
    if ($DEBUG >= 2) { foreach my $k (keys %params) { print "$k\t$params{$k}\n" } }
    
    my $AR_ret = _munge_array_to_scalar('dna_to_protein',$self,$AR_seqs,%params);
    return _munge_seq_return($AR_ret,$return_type,%params);
} # dna_to_protein

sub _munge_array_to_scalar {
    my $caller = shift; # can I get this from caller or such?
    my $self = shift;
    my $AR_seqs = shift;
    my %params = @_ ? @_ : ();
    
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
        foreach my $seq (@$AR_seqs) {
            my ($id,$dna) = split($seq);
            my $ret = $id . ${$$self{cbc}->$caller(\$dna,%params)};
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
        
        if (-s $seq) {
            unless (open(DAT,$seq)) {
                _error("Can't read table file $seq: $!");
                return;
            } # can't open DAT
            # 6 lines @ 80 char/line or enough of a tbl, raw or cdna seq to be fairly sure
            read(DAT,$tmparray[0],480);
            my $guess = $$self{cbc}->check_type(\@tmparray,%params);
            if ($guess eq "UNKNOWN") { _error("Can't determine sequence format") }
            else {
                tie @$AR_seqs,"CompBio::Simple::ArrayFile",$guess;
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
        unless (open(OUT,"> $params{OUTFILE}")) {
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
    if (($die && $DEBUG) || $DEBUG >= 4) {
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
for dna_to_protein to own _munge type internal method. Added RETURN_TYPE parameter to
_munge input so user can pick. Created the CompBio::Simple::ArrayFile module to TIE
a sequence file to an array. Added six_frame using same _munge as dna_to_protein.

=back

=head1 TO DO

_munge_seq_imput needs to also detect when it just has a list of loci and hand
request to DB to attempy to fulfil, and/or accept a DB param.

check_type invoked from here should make a few type checks with random samples
from the dataset when possible.

=head1 COPYRIGHT

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
