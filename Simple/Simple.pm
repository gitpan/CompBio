package CompBio::Simple;

require 5.005_62;
use strict;
use warnings;
use Carp qw(cluck confess croak carp);
use CompBio;

our @ISA = qw(Exporter);
use vars qw($AUTOLOAD);

our $VERSION = '0.42';
our $DEBUG = 0;

our %Auto_array = (
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
} # six_frame

sub complement {
} # complement

sub dna_to_protein {
} # dna_to_protein

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
    elsif ($ref eq "SCALAR") {
        $AR_seqs = [split(/\n/,$$seq)];
        $return_type = "REFSCALAR";
    } # if scalar reference submitted (!?)
    elsif ($seq =~ /^.+[\t ][A-Z\.\*\!]+$/) {
        $AR_seqs = [(split(/\n/,$seq))];
        $return_type = "SCALAR";
    } # elsif first line looks like scalar sequences in table format
    elsif ($ref) { _error("Got an unusable reference ->$ref\n") }
    elsif (-s $seq) {
        unless (open(TBL,$seq)) {
            _error("Can't read table file $seq: $!");
            return;
        } # can't open TBL
        while (<TBL>) {
            chomp;
            my @fields = split; # in case some used spaces instead of tab
            push(@$AR_seqs,join("\t",@fields[0,1]));
        } # while
        close TBL or _error("TBL didn't close\n");
        $return_type = "ARRAY"; # since this gives no indication of preference
    } # eslif $seq appears to be a filename
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
    
    return($self,$AR_seqs,$return_type,(each %params));
} # _munge_seq_input

sub _munge_seq_return {
    my $AR_seqs = shift;
    my $return_type = shift;
    my %params = @_ ? @_ : ();
    
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

    if ($Auto_array{$program}) {
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

Copy over functions from original BMERC::bio, making small improvements to code,
mostly by removing lingering locale assumptions and (hopefully) improving
interface, and converting to OOP.

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
