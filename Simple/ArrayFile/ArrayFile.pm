package CompBio::Simple::ArrayFile;

require 5.005_62;
use strict;
use warnings;
use Carp;
use Tie::Array;
our @ISA = ('Tie::Array');
require Exporter;

our $VERSION = '0.20';

=head1 NAME

CompBio::Simple::ArrayFile - provide access to a sequence file through an array

=head1 SYNOPSIS

  use CompBio::Simple::ArrayFile;
  tie @$AR_seqs,"CompBio::Simple::ArrayFile","format";
  # where format is as defined in CompBio::check_type

=head1 DESCRIPTION

This allows an input file of sequence data to be tied to the array returned by
CompBio::Simple::_munge_seq_input. It provides a way to have all the processing
code use one format without needing to read the whole file into memory.

Please check your favorite perl reference for more information on the tie function.

This tie module provides only three methods, the TIEARRAY (basically a new in OO
syntax) to create the tied array, FETCH to read a sequence record from the file,
and DESTROY which basically just closes the filehandle;

=head2 EXPORT

None by default.

=cut

sub TIEARRAY {
    my $proto = shift;
    my $class = ref($proto) || $proto; # is this worthwhile in the context of tie?
    my $filename = shift;
    my $seqformat = shift;
    
    unless (open(IN,$filename)) {
        carp "Can't tie $filename: $!";
        return;
    } # warn and return void if we can't read
    my $fh = *IN;
    
    return bless [($fh,$seqformat,$filename)], $class;
} # create the tied array

# maybe I should just try to write a temp file
# in table format first?
my $last_call = ""; # kludge to reset position
sub FETCH {
    my $self = shift;
    confess "I ($self) am not a class method!" unless ref $self;
    my $fh = $$self[0];
    my $format = $$self[1];
    my $seqrec = "";
    my $call = join(":",((caller(1))[3,2]));
    if ($last_call ne $call) {
        $last_call = $call;
        seek $fh,0,0;
    } # called from different routine - reset
    
    if ($format =~ /TBL|RAW|CDNA/) {
        $seqrec = <$fh>;
    } # single line records
    elsif ($format eq "FA") {
        my $pos = tell $fh;
        my @lines;
        chomp($seqrec = <$fh>);
        push(@lines,$seqrec);
                
        while (<$fh>) {
            chomp;
            if (/^>/ && $lines[-1] !~ /^>/) {
                seek $fh,$pos,0;
                last;
            } # we've reached the next record
            push(@lines,$_);
            $pos = tell $fh;
        } # get lines until we drop out at eof, or we see the next > indicator
        $seqrec = join("\n",@lines);
    } # get next fasta record
    elsif ($format eq "IG") {
        my $pos = tell $fh;
        my @lines;
        chomp($seqrec = <$fh>);
        push(@lines,$seqrec);
                
        while (<$fh>) {
            chomp;
            if (/^;/ && $lines[-1] !~ /1\s*$/) {
                seek $fh,$pos,0;
                last;
            } # we've reached the next record
            push(@lines,$_);
            $pos = tell $fh;
        } # get lines until we drop out at eof, or we see the next ; indicator
        $seqrec = join("\n",@lines);
    } # get next ig record
} # read data

sub FETCHSIZE {
    my $self = shift;
    confess "I ($self) am not a class method!" unless ref $self;
    my $fh = $$self[0];
    my $format = $$self[1];
    
    if ($format =~ /TBL|RAW|CDNA/) {
        my $pos = tell $fh;
        seek $fh,0,0;
        my $c = 0;
        while (<$fh>) { $c++ }
        seek $fh,$pos,0;
        return $c;
    } # single line records
    elsif ($format eq "FA") {
        my $pos = tell $fh;
        seek $fh,0,0;
        my $c = 0;
        while (my $l = <$fh>) { $c++ if $l =~ /^>/ } # will break on concurrent > lines
        seek $fh,$pos,0;
        return $c;
    } # elsif fasta
    elsif ($format eq "IG") {
        my $pos = tell $fh;
        seek $fh,0,0;
        my $c = my $on = 0;
        while (<$fh>) {
            if (/^\s*;/) {
                $on = 1;
            } # comment/seperator lines found
            elsif ($on) {
                $c++;
                $on = 0;
            } #  should be id line after comment/seperator lines
        } # while counting records
        seek $fh,$pos,0;
        return $c;
    } # elsif ig
} # FETCHSIZE


sub DESTROY {
    my $self = shift;
    confess "I ($self) am not a class method!" unless ref $self;
    my $fh = $$self[0];
    close $fh or carp "Couldn't properly close tied file $$self[1]: $!";
} # close filehandle

1;
__END__

=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC -n CompBio::Simple::ArrayFile

=head1 TO DO

Take a close look at whether or not the method for adding sequences should
be developed.

Consider ways to emulate array slices and specific index selects. Just
added random member testing for check_type and I'm sure this will break
when used there.

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

L<perl(1)>, L<CompBio(3)>, L<CompBio::Simple>.

=cut
