package CompBio::Simple::ArrayFile;

require 5.005_62;
use strict;
use warnings;
use Carp;

require Exporter;

#our @ISA = qw(Exporter);
#our %EXPORT_TAGS = ( 'all' => [ qw( ) ] );
#our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
#our @EXPORT = qw( );

our $VERSION = '0.01';

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

sub FETCH {
    my $self = shift;
    confess "I ($self) am not a class method!" unless ref $self;
    my $fh = $$self[0];
    my $format = $$self[1];
    my $seqrec = "";
    
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

sub DESTROY {
    my $self = shift;
    confess "I ($self) am not a class method!" unless ref $self;
    my $fh = $$self[0];
    close $fh or carp "Couldn't properly close tied file $$self[1]: $!";
} # close filehandle

1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

CompBio::Simple::ArrayFile - Perl extension for blah blah blah

=head1 SYNOPSIS

  use CompBio::Simple::ArrayFile;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for CompBio::Simple::ArrayFile, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC
	-n
	CompBio::Simple::ArrayFile

=back


=head1 AUTHOR

A. U. Thor, a.u.thor@a.galaxy.far.far.away

=head1 SEE ALSO

perl(1).

=cut
