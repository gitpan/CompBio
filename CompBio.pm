package CompBio;

require 5.005_62;
use strict;
use vars qw($AUTOLOAD);

require Exporter;

our @ISA = qw(Exporter);

our $CONF_FILE = $ENV{COMPBIO_CONF} || '/etc/CompBio.conf';
our $CPUSERVER = $ENV{HOST};

our %EXPORT_TAGS = ( 'all' => [ qw(check_type tbl_to_fa tbl_to_ig fa_to_tbl ig_to_tbl
    dna_to_aa complement six_frame aa_hash wu_blast) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw($CONF_FILE $CPUSERVER);

our $VERSION = '0.469';
our $DEBUG = 0;
our $RET_CODE = sub {
        my $ar_ret = shift;
        my $ret = shift;
        print "ref?->$ar_ret\nseq_in->$ret\n" if $DEBUG >= 3;
        push(@$ar_ret,$ret);
    }; # default return method

=head1 NAME

CompBio - Core library for some basic methods useful in computational biology/bioinformatics.

=head1 SYNOPSIS

use CompBio;
  
my $cbc = new->CompBio;

$AR_faseqs = $cbc->tbl_to_fa(\@tblseqs);

or

use CompBio qw(tbl_to_fa);

$AR_faseqs = tbl_to_fa(\@tblseqs);

=head1 DESCRIPTION

The CompBio module set is being developed as a new implementation of the code
base originally developed at the BioMolecular Engineering Research Center
(http://bmerc-www.bu.edu). CompBio.pm is intended to take a number
of small, commonly used methods for supporting bioinformatics research, and
make them into a single package. Many of the utils included with this
distribution are just command line interfaces to the methods contained herein.

The CompBio module set is _not_ intended to replace the bioperl project
(http://www.bioperl.org/). Although I do welcome suggestions for improving
or adding to the methods available in these modules, particularly I would
love any help with things on the TO DO list, these modules are not intended
to provide the depth that the bioperl suite can provide.

CompBio has a limited API. It expects it's input to be in specific
formats, as described in each methods description, and it's output
is in a format that makes the most sense to me for that method. It does
no error checking by and large, so incorrect input could cause bizzare
behavior and/or a noisy death. If you want a flexible interface with lots
of error checking and deep levels of verbosity, use L<CompBio::Simple> - that's
its job.

Latest version is available through SourceForge
http://sourceforge.net/projects/compbio/, or on the CPAN
http://www.cpan.org/authors/id/S/SE/SEANQ/CompBio-0.461.tar.gz.

Thanks!

Other modules available (or that will be available) in the CompBio set are:

The DB module will only be imediately useful if you import the databases as
used by us here at the BMERC(ftp://mcclintock.bu.edu/BMERC/mysql/) or develop
your own on the same basic design scheme. Otherwise I hope you find it useful
as a source of design ideas for rolling your own. Please note however that I
intend to expand the methods in CompBio/Simple.pm to allow seemless access to
data through the DB module. Although at no time will including the DB module be
required for Simple to work, I think that if you have the space for it, you will
find having the data locally and adding this module (or adapting it to work with
databases already installed) will have an imense and imediate benifiacial
impact. It did for us at least! :)

The Profile module was designed to work with our PIMA-II sofware and the
PIMA modules. The PIMA suite is available for license from Boston University
for a nominal fee, and free for academic use. For examples and more info
see http://bmerc-www.bu.edu/PIMA/. Unless you have or are interested in that
package, this module will have no functional value.

=head1 Methods

You may note that the majority of the methods here are for converting
sequences from one format to another. Mainly this is for converting other
formats to table format, which is used by most of these programs. This is
not meant to be a comprhensive collection of format guessing and
transformation methods. If you are looking for a converter for a format
CompBio doesn't handle, I suggest you look into bioperls SeqIO package or
the READSEQ program (java), found at
http://iubio.bio.indiana.edu/soft/molbio/readseq/

=head2 new

Construct an object for invoking methods in CompBio.

=cut
sub new {
    my ($proto,%parameters) = @_;
    my $class = ref($proto) || $proto;
    my $self = {};

    #handle params as nec. such as setting debug or changing env. variables
    $DEBUG = $parameters{'DEBUG'} || 0;
    (require warnings && import warnings) if $DEBUG; #??
    (require diagnostics && import diagnostics) if $DEBUG >= 2; #??
    $self->{'_created'} = 1;
    print "Created = ",$self->{'_created'},"\n" if $DEBUG >= 3;
    
    return bless ($self,$class);
} # new

=head2 help

Quits current application and uses perldoc to display the POD.

=cut
sub _help() {
    exec 'perldoc $0';
} # perldoc cheat


=head2 check_type

Checks a given sequence or set of sequences for it's type. Currently groks 
fasta(.fa), table(.tbl), raw genome(.raw), intelligenics[?](.ig) and coding
dna sequence(.cdna) types. Each index of the referenced array should
be an entire sequence record. * It is however B<not> recomended that you load
up an entire raw genome into memory to do this test - see L<perlfunc read> *

C<$type = check_type(\@seqs,%parameters);>
        
Possible return types are CDNA, TBL, FA, IG, RAW and UNKNOWN.

Be warned, this is intended only as a quick check and only uses as many
records as necisarry in the reference provided, stating with the first.
check_type assumes the rest of the records look the same and does not do any
kind of deep QA on the set. If you are not sure, invoke check_type with
a few random samples from your set, or use L<CompBio::Simple>, which does
that by default.

=cut
sub check_type {
    my $self = shift;
    my $seq = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($seq) && $seq eq "help");
    _error("Not given an array reference!\n",1) unless ref($seq) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    print "In check_type\n" if $DEBUG >= 2;
    print "Checking type for:\n",join("\n",@$seq[0..2]),"\n" if $DEBUG >= 3;
    
    # TBL should only accept aa codes
    # TBL and CDNA listed as multiline because of how sequences loaded to
    # test a file - given the way these are written now, should this even be an
    # array - none use $$seq[1]. Or keep just for format compat?
    if ($$seq[0] =~ /^.+?\t[CTUAGctuag]+\*?(\t|\n|$)/m) { return "CDNA" }
    elsif ($$seq[0] =~ /^.+?\t[A-Za-z!\.]+\*?(\t|$)/m) { return "TBL" }
    elsif ($$seq[0] =~ /^>.+\n[A-Za-z!\.]+\*?$/m) { return "FA" }
    elsif ($$seq[0] =~ /^(;.*\n)+\S+\n[A-Za-z!\.\*]+1?$/m) { return "IG" }
    elsif ($$seq[0] =~ /^[CTUAGMRWSYKVHDBXNctuagmrwsykvhdbxn]+\*?$/m) { return "RAW" }
    else { return "UNKNOWN" }
} # check_type

=head2 tbl_to_fa

Converts a sequence record in table (tab delimited, usually .tbl file extension)
format to fasta format.

C<$aref_faseqs = $cbc->tbl_to_fa(\@seqdat,%params);>

Each index in the @seqdat array must contain entire record (loci\tsequence) for
single sequence. Return is an array reference, still one sequence per index.

Extra data fields in the table format will be added to the annotation line (>)
in the fasta output.

=cut
sub tbl_to_fa {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    print "In tbl_to_fa\n" if $DEBUG >= 2;
    my @ret = ();

    foreach (@$aref_seqs) {
        my $seq = $_;
        chomp $seq;
        # using $seq is complete black magic. Somehow FETCH was being called
        # for ArrayFile when this split occurs, but not using the aliased $_
        # aliviates this. <SIGH> Someday I hope someone will explain this to me
        my @fields = split(/\t/,$seq);
        print "ID: $fields[0]\nSequence: $fields[1]\n" if $DEBUG >= 2;
        my $str = ">$fields[0]";
        if (@fields > 2) { $str .= " " . join(" ",@fields[2 .. $#fields]) }
        # generate fasta sequence lines at 80 char per line (including newline)
        my $tmpl = "a79" x ((length($fields[1])/79) + 1);
        $str .= "\n" . join("\n",(unpack($tmpl,$fields[1])));
        $str =~ s/\n$//;
        &$RET_CODE(\@ret,$str);
    } # foreach submitted sequence
    return \@ret if @ret;
    return 0;
} # tbl_to_fa

=head2 tbl_to_ig

Converts a sequence in table (tab delimited) format to .ig format. Accepts
sequences in a referenced array, one record per index.

C<$aref_igseqs = $cbc-&gt;tbl_to_ig(\@tbl_seqs,%params);>

Extra data fields in the table format will be placed on a single comment
line (;) in the ig output.

=cut
sub tbl_to_ig {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    print "In tbl_to_ig\n" if $DEBUG >= 2;
    my @ret = ();
    
    foreach (@$aref_seqs) {
        my $seq = $_;
        chomp $seq;
        my @fields = split(/\t/,$seq);
        my $str = ";";
        if (@fields > 2) { $str .= "\n;" . join("\n;",@fields[2 .. $#fields]) }
        $str .= "\n$fields[0]\n"; # seperators, comment & id lines
        $fields[1] .= "1";
        my $tmpl = "a79" x ((length($fields[1])/79) + 1);
        $str .= join("\n",(unpack($tmpl,$fields[1])));
        $str =~ s/\n$//;
        &$RET_CODE(\@ret,$str);
    } # foreach submitted sequence

    return \@ret if @ret;
    return 0;
} # tbl_to_ig

=head2 fa_to_tbl

Accepts fasta format sequences in a referenced array, one complete sequence
record per index. This method returns the sequence(s) in table(.tbl) format
contained in a referenced array.

C<$aref_faseqs = $cbc->fa_to_tbl(\@fa_seq);>

Extra data fields in the fasta format will be placed as extra tab delimited
values after the sequence record.

Options:

CLEAN: Reduces signifier to the first strech of non white space characters
found at least 4 characters long with no of the characters (| \ ? / ! *)

=cut
sub fa_to_tbl {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    print "In fa_to_tbl\n" if $DEBUG >= 2;
    my @ret = ();

    # traverse referenced array and convert
    foreach (@$aref_seqs) {
        #my $seqrec = $_; # to avoid calling non-existant STORE on tied array
        chomp;
        my $tbl = my $rem = "";
        
    	foreach (split(/[\n\r]+/)) {
    	    my $test = $_;
            if ($test =~ /^\s*>(\S+)\s*(.*)/) {
                my $sig = $1;
                print "sig = $sig\n" if $DEBUG >= 3;
                $rem = $2;
                if ($sig =~ s/^\w+\|([^\|\s:;.]+)(.+)/$1/) {
                    # Since this appears to be a genbank type id portion, try
                    # to isolate the second field, which is usually the gi#
                    print "sig2 = $sig\n" if $DEBUG >= 3;
                    my $temp = $2;
                    $temp =~ s/\|+/\t/g;
                    $rem .= " [$temp]";
                } # turn those annoying genbank pipes into tabs
                
                # maybe a better keyword than CLEAN?
                if ($params{'CLEAN'} && $sig =~ /([^\s\|\\?\/!*]{4,})/) {
                    $sig = $1;
                } # elsif
                print "In fa_to_tbl using $sig as tbl id field\n" if $DEBUG >= 2;
    	        $tbl = "$sig\t";
    	    } # if id line

    	    else {
	        $test =~ s/[\*1\!\.]$//;
    	        $tbl .= $test;
    	    } # else
        } # foreach line in seqrec
        $tbl .= "\t$rem" if $rem;
        &$RET_CODE(\@ret,$tbl);
    } # foreach seqrec
    
    return \@ret if @ret;
    return 0;
} # fa_to_tbl

=head2 ig_to_tbl

Accepts ig format sequences in a referenced array, one complete sequence
record per index. This method returns the sequence(s) in table(.tbl) format
contained in a referenced array.

C<$aref_igseqs = $cbc->ig_to_tbl(\@fa_seq);>

Extra comment lines in the ig format will be placed as extra tab delimited
values after the sequence record.

=cut
sub ig_to_tbl {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    my @ret = ();
    print "In ig_to_tbl\n" if $DEBUG >= 2;

    foreach (@$aref_seqs) { # traverse referenced array and convert
    	chomp;
        my $tbl = "";
        my @comments = ();
        
    	foreach (split(/[\n\r]+/)) {
            next if /^;\s*$/;
            if (/^;\s*(.+)$/) {
                push(@comments,$1);
            } # get comments
            elsif (! $tbl) {
                /^\s*([^\t]+)/;
                my $sig = $1;
                # maybe a better keyword than CLEAN?
                if ($params{'CLEAN'} && $sig =~ /(\S+)/) {
                    $sig = $1;
                } # if user want and we can, get a better id
                # Also, may want to add to bad_characters in []
                elsif ($params{'REALCLEAN'} && $sig =~ /(\S{3,})[\|\!*\-]/) {
                    $sig = $1;
                } # elsif
    	        $tbl .= "$sig\t";
            } # if first non comment line
            else {
                $tbl .= $_;
            } # we assume its a sequence line
        } # foreach line in seqrec
        $tbl =~ s/1$//;
        $tbl = join("\t",($tbl,@comments));
        $tbl =~ s/\s+$//;
        &$RET_CODE(\@ret,$tbl);
    } # foreach seqrec

    return \@ret if @ret;
    return 0;
} # ig_to_tbl

=head2 dna_to_aa

Converts a dna sequence, containing no whitespace, submited as a scalar reference,
to the amino acid residues as coded by the standard 'universal' genetic code.
Return is a reference to a scalar. dna sequence may contain standard special
characters (i.e. R,S,B,N, ect.).

C<$aa = dna_to_aa(\$dna_seq,%params);>

Options:

C: Set to a true value to indicate dna should be converted to it's complement
before translation.

ALTCODE: A reference to a hash containing alternate coding keys where the value 
is the new aa to code for. Stop codons are represented by ".". Multiple
allowable values for the final dna in the codon may be provided in the
format AT[TCA].

SEQFIX: Set to true to alter first position, making V or L an M, and
removing stop in last position.

=cut
sub dna_to_aa {
    my $self = shift;
    my $sref_seq = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    _help() if (! ref($sref_seq) && $sref_seq eq "help");
    _error("Not given a scalar reference!\n",1) unless ref($sref_seq) eq "SCALAR";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    my @ret = ();

    my %AA = aa_hash();
    my $dna = uc($$sref_seq);
    $dna =~ tr/U/T/;
    print "in dna_to_aa using:\n",$dna,"\n" if $DEBUG > 2;
    
    if($params{'ALTCODE'}) {
    	while ((my $codon,my $aar) = each %{$params{'ALTCODE'}}) {
    	    if ($codon =~ /^(\w\w)\[([A-Z]+)\]$/) {
                foreach (split("",$2)) { $AA{"$1$_"} = $aar }
            } # multiple options in third position
            else { $AA{$codon} = $aar }
    	} # for each alternate codon
    } # if

    if ($params{'C'}) { # convert all characters to thier complement
    	complement(\$dna);
#        die "Recieved C";
    } # if
    
    my $ret = "";
    # we allow any non-space character in match in case someone wants odd alternate coding
    while ($dna =~ /(\S{3})/g) {
        $ret .= $AA{$1} || "X";
    } # while translating sequence
    
    if ($params{'SEQFIX'}) {
        $ret =~ s/^[VL]/M/;
        $ret =~ s/\.$//;
    } # clean up ends for use as aa seq
    
    print "dna_to_aa produced $ret\n" if $DEBUG >= 3;
    return \$ret;
} # dna_to_aa

=head2 complement

Converts dna to it's complementary strand. DNA sequence is submitted as scalar
by reference. There is no return as sequence is modified. To maintain original
sequence, send a reference to a copy (or see documentation for 
CompBio::Simple::complement).

C<complement(\$dna);>

=cut
sub complement {
    my $self = shift;
    my $sref_seq = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    _help() if (! ref($sref_seq) && $sref_seq eq "help");
    _error("Not given a scalar reference!\n",1) unless ref($sref_seq) eq "SCALAR";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    my @ret = ();

    $$sref_seq =~ tr/[ACTUGMRYKVHDBkyrmbdhv]/[TGAACkyrmbdhvMRYKVHDB]/;
    $$sref_seq = reverse($$sref_seq); # should this uppercase before returning?
} # complement

=head2 six_frame

Converts a submitted dna sequence into all 6 frame translations to aa. Value
submitted may be a reference to a scalar containing the dna or a filename.
Either way dna must be in RAW(.raw) format (Not whitespace within sequence,
only one newline allowed at end of record). Output id's have start and stop
positions encoded; if first value is larger, translated from anti-sense.

C<$result = six_frame($raw_file,$id,$seq_len,$out_file);>
	
Options:

ALTCODE: A reference to a hash containing alternate coding keys where the value 
is the new aa to code for. Stop codons are represented by ".". Multiple
allowable values for the final dna in the codon may be provided in the
format AT[TCA].

ID: Prefix for translated segment identifiers, default is SixFrame.

SEQLEN: Minimum length of aa sequence to return, default is 10.

OUTPUT: A filename to pipe results to. This is recomended for large dna
sequences such as large contigs or whole chromosomes, otherwise results are
stored in memory until process complete. If the value of OUTFILE is STDOUT
results will be sent directly to standard out.

=cut
sub six_frame {
    my $self = shift;
    my $sref_raw = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    my $ref = ref($sref_raw);
    _help() if (! $ref && $sref_raw eq "help");
    _error("Not given a valid filename or scalar reference!\n",1)
        unless -e $sref_raw || $ref eq "SCALAR";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    my @ret = ();
    my %AA = aa_hash();

    if($params{'ALTCODE'}) {
    	while ((my $codon,my $aar) = each %{$params{'ALTCODE'}}) {
    	    if ($codon =~ /^(\w\w)\[([A-Z]+)\]$/) {
                foreach (split("",$2)) { $AA{"$1$_"} = $aar }
            } # multiple options in third position
            else { $AA{$codon} = $aar }
    	} # for each alternate codon
    } # if

    $params{'ID'} ||= "SixFrame";
    $params{'SEQLEN'} ||= 10;
    $params{'OUTFILE'} ||= "";
    my $fh_out = "";

# this section should now be hadled from Simple using RET_CODE
    if ($params{'OUTFILE'} eq "STDOUT") {
        $fh_out = *STDOUT;
    } # send to stdout
    elsif ($params{'OUTFILE'}) {
        my $out_type = $params{'CAT'} ? ">>" : ">";
    	open (OUT,"$out_type $params{'OUTFILE'}")
            or return "Cat: Can't open OUT $params{'OUTFILE'}: $!";
        $fh_out = *OUT;
    } # if concat onto existing file requested
    
    my $dna_len = my $pos = 0;
    my $read_code = "";
    
    if ($ref eq "SCALAR") {
        $read_code = sub { $$sref_raw =~ /\G(\w{3})/gc; return $1 };
        $dna_len = length($$sref_raw);
        die "$dna_len" unless $dna_len;
    } # else if cdna sequence
    elsif ($dna_len = -s $sref_raw) {
    	open (RAW,$sref_raw) or return "Can't open RAW $sref_raw to read: $!";
	$read_code = sub { read RAW,$_,3; return $_; };
    } # input is a raw file
    else { return "Don't understand input, or input empty\n" }

    my ($seq1,$seq2,$seq3,$cseq1,$cseq2,$cseq3,$ret);
    my $codon = &$read_code;
    $pos += 3;
    $seq1 = $AA{$codon} || "X";
    complement(\$codon);
    $cseq1 = $AA{$codon} || "X";
    my $last = $codon;

    my @codon;
    while ($codon = &$read_code) {
    	my $six = $last . $codon;
    	$last = $codon;
    	$pos+=3;

    	@codon = unpack "aa3X2a3",$six; # why get $_ frame into @codon?
    	$seq1 .= $AA{$codon} || "X";
        complement(\$codon);
        $cseq1 .= $AA{$codon} || "X";
        $seq2 .= $AA{$codon[1]} || "X";
        complement(\$codon[1]);
    	$cseq2 .= $AA{$codon[1]} || "X";
    	$seq3 .= $AA{$codon[2]} || "X";
        complement(\$codon[2]);
        $cseq3 .= $AA{$codon[2]} || "X";
    	    
    	my $out = _stop(\$seq1,$pos,$params{'SEQLEN'},$params{'ID'}) if $seq1 =~ /\.$/;
    	$out .= _stop(\$cseq1,$pos,$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq1 =~ /\.$/;
    	$out .= _stop(\$seq2,($pos-2),$params{'SEQLEN'},$params{'ID'}) if $seq2 =~ /\.$/;
    	$out .= _stop(\$cseq2,($pos-2),$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq2 =~ /\.$/;
    	$out .= _stop(\$seq3,($pos-1),$params{'SEQLEN'},$params{'ID'}) if $seq3 =~ /\.$/;
    	$out .= _stop(\$cseq3,($pos-1),$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq3 =~ /\.$/;
	
    	if ($params{'OUTFILE'} && $out) { print $fh_out $out }
	elsif ($out) { $ret .= $out }
    	undef($out);
    } # while reading dna
    
    my $out = "";
    $out .= _stop(\$seq1,$pos,$params{'SEQLEN'},$params{'ID'}) if $seq1;
    $out .= _stop(\$seq2,($pos-2),$params{'SEQLEN'},$params{'ID'}) if $seq2;
    $out .= _stop(\$seq3,($pos-1),$params{'SEQLEN'},$params{'ID'}) if $seq3;
    $out .= _stop(\$cseq1,$pos,$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq1;
    $out .= _stop(\$cseq2,($pos-2),$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq2;
    $out .= _stop(\$cseq3,($pos-1),$params{'SEQLEN'},$params{'ID'},$dna_len) if $cseq3;

    if ($out && $params{'OUTFILE'}) { print $fh_out $out }
    else { $ret .= $out }
    
    return $ret || 0; # this should really return a reference!
} # six_frame

# stop is used by six frame when a stop '.' is encountered to calculate the
# coding position and return the .tbl format id and sequence.
sub _stop {
    my $sref_seq = shift;
    my $loc = shift;
    my $min = shift;
    my $id = shift;
    my $raw_len = shift; # comp strand seq if provided
    my $start = my $end = my $ret = "";
    my $len = length($$sref_seq);
    if ($len < $min) { undef($$sref_seq);return }
    
    $$sref_seq =~ s/[X\.]$//;
    $$sref_seq =~ s/^\.// if $raw_len;
    
    if ($raw_len) {
    	$start = $loc;
    	$end = $loc - ($len * 3) + 1;
    } # if seq from complementary strand
    else {
    	$start = $loc - ($len * 3) + 1;
	$end = $loc;
    }
    
    $$sref_seq = reverse($$sref_seq) if $raw_len;
    $ret = $id . "_" . $start . "_" . $end . "\t$$sref_seq\n";
    $$sref_seq = "";
    return $ret;       
} # _stop

=head2 wu_blast

A simple interface to the Washington University distribution of blast. Only recently
tested with full 2.0MP-WashU release. I do however stick to using the compatability
features so it it should work with both the licensed and free releases of WUblast.

wu_blast requires two arguments, a sequence database filename or a reference to
an array containing the sequence(s) to search, and the query sequence filename or
a reference to an array containing the sequence to query with. Sequence files need
to be in fasta format and sequences contained in arrays need to be in table format.
For more versatility, use CompBio::Simple::blast. One thing wu_blast will do is
prepare the db for searching if it isn't set up in advance; this does require
you have write permision in the location of the base sequence file, and sufficient
room in your quota (if used).

Options:

METHOD: Alternate method to use for sequence comparison. Default is blastp.

PARAMS: Optional parameters to hand to the blast method. See the your local
documentation for the method chosen to see what options are available.

SERVER: Machine to use as the blast server. Defaults to the global $CPUSERVER.
rsh is currently used to launch the proccess.

BLASTPREP: Alternate method to use to prepare a sequence database for being
searched using the chosen blast method. Default is setdb.

C<$blast_out = blast($fa_file,$queryfile,(PARAMS =&gt; "-B=0 -matrix=BLOSUM62"));>

=cut
sub wu_blast {
    my $self = shift;
    my $db = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    my $refDb = ref($db);
    _help() if (! $refDb && $db eq "help");
    _error("Not given a valid filename or scalar reference for database!\n",1)
        unless -e $db || $refDb eq "ARRAY";
    
    my $query = shift;
    my $refQuery = ref($query);
    _error("Not given a valid filename or scalar reference for database!\n",1)
        unless -e $query || $refQuery eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    local $RET_CODE = exists $params{RET_CODE} ? $params{RET_CODE} : $RET_CODE;
    # change defaults set here to suit local preferences
    # make sure BLASTPREP is set appropriately if you change default method
    my $method = $params{METHOD} || "blastp";
    my $parameters = $params{PARAMS} || '';
    local $CPUSERVER = exists $params{SERVER} ? $params{SERVER}  : $CPUSERVER;
    # change this if your blast programs reside somewhere else
    my $base_dir = "/usr/local/bin/";
    # if they are in multiple locations, good luck!
    my $result = '';
    
#    $ENV{"BLASTFILTER"} = "/usr/local/ncbi/blast/filter/";
#    $ENV{"BLAST_FILTER"} = "/usr/local/ncbi/blast/filter/";
    
    print join(", ",($base_dir,$method,$db,$query,$parameters,$CPUSERVER)),
        "\n" if $DEBUG;
    
    if ($refDb eq "ARRAY") {
    
    } # make .fa file of db
    if ($refQuery eq "ARRAY") {
    
    } # make .fa file of query sequence
    
    unless (4) { # ? just look for 4 files with the $db base name?
    
    } # prepare database for blast if not already done
    
    if ($ENV{HOST} && $ENV{HOST} =~ /^$CPUSERVER/) {
    
    } # no need to set up for and use rsh
    else {
        # make sure the names contain the full path
        if ($db !~ m{/\w+/.+}) { $db = $ENV{PWD} . "/" . $db }
        if ($query !~ m{/\w+/.+}) { $query = $ENV{PWD} . "/" . $query }
        # change this to the location of blast methods
        my $executable = $base_dir . $method;
        
        # this currently only works for tcsh I think, and obviouslt /etc/cshrc.aliases
        # is a BMERC convention. So, a way needs to be found to handle this part
        # either by finding a way to make rsh DTRT and run login stuff to set 
        # ENV, or I'll need to have rsh run a utility of my design that gets and
        # sets nec. stuff and runs blast - what a waste!
        my $envstring = "setenv WUBLASTMAT = $ENV{WUBLASTMAT};source /etc/cshrc.aliases";
        
        my $exestring = "/usr/bin/rsh -n $CPUSERVER \"$envstring;$executable " .
            "$db $query $parameters";
        print $exestring if $DEBUG;
        
        # if not debugging, supress blasts normal warning, whch includes
        # sillyness like not reporting some alignments due to the -B value
        if (! $DEBUG) {
            $result = `$exestring 2>/dev/null`;
        } # I don't think this will work off *nix - how else to do?
        else { $result = `$exestring` }
    } # use rsh to launch blast method
        
    return $result;
} # blast

# debug context sensitive warn/die
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



=head2 AA_HASH

Creates a hash table with amino acid lookups by codon. Includes all cases where
even an alternate na code (such as M for A or C) would return an unambiguous aa.
Also consistent with the complement method in this package, ie, lower cases in
some contexts, for ease of use with six_frame.

 C<%aa_hash = aa_hash;>

=cut
sub aa_hash {
my %aa_hash =
      ( TTT => "F",
    	TTC => "F",
	TTY => "F",
	TTr => "F",
	TTA => "L",
	TTG => "L",
	TTR => "L",
	TTy => "L",
	
	TCT => "S",
	TCC => "S",
	TCA => "S",
	TCG => "S",
	TCX => "S",
	TCN => "S",	
	TCM => "S",
	TCR => "S",
	TCW => "S",
	TCS => "S",
	TCY => "S",
	TCK => "S",
	TCV => "S",
	TCH => "S",
	TCD => "S",
	TCB => "S",
	TCm => "S",
	TCr => "S",
	TCw => "S",
	TCs => "S",
	TCy => "S",
	TCk => "S",
	TCv => "S",
	TCh => "S",
	TCd => "S",
	TCb => "S",
	
	TAT => "Y",
	TAC => "Y",
	TAY => "Y",
	TAr => "Y",
	TAA => ".",
	TAG => ".",
	TAR => ".",
	TAy => ".",

	TGT => "C",
	TGC => "C",
	TGY => "C",
	TGr => "C",
	TGA => ".",
	TGG => "W",

	CTA => "L",
	CTC => "L",
	CTG => "L",
	CTT => "L",
	CTX => "L",
	CTN => "L",
	CTM => "L",
	CTR => "L",
	CTW => "L",
	CTS => "L",
	CTY => "L",
	CTK => "L",
	CTV => "L",
	CTH => "L",
	CTD => "L",
	CTB => "L",
	CTm => "L",
	CTr => "L",
	CTw => "L",
	CTs => "L",
	CTy => "L",
	CTk => "L",
	CTv => "L",
	CTh => "L",
	CTd => "L",
	CTb => "L",

    	CCA => "P",
	CCC => "P",
	CCG => "P",
	CCT => "P",
	CCX => "P",
	CCN => "P",
	CCM => "P",
	CCR => "P",
	CCW => "P",
	CCS => "P",
	CCY => "P",
	CCK => "P",
	CCV => "P",
	CCH => "P",
	CCD => "P",
	CCB => "P",
	CCm => "P",
	CCr => "P",
	CCw => "P",
	CCs => "P",
	CCy => "P",
	CCk => "P",
	CCv => "P",
	CCh => "P",
	CCd => "P",
	CCb => "P",

	CAT => "H",
	CAC => "H",
	CAY => "H",
	CAr => "H",
	CAA => "Q",
	CAG => "Q",
	CAR => "Q",
	CAy => "Q",

	CGT => "R",
	CGC => "R",
	CGA => "R",
	CGG => "R",
	CGX => "R",
	CGN => "R",	
	CGM => "R",
	CGR => "R",
	CGW => "R",
	CGS => "R",
	CGY => "R",
	CGK => "R",
	CGV => "R",
	CGH => "R",
	CGD => "R",
	CGB => "R",
	CGm => "R",
	CGr => "R",
	CGw => "R",
	CGs => "R",
	CGy => "R",
	CGk => "R",
	CGv => "R",
	CGh => "R",
	CGd => "R",
	CGb => "R",

	ATT => "I",
	ATC => "I",
	ATA => "I",
	ATY => "I",
	ATr => "I",
	ATW => "I",
	ATM => "I",
	ATk => "I",
	ATH => "I",
	ATd => "I",
	ATG => "M",

	ACT => "T",
	ACC => "T",
	ACA => "T",
	ACG => "T",
	ACX => "T",
	ACN => "T",	
	ACM => "T",
	ACR => "T",
	ACW => "T",
	ACS => "T",
	ACY => "T",
	ACK => "T",
	ACV => "T",
	ACH => "T",
	ACD => "T",
	ACB => "T",
	ACm => "T",
	ACr => "T",
	ACw => "T",
	ACs => "T",
	ACy => "T",
	ACk => "T",
	ACv => "T",
	ACh => "T",
	ACd => "T",
	ACb => "T",

	AAT => "N",
	AAC => "N",
	AAY => "N",
	AAr => "N",
	AAA => "K",
	AAG => "K",
	AAR => "K",
	AAy => "K",

	AGT => "S",
	AGC => "S",
	AGY => "S",
	AGr => "S",
	AGA => "R",
	AGG => "R",
	AGR => "R",
	AGy => "R",

	GTT => "V",
	GTC => "V",
	GTA => "V",
	GTG => "V",
	GTX => "V",
	GTN => "V",	
	GTM => "V",
	GTR => "V",
	GTW => "V",
	GTS => "V",
	GTY => "V",
	GTK => "V",
	GTV => "V",
	GTH => "V",
	GTD => "V",
	GTB => "V",
	GTm => "V",
	GTr => "V",
	GTw => "V",
	GTs => "V",
	GTy => "V",
	GTk => "V",
	GTv => "V",
	GTh => "V",
	GTd => "V",
	GTb => "V",

	GCT => "A",
	GCC => "A",
	GCA => "A",
	GCG => "A",
	GCX => "A",
	GCN => "A",	
	GCM => "A",
	GCR => "A",
	GCW => "A",
	GCS => "A",
	GCY => "A",
	GCK => "A",
	GCV => "A",
	GCH => "A",
	GCD => "A",
	GCB => "A",
	GCm => "A",
	GCr => "A",
	GCw => "A",
	GCs => "A",
	GCy => "A",
	GCk => "A",
	GCv => "A",
	GCh => "A",
	GCd => "A",
	GCb => "A",

	GAT => "D",
	GAC => "D",
	GAY => "D",
	GAr => "D",
	GAA => "E",
	GAG => "E",
	GAR => "E",
	GAy => "E",
	
	GGT => "G",
	GGC => "G",
	GGA => "G",
	GGG => "G",
	GGX => "G",
	GGN => "G",	
	GGM => "G",
	GGR => "G",
	GGW => "G",
	GGS => "G",
	GGY => "G",
	GGK => "G",
	GGV => "G",
	GGH => "G",
	GGD => "G",
	GGB => "G",
	GGm => "G",
	GGr => "G",
	GGw => "G",
	GGs => "G",
	GGy => "G",
	GGk => "G",
	GGv => "G",
	GGh => "G",
	GGd => "G",
	GGb => "G",	
	);
return %aa_hash;
} # aa_hash

sub AUTOLOAD {
    my $program = $AUTOLOAD;
    $program =~ s/.*:://;
    print "AUTOLOAD being used to call $program\n" if $DEBUG >= 3;
    
    if ($program !~ /tbl/ && $program =~ s/cdna/tbl/) {
        no strict 'refs';
        print "Trying conversion as $program\n" if $DEBUG >= 2;
        return $program(@_);
    } # method didn't exist
} # AUTOLOAD

1;
__END__

=head1 EXPORT

None by default.

=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC -n CompBio

=item 0.20

Begin porting core function code for most initial methods from various
sources at BMERC. Write protocode and placeholders for what I want in
initial 'complete' version.

=item 0.44

Almost everything eneded up being rewritten practically from scratch.
Removed lingering locale assumptions and (hopefully) improving
interface (added use of %params mainly).
Added OOP useability so this module should now work either way.

=item 0.45

Modifications to Simple primarilly.

=item 0.451

Fixed faulty statement in MANIFEST.SKIP that excluded Makefile.PL from
distribution (thanks to Andreas Riechert for pointing this out to me almost
imediately).

=item 0.46

Added some documentation, mostly on options.
Added the ALTCODE option to six_frame.
Finished modifying tbl_to* converters to make a pass at handeling extra
data feilds from table format and changing *_to_tbl converters to keep
extra data besides id in extra fields in table format.

=item 0.461

Started working on bringing the utility programs included into using
CompBio::Simple, using dna_to_aa as first trial. Also some small changes
to docs.

=item 0.464

Major updates and alterations. Got most of the CompBio utils at least minimally
working. Added $RET_CODE usage to CompBio, setting old behavior as default, but
allowing Simple to provide as a parameter, so that return formating (such as
setting dna_to_aa to return in fasta) and file/STDOUT behavior is done during
process, reducing memory usage. Matching major alterations to Simple.

=item 0.466

NOTE .464 was acidently never fully distributed as such. Sorry!

Project has been added to SourceForge and the modules list on CPAN.

Worked out some more bugs, mostly cropping up from interactions with utils
such as tbl_to_fa, which interface through Simple. Mainly found places where
lazy programming broke using the files tied to an array, so had to use less of
$_. Things are shaking out nicely, and I think all the core bits are now in
place. Once the utils are done, I'll do a major code cleanup and comment run
and step up to 0.5 and push release. Finishing adding converters and fixing any
reported bugs (new converters shouldn't add any new bugs, except inside the new
methods) will benchmark going to 0.6.

=item 0.467

Minor fixes, some work on getting utils ready.

=item 0.468

Testing utils and started working on getting them to state. Improved fa_to_tbl
basic genbank | parsing without using CLEAN, but still needs some work. I despair
ever resolving this to my satisfaction, but at least now almost all the stuff not
used is kept in the extra field allowed in the latest table file spec. .47 line
will be working on getting the utils all working and writing tests for them. The
test script might end up as long as the module at this rate!

=item 0.469

More minor fixes to CompBio. Retruned to the version of Simple.pm from the .466
dist., as it was the last to complete all it's tests. I don't know how soon I'll
be able to resolve those problems. But the base version of the DB module is now
in, and overall I'm very pleased. It hasn't strayed far from the original BMERC
module as far as purpose of methods, but the code is a lot nicer. I've also
included the sql scripts for generating the base set of MySQL databases it's
designed to work with, and the .MRG files I use.

=back

=head1 TO DO - in no particular order

I like the basic design so far - except where large sequence sets are to be munged.
There this becomes a seriouse memory hog. Some faster & more efficient way
needs to be provided when dealing with files & large data sets.

complement and dna_to_aa need to be loop enabled, also alowing
&$RET_CODE(\@ret,$str); to be used. Why am I using _munge_array?!?

**Got it - needs work!** Get the full release of WashU Blast and make sure the
blast stuff works with it as well.

OK, rethink blast method(s) completely. Old and kludged and needs to die.

Add a method for handeling ncbi blast. CompBio::Simple should only have a blast
method that will DTRT and use the appropriate core methods. Need to revisit
only loading modules as needed, as this method should use the XML output
preferentially.

Stop using lowercase abiguity codes after complementing! (i.e. out of %AA)

If make install isn't putting right perl line in executibles, run perl config
myself?

Add DNA* and GCG format handelers

Add handler for the genbank report type format:
1   attgc   gtgct
11  gtgtg   gacaa
Which, though annoying, seems a certain candidate for recieving in cut and paste
operations, particularly through the web.

Better way to handle CPUSERVER. There must be some way to allow the user to
define (presumably during ./config?) how to submit jobs for more intensive
computational tools such as Blast and PIMA. Things like rsh, submitting to a
batch queue, etc should all be definable at install somehow.

write configure script to populate these globals as part of install
process.

modify tests to look for inclusion of Profile.pm or such and skip testing
where apropriate if the PIMA suite is not installed.

Find out if there is a way using OOP to allow a user to only need to include
this module and create new objects for the submodules, or even better just
DTRT. Can this be done just through inheritence? Should it be triggered by
requested export (like use CompBio qw(Simple DB)? Or through an AUTOLOAD
type interface, returning the correct object ($cbs = CompBio->Simple("new")
or some such)? I think that would be far more desirable than _having_ to
use a bunch of modules all in the CompBio namespace.

_error (all packages) needs to correctly report line where error occured.
Can this be done through caller or do I need to pass manually?

hmmm, how about a simple tm calculator for primers? I seem to be able to
find lots of them that work through web interface, or more complex primer
selection software, but none that plug and play to a perl program, and I'm
needing one. Not a primer selecter, just a _solid_ (more than just GC content),
straight foward Tm calculator.

=head1 COPYRIGHT

Developed at the BioMolecular Engineering Research Center at Boston
University under the NHLBI's Programs for Genomic Applications grant.

Copyright Sean Quinlan, Trustees of Boston University 2000-2002.

All rights reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Sean Quinlan, seanq@darwin.bu.edu

Please email me with any changes you make or suggestions for changes/additions.
Latest version is available through SourceForge
L<http://sourceforge.net/projects/compbio/>, or on the CPAN
L<http://www.cpan.org/authors/id/S/SE/SEANQ/>.
Thank you!

I would like to thank the staff at the BMERC for being my guinee pigs for
the past few years, Jim Freeman who got me into this and left me with flint
and tinder rather than the ember in a clay pot he started with, and my boss
and (tor)mentor Dr. Temple F. Smith! I would never have gotten this far or
learned this much without them.

=head1 SEE ALSO

perl(1), CompBio::Simple(1).

=cut
