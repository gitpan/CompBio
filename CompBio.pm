package CompBio;

require 5.005_62;
use strict;

require Exporter;

our @ISA = qw(Exporter);

# these are all expected to become unecisary
our $GENOME_HOME = "/seq/genome"; # base directory for genome databases
our $DBSERVER = "mysql"; # server type, as used for calling DBI::DBD driver
our $DBHOST = "sigler.bu.edu"; # location of server
# default machine to do heavy work like blast and PIMA on
our $CPUSERVER = "vavilov";

our %EXPORT_TAGS = ( 'all' => [ qw(check_type tbl_to_fa tbl_to_ig fa_to_tbl ig_to_tbl
    dna_to_protein complement six_frame aa_hash) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw($GENOME_HOME $DBSERVER $DBHOST $CPUSERVER);

our $VERSION = '0.44';
our $DEBUG = 0;

=head1 NAME

CompBio - Core library for some basic methods useful in computational biology/bioinformatics.

=head1 SYNOPSIS

use CompBio;
  
my $cbc = new->CompBio;

=head1 DESCRIPTION

The CompBio module set is _not_ intended to replace the bioperl project
(http://www.bioperl.org/). Although I welcome suggestions for improving
or adding to the methods available in these modules, particularly I would
love any help with things on the TO DO list, these modules are not intended
to provide the debth that the bioperl suite can provide. These modules
grew out of a set we have used at the BMERC(http://bmerc-www.bu.edu) and
worked with for years.

Originally developed at the BioMolecular Engineering Research Center
(http://bmerc-www.bu.edu), these modules and utilities grew out of a set we
have used and worked with for years. CompBio.pm is intended to take a number
of small, commonly used methods, and make them into a single package. Many
of the utils are just command line interfaces to these methods. To get the
most out of this set I highly recomend installing CompBio::DB.pm and either
importing our databases (ftp://mcclintock.bu.edu/BMERC/mysql/) or adapting
it to your local needs. Suggestions for improving portability are welcome!

The early versions of this module assumed installation on our local system.
Although I have tried to correct this in the current version, you may find
this package requires a litle twidling to get working. I'll try to leave
comments where I think it is most likely, but hopefully use of a relational
database and local setting changes in the globals will have
taken care of it. If not _please_ email me at seanq@darwin.bu.edu with the
details.

CompBio has a limited API. It expects it's input to
be in specific formats, as described in each methods docs, and it's output
is in a format that makes the most sense to me for that method. It does
no error checking by and large, so incorrect input could cause bizzare
behavior and/or a noisy death. If you want a flexible interface with lots
of error checking and deep levels of vebosity, use L<CompBio::Simple> - that's
its job.

Thanks!

Other modules available (or that will be available) in the CompBio set are:

The DB module will only be imediately useful if you import the databases as
used by us here at the BMERC(ftp://mcclintock.bu.edu/BMERC/mysql/) or develop
your own on the same basic design scheme. Otherwise I hope you find it useful
as a source of design ideas for rolling your own.

The Profile module was designed to work with our PIMA-II sofware and the
PIMA modules. The PIMA suite is available for license from Boston University
for a nominal fee, and free for academic use. For examples and more info
see http://bmerc-www.bu.edu/PIMA/. Unless you have or are interested in that
package, this module will have no functional value.

=head1 Methods

A couple of important general notes. All methods will describe the key required
arguments and will also indicate which arg is the optional %params hash. %params
can be passed as a hash or hash reference. Also, most arguments that handle sequences
or ids accept input as scalar, arrays (by reference) or hashes (by reference).
If I miss documenting it, try it just in case. Unless a parameter option is
available to signify otherwise, I assume 'doing the right thing' is to return the
same type of data construct as recieved.

=head2 new

Construct an object for invoking methods in CompBio.

=cut
sub new($%) {
    my ($proto,%parameters) = @_;
    my $class = ref($proto) || $proto;
    my $self = {};

    #handle params as nec. such as setting debug or changing env. variables
    $DEBUG = $parameters{'DEBUG'};
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
    print "Checking type for:\n",join("\n",@$seq[0..2]),"\n" if $DEBUG >= 3;
    
    if ($$seq[0] =~ /^.+?\t[CTUAGctuag]+$/) { return "CDNA" }
    elsif ($$seq[0] =~ /^.+?\t[A-Za-z!\.\*]+$/) { return "TBL" }
    elsif ($$seq[0] =~ /^>.+\n[A-Za-z!\.\*]+[\.\*]?$/m) { return "FA" }
    elsif ($$seq[0] =~ /^;\s*\n;\s*\n;\s*\n.+?\n[A-Za-z!\.\*]+1?S/) { return "IG" }
    elsif ($$seq[0] =~ /^[CTUAGMRWSYKVHDBXNctuagmrwsykvhdbxn]+\*?\n?$/) { return "RAW" }
    else { return "UNKNOWN" }
} # check_type

=head2 tbl_to_fa

Converts a sequence record in table (tab delimited, usually .tbl file extension)
format to fasta format.

C<$aref_faseqs = $cbc->tbl_to_fa(\@seqdat,%params);>

Each index in the @seqdat array must contain entire record (loci\tsequence) for
single sequence. Return is an array reference, still one sequence per index.

=cut
sub tbl_to_fa {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();

    foreach (@$aref_seqs) {
        chomp;
        my @fields = split(/\t/);
        my $str = ">$fields[0]\n";
        # generate fasta sequence lines at 80 char per line (including newline)
        my $tmpl = "a79" x ((length($fields[1])/79) + 1);
        $str .= join("\n",(unpack($tmpl,$fields[1])));
        $str =~ s/\n$//;
        push(@ret,$str);
    } # foreach submitted sequence
    
    return \@ret;
} # tbl_to_fa

=head2 tbl_to_ig

Converts a sequence in table (tab delimited) format to .ig format. Accepts
sequences in a referenced array, one record per index.

C<$aref_igseqs = $cbc->tbl_to_ig(\@tbl_seqs,%params);>

=cut
sub tbl_to_ig {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();
    
    foreach (@$aref_seqs) {
        chomp;
        my @fields = split(/\t/);
        my $str = ";\n;\n;\n$fields[0]\n"; # seperators, comment & id lines
        $fields[1] .= "1";
        my $tmpl = "a79" x ((length($fields[1])/79) + 1);
        $str .= join("\n",(unpack($tmpl,$fields[1])));
        $str =~ s/\n$//;
        push(@ret,$str);
    } # foreach submitted sequence

    return \@ret;
} # tbl_to_ig

=head2 fa_to_tbl

Accepts fasta format sequences in a referenced array, one complete sequence
record per index. This method returns the sequence(s) in table(.tbl) format
contained in a referenced array.

C<$aref_faseqs = $cbc->fa_to_tbl(\@fa_seq);>

=cut
sub fa_to_tbl {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();

    # traverse referenced array and convert
    foreach (@$aref_seqs) {
        chomp;
        my $tbl = "";
        
    	foreach (split(/[\n\r]+/)) {
    	    if (/^\s*>\s*([^\t]+)/) {
                my $sig = $1;
                # maybe a better keyword than CLEAN?
                # Also, may want to add to bad_characters in []
                if ($params{'CLEAN'} && $sig =~ /(\S+)/) {
                    $sig = $1;
                } # if user want and we can, get a better id
                elsif ($params{'REALCLEAN'} && $sig =~ /(\S{3,})[\|\!*\-]/) {
                    $sig = $1;
                } # elsif
    	        $tbl .= "$sig\t";
    	    } # if id line

    	    else {
	        s/[\*1\!\.]$//;
    	        $tbl .= $_;
    	    } # else
        } # foreach line in seqrec
        push(@ret,$tbl);
    } # foreach seqrec

    return \@ret;
} # fa_to_tbl

=head2 ig_to_tbl

Accepts ig format sequences in a referenced array, one complete sequence
record per index. This method returns the sequence(s) in table(.tbl) format
contained in a referenced array.

C<$aref_igseqs = $cbc->ig_to_tbl(\@fa_seq);>

=cut
sub ig_to_tbl {
    my $self = shift;
    my $aref_seqs = (ref($self) && ref($self) ne "ARRAY") ? shift : $self;
    _help() if (! ref($aref_seqs) && $aref_seqs eq "help");
    _error("Not given an array reference!\n",1) unless ref($aref_seqs) eq "ARRAY";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();

    foreach (@$aref_seqs) { # traverse referenced array and convert
    	chomp;
        my $tbl = "";
        
    	foreach (split(/[\n\r]+/)) {
            next if /^;/;
            if (! $tbl) {
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
        push(@ret,$tbl);
    } # foreach seqrec

    return \@ret;
} # ig_to_tbl

=head2 dna_to_protein

Converts a dna sequence, containing no whitespace, submited as a scalar reference,
to the amino acid residues as coded by the standard 'universal' genetic code.
Return is a reference to a scalar. dna sequence may contain standard special
characters (i.e. R,S,B,N, ect.). Default behavior is to trim a final stop, if
present, and to substitute an M for an I or L in the first position - this is
usually correct when translating whole sequences from a coding DNA sequence.
A hash containing optional parameters may be passed as the second argument.

Options allowed are:

C: Set to a true value to indicate dna should be converted to it's compliment
before translation.

ALTCODE: A reference to a hash containing alternate coding keys where the value 
is the new aa to code for. Stop codons are represented by ".".

SEQFIX: Set to true to alter first position, making V or L an M, and
removing stop in last position.

C<$aa = dna_to_protein(\$dna_seq,%params);>

=cut
sub dna_to_protein {
    my $self = shift;
    my $sref_seq = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    _help() if (! ref($sref_seq) && $sref_seq eq "help");
    _error("Not given a scalar reference!\n",1) unless ref($sref_seq) eq "SCALAR";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();

    my %AA = aa_hash();
    $$sref_seq = uc($$sref_seq);
    $$sref_seq =~ tr/U/T/;
    print $$sref_seq,"\n" if $DEBUG > 1;
    
    if($params{'ALTCODE'}) {
    	while ((my $codon,my $aar) = each %{$params{'ALTCODE'}}) {
    	    if ($codon =~ /^(\w\w)\[([A-Z]+)\]$/) {
                foreach (split("",$2)) { $AA{"$1$_"} = $aar }
            } # multiple options in third position
            else { $AA{$codon} = $aar }
    	} # for each alternate codon
    } # if

    if ($params{'C'}) { # convert all characters to thier compliment
    	$sref_seq = complement($sref_seq);
    } # if
    
    my $ret = "";
    # we allow any non-space character in match in case someone wants odd alternate coding
    while ($$sref_seq =~ /(\S{3})/g) {
        $ret .= $AA{$1} || "X";
    } # while translating sequence
    
    if ($params{'SEQFIX'}) {
        $ret =~ s/^[VL]/M/;
        $ret =~ s/\.$//;
    } # clean up ends for use as aa seq
    
    return \$ret;
} # dna_to_protein

=head2 complement

Converts dna to it's complimentary strand. DNA sequence is submitted as scalar
by reference. There is no return as sequence is modified. To maintain original
sequence, send a reference to a copy.

C<compliment(\$dna);>

=cut
sub complement {
    my $self = shift;
    my $sref_seq = (ref($self) && ref($self) ne "SCALAR") ? shift : $self;
    _help() if (! ref($sref_seq) && $sref_seq eq "help");
    _error("Not given a scalar reference!\n",1) unless ref($sref_seq) eq "SCALAR";
    
    my %params = @_ ? @_ : ();
    local $DEBUG = exists $params{DEBUG} ? $params{DEBUG} : $DEBUG;
    my @ret = ();

    $$sref_seq =~ tr/[ACTUGMRYKVHDBkyrmbdhv]/[TGAACkyrmbdhvMRYKVHDB]/;
    $$sref_seq = reverse($$sref_seq);

    return  $sref_seq;
} # complement

=head2 six_frame

Converts a submitted dna sequence into all 6 frame translations to aa. Arguments are the file
containing the raw dna sequence (in .raw format! no whitesace), the id to prefix the output,
the min length of amino acid sequences to be recorded, and the file to output. Note that the output
file will be truncated. Output id's have strand, frame, start and stop positions encoded.

C<$result = six_frame($raw_file,$id,$seq_len,$out_file);>
	
Note: six_frame returns a result of 0 on success, or an error message on failure.

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
    my @ret = ();

    $params{'ID'} ||= "SixFrame";
    $params{'SEQLEN'} ||= 10;
    $params{'OUTFILE'} ||= "";
    my %AA = aa_hash();
    my $fh_out = "";

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
    $cseq1 = $AA{${complement(\$codon)}} || "X";
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
    
    return $ret || 0;
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
    } # if seq from complimentary strand
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

A simple interface to the Washington University distribution of blast. Currently
has only been tested with the 2.0a8MP release.


 Takes series of arguments qw($method $database $queryfile 
"parameter list","noise",$cpu_server). Method can be any allowed blast method. Database is the blastable fasta format 
database to be searched. queryfile is the file containg the fasta sequence you wish to use in the 
blast comparison. "parameter list" is a quote inclosed set of parameter arguments to hand blast.
Optional final arguments are the noise or debig level to use (default is quiet), "quiet" tells blast
to operate silently, anything else prints some debugging messages and allows any errors from blast
to be printed. Final argument is for the 'cpu server' where blast is to execute and must follow
a noise level argument.

	Only blastp and blastx  have default parameters at moment:
		blastp  B=10  e=1e-3  filter=seg+xnu
		blastx  B=25  e=1e-3  filter=seg+xnu
	If any parameters are supplied method will not use any default values!

C<$blast_out = blast('blastp',$fa_file,$queryfile,"-B=0 -matrix=BLOSUM62","quiet"); >
	
=cut
sub wu_blast {
    my ($method,$database,$query_file,$parameters,$noise_level,$cpu) = @_;
    if ($method eq "help") { &help }
    $noise_level = "quiet" unless $noise_level;

    print "$method , $database , $query_file , $parameters , $cpu\n" if $noise_level ne "quiet";
    if ($database !~ m{/\w+/.+}) { $database = $ENV{PWD} . "/" . $database }
    if ($query_file !~ m{/\w+/.+}) { $query_file = $ENV{PWD} . "/" . $query_file }

    my $envstring = "";
    my $cpuserver = $cpu || "vavilov";
    my $executable = "";
    
    if (-e "/usr/local/bin/$method") { $executable = "/usr/local/bin/$method" }
    elsif ($ENV{HOST} && $ENV{HOST} eq $cpuserver) {
    	$executable = "/seq/genome/blastserver/$ENV{OSTYPE}/" . $method;
    } # if
    else {
    	my $exe = qq(/usr/bin/rsh -n $cpuserver "env | grep OSTYPE");
	print "$exe\n" if $noise_level ne "quiet";
    	chomp(my $temp = `$exe`);
	$temp =~ s/OSTYPE=//;
	unless ($temp eq "osf1" || $temp eq "solaris") {
	    die "Can't determine where blast executable is for $temp.";
	} # if
	$executable = "/seq/genome/blastserver/$temp/" . $method;
    } # else

    $ENV{"BLASTFILTER"} = "/usr/local/ncbi/blast/filter/";
    $ENV{"BLAST_FILTER"} = "/usr/local/ncbi/blast/filter/";

    if($method =~ /blast/) {
    	# these are not formated properly for a simple rsh sesion
	$envstring = "setenv BLASTFILTER /usr/local/ncbi/blast/filter/;source /etc/cshrc.aliases";
    } # if
    
    if (! defined $parameters) {
    	open (PARAMS,"$GENOME_HOME/blastserver/presets.tbl") or die "Can't open PARAMS for blast: $!\n";
    
    	while (defined(<PARAMS>)) {
	    if (/$method\t(\w+)\t(.+)/) { $parameters .= "-$1=$2 " }
    	} # while
    	close PARAMS or die "can't close PARAMS in blast: $!\n";
    } # if
    
    my $exestring = qq(/usr/bin/rsh -n $cpuserver "$envstring;$executable $database $query_file $parameters");
    if ($ENV{HOST} && $ENV{HOST} eq $cpuserver) { $exestring = "$executable $database $query_file $parameters" }
    my $result = "";
    if ($noise_level eq "quiet") { $result = `$exestring 2>/dev/null` }
    else { $result = `$exestring` }
        
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
even an alternate na code (such as M for A or C) would return an unambiguos aa.
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

1;
__END__

=head1 EXPORT

None by default.

=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC -n CompBio

=item 0.44

Copy over functions from original BMERC::bio (ver 0.74), making improvements to code,
mostly by removing lingering locale assumptions and (hopefully) improving
interface, and converting to OOP.

=back

=head1 TO DO - in no particular order

I like the basic design so far - except where large sequence sets are to be munged.
There this becomes a seriouse memory hog. Some faster & more efficient way
needs to be provided when dealing with files & large data sets.

**Got it - needs work!** Get the full release of WashU Blast and make sure the
blast stuff works with it as well.

OK, rethink blast method(s) completely. Old and kludged and needs to die.

Add a method for handeling ncbi blast. CompBio::Simple should only have a blast
method that will DTRT and use the appropriate core methods.

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

=head1 COPYRIGHT

Copyright Sean Quinlan, Trustees of Boston University 2000-2001.

All rights reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Sean Quinlan, seanq@darwin.bu.edu

Please email me with any changes you make or suggestions for changes/additions.
Latest version is available under ftp://mcclintock.bu.edu/BMERC/perl/. Thank you!

=head1 SEE ALSO

perl(1), CompBio::Simple, CompBio::DB.

=cut
