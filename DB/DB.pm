package CompBio::DB;

require 5.005_62;
use strict;
use DBI;
use Carp qw(cluck confess croak carp);
use CompBio;
use diagnostics;
BEGIN {
    disable diagnostics;
} # BEGIN - turn off diagnostics unless called for

require Exporter;
our @ISA = qw(Exporter);
use vars qw($AUTOLOAD);

our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(db_help);

our $VERSION = '0.2';
our $DEBUG = 0;

my $DB_HOST     = 'localhost';
my $DB_SERVER   = 'mysql';
my $DB_DATABASE = 'CompBio';
my $DB_USER     = 'nobody';
my $DB_PASSWORD = '';

if (-s $CONF_FILE) {
    open(CONF,$CONF_FILE) or die "Can't read configuration file $CONF_FILE: $!";
    while (<CONF>) {
        if (/DB_DATABASE\s*=\s*(\S+)/)  { $DB_DATABASE = $1 }
        elsif (/DB_SERVER\s*=\s*(\S+)/) { $DB_SERVER   = $1 }
        elsif (/DB_HOST\s*=\s*(\S+)/)   { $DB_HOST     = $1 }
        elsif (/DB_USER\s*=\s*(\S+)/)   { $DB_USER     = $1 }
        elsif (/DB_HOST\s*=\s*(\S+)/)   { $DB_PASSWORD = $1 }
    } # read
} # read config

$DB_HOST     = $ENV{COMPBIO_DB_HOST}     if $ENV{COMPBIO_DB_HOST};
$DB_SERVER   = $ENV{COMPBIO_DB_SERVER}   if $ENV{COMPBIO_DB_SERVER};
$DB_DATABASE = $ENV{COMPBIO_DB_DATABASE} if $ENV{COMPBIO_DB_DATABASE};
$DB_USER     = $ENV{COMPBIO_DB_USER}     if $ENV{COMPBIO_DB_USER};
$DB_PASSWORD = $ENV{COMPBIO_DB_PASSWORD} if $ENV{COMPBIO_DB_PASSWORD};

=head1 NAME

CompBio::DB - Methods for accessing data stored acording to the CompBio base schema.

=head1 SYNOPSIS

  use CompBio::DB;
  my $cbdb = CompBio::DB->new({host => "foo.bar.edu");

  my $AR_defs = $cbdb->get_annotation([keys %seqs]);

Example code for parsing the return for a simple print:

 my $AR_result = $cbdb->get_aa_seq(\@id_list);
 foreach my $AR_row (@$AR_result) {
    next unless @$AR_row;
    print join("\t",@$AR_row) , "\n";
 } # return list

=head1 DESCRIPTION

This module and the related database schemas where developed to be used with a
MySQL server. Although I would certainly prefer it to be portable, that was not
a priority at this time. However, I expect the manual alterations that may be
necisary to use this module should be very resonable for someone sufficiently
familiar with the new target database server, and I would be happy to collaborate
with anyone who wants to work on this!

Most of the methods in this database fetch a type of data from a given id or
list of ids. Unless otherwise stated for a specific method, only two arguments
are ever used; the id or an array reference to a list of ids, and
a hash reference containing any optional arguments. The request will be made
against the default database defined when the CompBio::DB object was created,
but most methods accept a "database" option to use a different database on the
same server.

All methods that return query results return a reference to a 2D array. See the
documentaion for specific methods for the order of returned fields.

=head2 EXPORT

db_help

=head2 db_help

Call to run perldoc on CompBio::DB. Uses exec to launch perldoc, so should be last
ditch effort in place of die, perhaps when caller can't negotiate proper usage and
will fail regardless.

=cut
sub db_help {
    exec 'perldoc CompBio::DB';
} # db_help

=head2 new

Creates the CompBio::DB object. This requires succesfully connecting to the
default database. There are no required arguments. Below
are the optional parameters. Defaults may be changed by setting them in
environment variables, such as COMPBIO_DBUSER, or by altering the defaults in
CompBio.conf.

B<OPTIONS:>

=over 2

=item *

host

Name or IP of the host on which the database server is located.
Defaults to localhost. Environment variable is COMPBIO_DB_HOST.

=item *

server

Type of database server to connect to (determines DBI driver).
Defaults to mysql. Environment variable is COMPBIO_DB_SERVER.

=item *

database

Name of the database to have handle connected to.
Defaults to Key. Environment variable is COMPBIO_DB_DATABASE.

=item *

user

Username to use in connection.
Defaults to nobody. Environment variable is COMPBIO_DB_USER.

=item *

password

Password to use in connection.
Defaults to empty string. Environment variable is COMPBIO_DB_PASSWORD.

=back

=cut
sub new {
    my $proto = shift;
    my $params = shift;
    my $class = ref($proto) || $proto;
    my $self = {};

    my %params = ();
    if (ref $params eq "HASH") { %params = %{$params} }
    elsif (@_) { %params = ($params,@_) }
    
    #handle params as nec. such as setting debug or changing env. variables
    $DEBUG = $params{DEBUG} || 0;
    if    ($DEBUG >= 3) { enable diagnostics }
    elsif ($DEBUG)      { $^W = 1 }
    
    $self->{'_created'} = 1;
    print "Created = ",$self->{'_created'},"\n" if $DEBUG >= 3;
    
    my $db       = $params{database} || $DB_DATABASE;
    my $user     = $params{user}     || $DB_USER;
    my $host     = $params{host}     || $DB_HOST;
    my $password = $params{password} || $DB_PASSWORD;
    my $server   = $params{server}   || $DB_SERVER;
    
    my $dsn = "DBI:$server:$db:$host";
    
    my $attemp_count = 1;
    my $dbh = "";
    
    CONNECT: {
    $dbh = DBI->connect($dsn,$user,$password);
    _error($DBI::errstr) if $DBI::err && $DEBUG >= 2;
    unless ($dbh) {
    	_error("Have no connection to DB, $dsn; retrying in 3") if $DEBUG;
	sleep(3);
	++$attemp_count;
	redo CONNECT unless $attemp_count > 10;
    } # no connection
    } # CONNECT control block
    
    $self->{dbh} = $dbh;
    return bless ($self,$class) if $dbh;
    _error("Failed to store database connection",1);
} # new

##
## Basic data fetch statements
##

=head2 get_aa_seq

Returns the peptide sequence for the supplied id(s).

B<OPTIONS:>

=over 2

=item *

database <database>

Specify a specific database to search. Default is the merged table.

=item *

ISKEY (0/1)

Signals that the id(s) supplied are already unique database keys. This is
preffered when speed is desired.

=item *

include_suffix (0/1)

Combines the suffix with the id supplied in the return, such as a genbank
accesion number's version. This argument is ignored when used with ISKEY.
Please note that if the id has a suffix which is not used in matching with
the require_suffix option, and the id is repeated in the database, you will
get only one of the possible results returned with no predictable preference.

=item *

require_suffix <suffix>

Requires the suplied suffix to match as well as the id. The value may be NULL
for some ids as appropriate. This option is ignored when option ISKEY is true.

=back

=cut
sub gen_get_aa_seq {
    my %params = @_;
    
	my $statement = '';
	my $db = $params{database} ? "$params{database}." : '';
    if ($params{ISKEY}) {
		$statement = "SELECT ID, Sequence FROM $db"
			. "Peptide WHERE ID = ?";
	} # ids are numeric atabase keys
	
	else {
		my $suffix = $params{include_suffix} ? ', Suffix' : '';
		$statement = "SELECT CONCAT(Name$suffix), Sequence FROM $db"
			. "Peptide, ID WHERE Name = ? AND ID = DBKey";
	} # ids are text
	
	return $statement;
} # gen_get_aa_seq

=head2 gen_get_annotation


B<OPTIONS:>

=over 2

=item *

=back

=cut
sub gen_get_annotation {
    my %params = @_;
    
	my $statement = '';
	my $db = $params{database} ? "$params{database}." : '';
    if ($params{ISKEY}) {
		$statement = "SELECT ID, Definition, Supplemental, Submitted FROM $db"
			. "Annotation WHERE ID = ?";
	} # ids are numeric atabase keys
	
	else {
		my $suffix = $params{include_suffix} ? ', Suffix' : '';
		$statement = "SELECT CONCAT(Name$suffix), Definition, Supplemental, Submitted FROM $db"
			. "Annotation, ID WHERE Name = ? AND ID = DBKey";
	} # ids are text
	
	return $statement;
} # gen_get_annotation

=head2 gen_get_cds


B<OPTIONS:>

=over 2

=item *

=back

=cut
sub gen_get_cds {
    my %params = @_;
    
	my $statement = '';
	my $db = $params{database} ? "$params{database}." : '';
    if ($params{ISKEY}) {
		$statement = "SELECT ID,Start,End,Strand,TotalSegments,ThisSegment,"
			. "GenomicElement FROM $db"
			. "FeatureCoding WHERE ID = ?";
	} # ids are numeric atabase keys
	
	else {
		my $suffix = $params{include_suffix} ? ', Suffix' : '';
		$statement = "SELECT CONCAT(Name$suffix),Start,End,Strand,"
			. "TotalSegments,ThisSegment FROM $db"
			. "FeatureCoding, ID WHERE Name = ? AND ID = DBKey";
	} # ids are text
	
	return $statement;
} # gen_get_cds

sub AUTOLOAD {
    my $method = $AUTOLOAD;
    $method =~ s/.*:://;
    print "AUTOLOAD being used to call $method\n" if $DEBUG >= 3;
    
	if ($method =~ /get_/) {
    	my $self = shift;
    	my $ids = shift;
    	my $params = shift;
    
    	my %params = ();
    	if (ref $params eq "HASH") { %params = %{$params} }
    	elsif (@_) { %params = ($params,@_) }
    	
		my $statement = '';
		{ no strict "refs";
		my $gen_get = "gen_" . $method;
		$statement = &$gen_get(%params);
		} # strict refs
		
		if (ref $ids) {
			_error("Not an ARRAY reference",1) unless ref $ids eq "ARRAY";
			return _fetch_multiple_id($self,$ids,$statement);
		} # multiple ids
	
		else { return _fetch_single_id($self,$ids,$statement) }
	} # get method
	
} # AUTOLOAD

##
## Internal functions
##

sub _fetch_single_id {
	my $self = shift;
	my $id = shift;
	my $query = shift;
	
	# query should have been prepared the same as for multiple ids
	# so we'll replace the placeholder here
	#$query =~ s/ \? / $id /;
	
	my $AR_ret = $$self{dbh}->selectall_arrayref($query,{},$id);
	_error("Problem executing $query: $DBI::errstr") if $DBI::err;
	return $AR_ret;
} # _fetch_single

sub _fetch_multiple_id {
	my $self = shift;
	my $AR_ids = shift;
	my $query = shift;
	
	my $sth = $$self{dbh}->prepare($query);
	_error("Problem preparing $query: $DBI::errstr") if $DBI::err;
	
	my $AR_ret = [];
	foreach my $sig (@$AR_ids) {
		$sth->execute($sig);
		_error("Problem executing $query: $DBI::errstr") if $DBI::err;
		my $AR_row = $sth->fetchall_arrayref;
		_error("Problem fetching from $query: $DBI::errstr") if $DBI::err;
		# since multiple rows could be returned on certain queries, such as
		# for FeatureCoding, we need to push all the row references returned
		# for each id onto @$AR_ret.
		push(@$AR_ret,@$AR_row[0 .. $#{$AR_row}]);
	} # foreach id
	$sth->finish;
	return $AR_ret;
} # _fetch_multiple

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

1;
__END__

# stub sub

=head2 method

stub method pod

B<OPTIONS:>

=over 2

=item *

=back

=cut
sub gen_ {
    my %params = @_;
    
	my $statement = '';
	my $db = $params{database} ? "$params{database}." : '';
    if ($params{ISKEY}) {
		$statement = "SELECT  FROM $db"
			. " WHERE ID = ?";
	} # ids are numeric atabase keys
	
	else {
		my $suffix = $params{include_suffix} ? ', Suffix' : '';
		$statement = "SELECT CONCAT(Name$suffix) FROM $db"
			. ", ID WHERE Name = ? AND ID = DBKey";
	} # ids are text
	
	return $statement;
} # 

=head1 HISTORY

=over 8

=item 0.01

Original version; created by h2xs 1.20 with options

  -AXC -n CompBio::DB

=item 0.2 Sun Jun 16 23:08:45 EDT 2002

Stable base code. Two internal functions for getting results from database
for a given query. Could be one but avoids looping for single ids. Handles
multiple possible row returns for single ids. Three working get methods passed
basic tests. Most of the parts of the methods that were always the same are
now in AUTOLOAD, which calls an internal method to generate the appropriate
SQL statement.

=back

=head1 TO DO

All basic optional parameters that effect building the should be defined now,
at least as much as possible, so that the stub for new methods can be fairly
complete and stable.

Fix case discrepencies between CompBio ID table and actual database names.

How to handle using suffixes, in search and returns.

A seperate README is needed on creating the databases!

=head1 COPYRIGHT

Developed at the BioMolecular Engineering Research Center at Boston
University under the NHLBI's Programs for Genomic Applications grant.

Copyright Sean Quinlan, Trustees of Boston University 2000-2002.

All rights reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Sean Quinlan, seanq@molbio.mgh.harvard.edu

Please email me with any changes you make or suggestions for changes/additions.
Latest versions are available through SourceForge
L<http://sourceforge.net/projects/compbio/>, or you can get the latest release
version on the CPAN
L<http://www.cpan.org/authors/id/S/SE/SEANQ/>.
Thank you!

=head1 SEE ALSO

perl(1), CompBio.

=cut
