use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'          => 'CompBio::Simple',
    'VERSION_FROM'  => 'Simple.pm', # finds $VERSION
    'AUTHOR'        => 'Sean P. Quinlan <seanq@darwin.bu.edu>',
    'ABSTRACT'      => 'CompBio::Simple.pm  - Simple OO interface to some basic methods useful in bioinformatics.',
);
#    'ABSTRACT_FROM'     => 'Simple.pm'
