use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'CompBio::DB',
    'VERSION_FROM'	=> 'DB.pm', # finds $VERSION
    'PREREQ_PM'		=> {DBI => 1.14,
						DBD::mysql => 0}, # e.g., Module::Name => 1.1
    'EXE_FILES'         => [qw(bin/dbaccess bin/get_def bin/get_seq bin/get_org
        bin/get_blast)],
    'AUTHOR'            => 'Sean P. Quinlan <seanq@darwin.bu.edu>',
    'ABSTRACT_FROM'     => 'DB.pm',
);
