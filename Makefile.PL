use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'          => 'CompBio',
    'AUTHOR'        => 'Sean P. Quinlan <seanq@darwin.bu.edu>',
    'ABSTRACT'      => 'CompBio module suite developed at the BioMolecular Engineering Research Center <http://bmerc-www.bu.edu> to support or bioinformatics work',
    'PREREQ_PM'     => {Test::More => 0.07},
    'VERSION_FROM'  => 'CompBio.pm',
    'EXE_FILES'     => [qw(bin/tbl_to_fa bin/tbl_to_ig bin/fa_to_tbl
        bin/ig_to_tbl bin/dna_to_aa)],
);
