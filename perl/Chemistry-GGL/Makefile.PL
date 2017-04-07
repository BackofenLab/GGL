use ExtUtils::MakeMaker;
use Config;

use strict;
use warnings;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Chemistry::GGL',
    AUTHOR            => q{Heinz Ekker <hekker@tbi.univie.ac.at>},
    VERSION_FROM      => 'lib/Chemistry/GGL.pm', # finds $VERSION
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PREREQ_PM         => { 
        'Chemistry::Mol' => 0,
        'Chemistry::File::SMILES' => 0,
    }, # e.g., Module::Name => 1.1
    CC                => 'g++',
    CCFLAGS            => '-Wno-deprecated',
    LIBS              => [ '-lstdc++ -L/usr/local/lib', ], 
    INC               => '-I. -I../../src -I../../src/vflib2-080701 -I../../src/sgd-1.3 -I../../src/nspdk', 
    MYEXTLIB          => '../../src/bin/toyChemUtil.o ../../src/sgm/*.o ../../src/ggl/*.o ../../src/vflib2-080701/*.o ../../src/nspdk/*.o ../../src/sgd-1.3/*.o',
    # OBJECT            => '$(O_FILES)', # link all the C files too
    dist              => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean             => { FILES => 'Chemistry-GGL-*' },
);
