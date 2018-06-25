#!/bin/sh

# This script is for testing if the bundled source and test files can
# be unzipped to form an installation of escript. 

# This is the testing directoryi - set to an appropriate destination.
export SANDBOX="$HOME/sandbox"

# ESYS_ROOT should be set by the environment already.

cd $ESYS_ROOT
# make the tarballs via the 'release' target in scons
scons release

# move to the test area, and remove anything there.
cd $SANDBOX
rm -rf *

# copy over the tarballs for testing
cp $ESYS_ROOT/release/*.tar.gz .

# untar and run scons commands
tar -xzvf escript_src.tar.gz
scons
tar -xzvf escript_tests.tar.gz
scons all_tests

