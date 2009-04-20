#!/bin/sh

# Create a source distribution file for Escript/Finley

# Before running choose the release name and Subversion revision number

# This script creates a subdirectory ./escript-2.0.0beta and makes a tar file from it

name='escript-2.0.0beta'	# Must be a valid directory name
revision='2098'

svn export -r $revision https://shake200.esscc.uq.edu.au/svn/esys13/trunk $name || exit 1

cd $name || exit 1

# Remove some files we don't want to distribute

/bin/rm -rf autodocs.sh autotest-crontab autotest-scons README_TESTS

cat << _EOF_ > README

Escript is a python-based environment for implementing mathematical models, in
particular those based on coupled, non-linear, time-dependent partial
differential equations. It implements the finite element method. The code has
been parallelized efficiently with both MPI and OpenMP.

This is release $name of Escript/Finley based on Subversion
revision $revision.

`date '+%A, %B %d, %Y.'`

The User Guide is available in doc/guide.pdf.

A complete guide for compiling and installing is available at

	https://shake200.esscc.uq.edu.au/twiki/bin/view/ESSCC/EsysInstallationGuide

Example python scripts are available in doc/examples.

Complete documentation is available on our wiki at

	https://shake200.esscc.uq.edu.au/twiki/bin/view/ESSCC/WebHome

_EOF_

# Include the User Guide

(cd doc; wget http://shake200.esscc.uq.edu.au/esys/esys13/nightly/user/guide.pdf)

cd ..

tar cf $name.tar $name
gzip $name.tar

echo ''
echo "Remember to set a tag 'RELEASE_$name' for this release"
echo ''

