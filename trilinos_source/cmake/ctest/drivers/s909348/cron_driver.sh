#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on s909348: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/Users/bmpersc/bin/eg
export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Nightly

# Trilinos settings:
#

#export CTEST_TEST_TYPE=Experimental

#export CTEST_DO_SUBMIT=FALSE

#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on s909348: `date`"
echo
