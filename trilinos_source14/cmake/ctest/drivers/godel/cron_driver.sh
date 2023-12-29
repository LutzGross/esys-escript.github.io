#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on godel: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/home/trilinos/install/bin/eg
export TDD_PARALLEL_LEVEL=4
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

# Ensure the tests can find the Intel libraries at runtime.
export LD_LIBRARY_PATH=/opt/intel/cc/10.1.015/lib

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on godel: `date`"
echo
