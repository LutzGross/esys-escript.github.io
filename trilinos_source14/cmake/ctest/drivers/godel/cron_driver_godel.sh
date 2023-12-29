#!/bin/bash

CTEST_EXE=/home/trilinos/install/bin/ctest
EG_EXE=/home/trilinos/install/bin/eg
BASEDIR=/home/trilinos/dashboards/development
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/godel
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"


if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

echo "_DAYOFWEEK=$_DAYOFWEEK"

if [ "$_DAYOFWEEK" == "Saturday" ] ; then
  _RUN_REGULAR_TESTS=1
  _RUN_COVERAGE_TESTS=1
  _RUN_MEMCHECK_TESTS=0
elif [ "$_DAYOFWEEK" == "Sunday" ] ; then
  _RUN_REGULAR_TESTS=1 
  _RUN_COVERAGE_TESTS=0
  _RUN_MEMCHECK_TESTS=0
else
  _RUN_REGULAR_TESTS=1
  _RUN_COVERAGE_TESTS=0
  _RUN_MEMCHECK_TESTS=0
fi

echo "_RUN_REGULAR_TESTS=$_RUN_REGULAR_TESTS"
echo "_RUN_COVERAGE_TESTS=$_RUN_COVERAGE_TESTS"
echo "_RUN_MEMCHECK_TESTS=$_RUN_MEMCHECK_TESTS"

#exit


echo
echo "Starting nightly Trilinos testing on godel: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

cd $BASEDIR
if [ -d Trilinos ]; then
  echo Doing an update of existing directory
  cd Trilinos
  $EG_EXE pull
  cd ..
else
  echo Cloning the repository because none exists yets
  $EG_EXE clone $TRILINOS_REPOSITORY_LOCATION
fi


#
# Weekday tests
#

if [ "$_RUN_REGULAR_TESTS" == "1" ] ; then

echo
echo "Running weekday tests ..."
echo

echo
echo "Doing serial performance build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_performance_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_performance_godel.out

#echo
#echo "Doing mpi optimized build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_optimized_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_optimized_godel.out
#
#echo
#echo "Doing serial debug build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_serial_debug_godel.out
#
#echo
#echo "Doing mpi optimized shared library build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_optimized_shared_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_optimized_shared_godel.out

#echo
#echo "Doing serial optimized implicit instantiation build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_opt_impl_instant_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_serial_opt_impl_instant_godel.out

#echo
#echo "Doing mpi debug build: `date`"
#echo
#
#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_debug_godel.cmake -VV \
#  &> $BASEDIR/ctest_linux_nightly_mpi_debug_godel.out

echo
echo "Doing mpi optimized zoltan c-only build (32bit and 64bit global IDs): `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_opt_zoltan_c_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_opt_zoltan_c_godel.out

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_opt_zoltan_c64_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_opt_zoltan_c64_godel.out

echo
echo "Doing serial debug intel build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_icpc_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_icpc_godel.out

echo
echo "Doing serial release intel build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_release_icpc_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_release_icpc_godel.out

fi


#
# Coverage tests
#

if [ "$_RUN_COVERAGE_TESTS" == "1" ] ; then

echo
echo "Running coverage tests ..."
echo

echo
echo "Doing serial debug coverage build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_coverage_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_coverage_godel.out

echo
echo "Doing mpi debug coverage build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_debug_coverage_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_debug_coverage_godel.out

fi


#
# Memcheck tests
#

#memory testing has been disabled due to hanging
if [ "$_RUN_MEMCHECK_TESTS" == "1" ] ; then

echo
echo "Running memcheck tests ..."
echo

echo
echo "Doing serial debug memcheck build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_memcheck_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_memcheck_godel.out

echo
echo "Doing mpi debug memcheck build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_debug_memcheck_godel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_debug_memcheck_godel.out

echo
echo "Kill remaining 'memcheck' processes: `date`"
echo

killall -s 9 memcheck

fi

#forcing all files/directories to be group accessible
chgrp -R trilinos-dev *
chmod -R g+w *


echo
echo "Ending nightly Trilinos testing on godel: `date`"
echo


#/home/rabartl/mailmsg.py "Finished nightly Trilinos tests godel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
