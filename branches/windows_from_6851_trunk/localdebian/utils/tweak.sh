#!/bin/bash

# set paths for each package and inform launcher it is operating with 
# tools in known locations

if [ $# -ne 4 ]
then
    echo "Incorrect parameter count" >&2
    exit 1
fi


WORK=$1
WORKM=$2
WORK3=$3
WORK3M=$4

# python-escript
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORK/bin/run-escript2

# python-escript-mpi
sed -i -e "s%ESCRIPT_ROOT=/usr/lib/python-escript%ESCRIPT_ROOT=/usr/lib/python-escript-mpi%" $WORKM/bin/run-escript2-mpi
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORKM/bin/run-escript2-mpi

# python-escript3
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORK3/bin/run-escript3
sed -i -e "s%ESCRIPT_ROOT=/usr/lib/python-escript%ESCRIPT_ROOT=/usr/lib/python3-escript%" $WORK3/bin/run-escript3
sed -i -e "s%PYTHON_CMD=python%PYTHON_CMD=python3%" $WORK3/bin/run-escript3

# python-escript3-mpi
sed -i -e "s%ESCRIPT_ROOT=/usr/lib/python-escript%ESCRIPT_ROOT=/usr/lib/python3-escript-mpi%" $WORK3M/bin/run-escript3-mpi
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORK3M/bin/run-escript3-mpi
sed -i -e "s%PYTHON_CMD=python%PYTHON_CMD=python3%" $WORK3M/bin/run-escript3-mpi

