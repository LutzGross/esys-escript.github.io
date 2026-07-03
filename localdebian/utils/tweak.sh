#!/bin/bash

# set paths for each package and inform launcher it is operating with 
# tools in known locations

#if [ $# -ne 2 ]
#then
#    echo "Incorrect parameter count" >&2
#    exit 1
#fi


WORK3=$1
WORK3M=$2

# python-escript3
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORK3/bin/run-escript
sed -i -e "s%ESCRIPT_ROOT=/usr/lib/python-escript%ESCRIPT_ROOT=/usr/lib/python3-escript%" $WORK3/bin/run-escript
sed -i -e "s%PYTHON_CMD=python%PYTHON_CMD=python3%" $WORK3/bin/run-escript

# python-escript3-mpi
sed -i -e "s%ESCRIPT_ROOT=/usr/lib/python-escript%ESCRIPT_ROOT=/usr/lib/python3-escript-mpi%" $WORK3M/bin/run-escript-mpi
sed -i -e "s%STDLOCATION=0%STDLOCATION=1%" $WORK3M/bin/run-escript-mpi
sed -i -e "s%PYTHON_CMD=python%PYTHON_CMD=python3%" $WORK3M/bin/run-escript-mpi

# documentation
#sed -i -e 's%<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>%%' $WORK3/release/doc/sphinx_api/*.html
#sed -i -e 's%<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>%%' $WORK3/release/doc/sphinx_api/*.html
