#!/bin/sh

# Create the finley wrapper script bin/finleypython
# by substituting @@VAR@@ for an env variable $VAR

PYTHON_CMD=`type python | sed -e 's/.* //'`

sed							\
  -e "s%@@ESCRIPT_ROOT@@%$ESCRIPT_ROOT%"		\
  -e "s%@@LD_LIBRARY_PATH@@%$LD_LIBRARY_PATH%"		\
  -e "s%@@PYTHONPATH@@%$PYTHONPATH%"			\
  -e "s%@@PYTHON_CMD@@%$PYTHON_CMD%"			\
  -e "s%@@PATH@@%$PATH%"				\

