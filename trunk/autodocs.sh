#!/bin/bash

# Build all escript documentation on shake71 and put it on the web server
# in the proper place so the ESSCC twiki can find it. The build command is
#	scons docs

# Run this nightly on shake71 via cron, such as at 1:00 AM every night, with
# 0      1    *    *     0-6      /home/Work/Documentation_Escript/autodocs.sh >& /home/Work/Documentation_Escript/log

# This will have an error occasionally because svn checkout fails. 
# The problem may be that the certificate has expired.  If that's the
# case try an svn checkout (or svn list) by hand and accept the
# certificate.  Then the problem will be resolved.

# Requires:
#	svn checkout from shake200
#	scp to shake200 to upload new documentation to web site
#	epydoc
#	doxygen
#	latex
#	latex2html

DIR="/home/Work/Documentation_Escript"

START=`date '+%Y/%m/%d %H:%M'`
RunDate=`date '+%Y_%m_%d'`

scons='/home/Work/scons-0.96.92/bin/scons'

finish () {
  # state will be 'FAILURE' or 'SUCCESS'
  state="$1"
  date
  # Clean up the sandbox
  cd $DIR
  ### /bin/rm -rf sandbox
  END=`date '+%Y/%m/%d %H:%M'`
  cat << END_MSG | mail -s "ESYS_TESTS docs $RunDate $state" k.steube@uq.edu.au
$2.
The tests ran from $START to $END
See the log file /home/Work/Documentation_Escript/log for info
This mail was sent by $0
running via cron as $USER on `hostname`.
END_MSG
  if [ "x$state" = "xFAILURE" ]; then
    exit 1
  fi
  exit 0
}

umask 022

cd $DIR			|| finish FAILURE "Could not cd to $DIR"

/bin/rm -rf sandbox
mkdir sandbox		|| finish FAILURE "Could not mkdir sandbox"
cd sandbox		|| finish FAILURE "Could not cd to sandbox"

echo "Checking out esys13/trunk"
svn checkout https://shake200.esscc.uq.edu.au/svn/esys13/trunk || finish FAILURE "Could not checkout esys13/trunk"

export PATH="/home/Work/latex2html-2002-2-1/bin:$PATH"
export LD_LIBRARY_PATH="$DIR/sandbox/trunk/lib:/home/Work/VTK-4.4.2/lib"
export PYTHONPATH="$DIR/sandbox/trunk"

# Generate documentation
echo "Generating documentation"

cd trunk							|| finish FAILURE "Could not cd to trunk"
mkdir release release/doc					|| finish FAILURE "Could not create release directory"
$scons dodebug=yes useMPI=no docs				|| finish FAILURE "Could not run scons docs"
scp -r release/doc/* shake200:/home/www/esys/esys13/nightly	|| finish FAILURE "Could not copy documentation to nightly area"

echo "Cleaning up"

finish SUCCESS "Successfully ran 'scons docs' on `hostname`"

