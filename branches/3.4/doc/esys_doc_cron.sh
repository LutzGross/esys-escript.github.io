#!/bin/bash

# keep a record of where we are
HERE=`pwd`

# set stuff up
. /opt/modules/default/init/bash
export MODULEPATH=$MODULEPATH:/raid2/tools/modulefiles/:/data/raid2/toolspp4/modulefiles/gcc-3.3.6

# load the relevant modules
module load vtk-4.2.1-MangledMesa_gcc-3.3.5

# get into position
cd $HOME/raid2/esys_doc_sandbox/esys13/trunk

# set up the environment
export ESYS_ROOT=`pwd`
export LD_LIBRARY_PATH=$ESYS_ROOT/lib:/raid2/tools/boost/lib:$LD_LIBRARY_PATH
export PATH=/raid2/tools/python-2.3.4/bin:$PATH
export PYTHONPATH=$ESYS_ROOT:$PYTHONPATH
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=TRUE
export OMP_NESTED=FALSE

# update the cvs
svn update > /var/tmp/svn_msgs.txt 2>&1

# make esys
module load intel_cc.80.055 
#echo "Making esys clean" > /var/tmp/scons_msgs.txt
#scons -c >> /var/tmp/scons_msgs.txt 2>&1
echo "Making esys install" >> /var/tmp/scons_msgs.txt
scons >> /var/tmp/scons_msgs.txt 2>&1

# run doxygen
module load doxygen
cd $HOME/raid2/esys_doc_sandbox/esys13/trunk/doc/
doxygen doxygen/doxygen_esys > /var/tmp/doxygen_msgs.txt 2>&1
#doxygen doxygen_esys_python > /var/tmp/doxygen_python_msgs.txt 2>&1

# run epydoc
cd $HOME/raid2/esys_doc_sandbox/esys13/trunk/
module load epydoc
epydoc --html -o doc/obj/epydoc -n esys esys > /var/tmp/epydoc_msgs.txt 2>&1

# copy the docs to the webserver
cd $HOME/raid2/esys_doc_sandbox/esys13/trunk/doc/obj/
echo "Copying python docs to the webserver" > /var/tmp/scp_msgs.txt 2>&1
scp -i ~/.cron-ess-rsync-key -r epydoc/ shake200:/home/www/esys/ >> /var/tmp/scp_msgs.txt 2>&1
echo "Copying C++ docs to the webserver" >> /var/tmp/scp_msgs.txt 2>&1
scp -i ~/.cron-ess-rsync-key -r doxygen/ shake200:/home/www/esys/ >> /var/tmp/scp_msgs.txt 2>&1

# make an html file to link to the errors if any
echo "<html>" > /var/tmp/docbuildlog.html
echo "<head>" >> /var/tmp/docbuildlog.html
echo "<title>Documentation autobuild log files</title>" >> /var/tmp/docbuildlog.html
echo "</head>" >> /var/tmp/docbuildlog.html
echo "<body>" >> /var/tmp/docbuildlog.html
DATE=`date`
echo "<b>Last Update: $DATE</b>" >> /var/tmp/docbuildlog.html
echo "<p>" >> /var/tmp/docbuildlog.html
echo "<a href=\"svn_msgs.txt\">svn_msgs.txt</a><br>" >> /var/tmp/docbuildlog.html
echo "<a href=\"scons_msgs.txt\">scons_msgs.txt</a><br>" >> /var/tmp/docbuildlog.html
echo "<a href=\"doxygen_msgs.txt\">doxygen_msgs.txt</a><br>" >> /var/tmp/docbuildlog.html
echo "<a href=\"epydoc_msgs.txt\">epydoc_msgs.txt</a><br>" >> /var/tmp/docbuildlog.html
echo "<a href=\"scp_msgs.txt\">scp_msgs.txt</a><br>" >> /var/tmp/docbuildlog.html
echo "</p>" >> /var/tmp/docbuildlog.html
echo "</body>" >> /var/tmp/docbuildlog.html
echo "</html>" >> /var/tmp/docbuildlog.html

scp -i ~/.cron-ess-rsync-key /var/tmp/svn_msgs.txt shake200:/home/www/esys
scp -i ~/.cron-ess-rsync-key /var/tmp/scons_msgs.txt shake200:/home/www/esys
scp -i ~/.cron-ess-rsync-key /var/tmp/doxygen_msgs.txt shake200:/home/www/esys
scp -i ~/.cron-ess-rsync-key /var/tmp/epydoc_msgs.txt shake200:/home/www/esys
scp -i ~/.cron-ess-rsync-key /var/tmp/scp_msgs.txt shake200:/home/www/esys

scp -i ~/.cron-ess-rsync-key /var/tmp/docbuildlog.html shake200:/home/www/esys/

# get back to "here"
cd $HERE
