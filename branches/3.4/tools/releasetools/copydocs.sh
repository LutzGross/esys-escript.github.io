#!/bin/bash

# this script builds the release documentation (according to scons docs) and copies it to the webserver.

# keep a record of where we are
HERE=`pwd`

# get into position
cd $HOME/src/svn/esys13/trunk

# set up the environment for shake47
export LD_LIBRARY_PATH="/home/elspeth/src/svn/esys13/trunk/lib"
export PYTHONPATH="${ESYS_ROOT}"
export ESYS_ROOT="/home/elspeth/src/svn/esys13/trunk"
export OMP_SCHEDULE="dynamic"
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=TRUE
export OMP_NESTED=FALSE

# update the cvs
svn update > /var/tmp/svn_msgs.txt 2>&1

# make esys
echo "Making esys install" >> /tmp/scons_msgs.txt
scons >> /tmp/scons_msgs.txt 2>&1

# doxygen stuff will happen later

# run doxygen
#cd $HOME/src/svn/esys13/trunk/doc/
#doxygen doxygen/doxygen_esys > /tmp/doxygen_msgs.txt 2>&1

# run epydoc
cd $HOME/src/svn/esys13/trunk/
scons docs

# copy the docs to the webserver
cd $HOME/src/svn/esys13/trunk/release/doc
echo "Copying python docs to the webserver" > /tmp/scp_msgs.txt 2>&1
scp -r epydoc/* shake200:/home/www/esys/esys13/release/epydoc >> /tmp/scp_msgs.txt 2>&1
echo "Copying user guide (pdf/html) to the webserver" >> /tmp/scp_msgs.txt 2>&1
scp -r user/html/* shake200:/home/www/esys/esys13/release/user/html >> /tmp/scp_msgs.txt 2>&1
scp user/guide.pdf shake200:/home/www/esys/esys13/release/user >> /tmp/scp_msgs.txt 2>&1
echo "Copying example zip and tarball files" >> /tmp/scp_msgs.txt 2>&1
scp escript_examples.tar.gz shake200:/home/www/esys/esys13/release/escript_examples.tar.gz >> /tmp/scp_msgs.txt 2>&1
scp escript_examples.zip shake200:/home/www/esys/esys13/release/escript_examples.zip >> /tmp/scp_msgs.txt 2>&1

#echo "Copying C++ docs to the webserver" >> /tmp/scp_msgs.txt 2>&1
#scp -i ~/.cron-ess-rsync-key -r doxygen/ shake200:/home/www/esys/ >> /tmp/scp_msgs.txt 2>&1

# make an html file to link to the errors if any
echo "<html>" > /tmp/docbuildlogrelease.html
echo "<head>" >> /tmp/docbuildlogrelease.html
echo "<title>Documentation autobuild log files</title>" >> /tmp/docbuildlogrelease.html
echo "</head>" >> /tmp/docbuildlogrelease.html
echo "<body>" >> /tmp/docbuildlogrelease.html
DATE=`date`
echo "<b>Last Update: $DATE</b>" >> /tmp/docbuildlogrelease.html
echo "<p>" >> /tmp/docbuildlogrelease.html
echo "<a href=\"svn_msgs.txt\">svn_msgs.txt</a><br>" >> /tmp/docbuildlogrelease.html
echo "<a href=\"scons_msgs.txt\">scons_msgs.txt</a><br>" >> /tmp/docbuildlogrelease.html
echo "<a href=\"doxygen_msgs.txt\">doxygen_msgs.txt</a><br>" >> /tmp/docbuildlogrelease.html
echo "<a href=\"epydoc_msgs.txt\">epydoc_msgs.txt</a><br>" >> /tmp/docbuildlogrelease.html
echo "<a href=\"scp_msgs.txt\">scp_msgs.txt</a><br>" >> /tmp/docbuildlogrelease.html
echo "</p>" >> /tmp/docbuildlogrelease.html
echo "</body>" >> /tmp/docbuildlogrelease.html
echo "</html>" >> /tmp/docbuildlogrelease.html

scp /tmp/svn_msgs.txt shake200:/home/www/esys/esys13/release
scp /tmp/scons_msgs.txt shake200:/home/www/esys/esys13/release
#scp /tmp/doxygen_msgs.txt shake200:/home/www/esys
scp /tmp/epydoc_msgs.txt shake200:/home/www/esys/esys13/release
scp /tmp/scp_msgs.txt shake200:/home/www/esys/esys13/release
scp /tmp/docbuildlogrelease.html shake200:/home/www/esys/esys13/release
#scp -i ~/.cron-ess-rsync-key /tmp/docbuildlog.html shake200:/home/www/esys/

# get back to "here"
cd $HERE
