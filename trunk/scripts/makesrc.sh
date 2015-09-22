
#Make the source tarball for debian release
#Run this from a clean checkout

SRCVERSION=`head -1 debian/changelog | cut -f2 -d- | cut -d\( -f2`

svnversion | grep -q :
if [ $? == 0 ]
then
    echo "This does not appear to be a clean checkout."
    echo "Exiting"
    exit 1
fi
svnversion > svn_version

tar -czf ../python-escript_$SRCVERSION.orig.tar.gz --exclude-vcs --exclude=debian --exclude=localdebian --exclude=scons/*options.py *

