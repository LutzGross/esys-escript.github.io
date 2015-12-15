
#Make the source tarball for debian release
#Run this from a clean checkout

SRCVERSION=`head -1 debian/changelog | tr -d '()' | tr -s '-' ' '| cut -d\  -f3`

svnversion | grep -q :
if [ $? == 0 ]
then
    echo "This does not appear to be a clean checkout."
    echo "Exiting"
    exit 1
fi
svnversion > svn_version

ls scons/*options.py > toexclude

tar -czf ../python-escript_$SRCVERSION.orig.tar.gz --exclude-vcs --exclude=debian --exclude=localdebian --exclude=toexclude --exclude-from toexclude *

rm toexclude
