#!/bin/bash


# To be run from an esys13 root directory which has had a full build done in it.
#This will make a directory of files which can be shipped elsewhere to test an install
# The script now accounts for a build directory being somewhere other than ./build


if [ $# -lt 2 ]
then
   echo "Usage: $0 build_directory targetdirectory"
   exit 1
fi

if [ -f $2 ]
then
   echo "Target exists and is not a directory"
   exit 2
fi

if [ ! -d $1 ]
then
    echo "Build dir either does not exist or is not a directory"
    exit 2
fi

bdir=$1
targetdir=$2

if [ "$dest" == ".." ]
then
   # coz if you call this from inside a directory called src, you
   # wipe out your working copy
   echo "Using .. as a target is a bad idea. Suggest ../tests"
   exit 2
fi

if [ ! -f itest.sh ]
then
   echo "itest.sh not found. Have you run a build?"
   exit 3
fi

if [ ! -d $targetdir/build ]
then
   mkdir -p $targetdir/build
fi

cp itest.sh $targetdir
find . -maxdepth 1 -type d -not -name '*debian' -not -name '.' -not -name build -not -name esys -not -name bin -not -name lib -not -name '.?*' -print0 | xargs -0 -I'{}' cp -r '{}' $targetdir

cp -r $bdir/* $targetdir/build
cd $targetdir || exit 4
    
find build -name '*.o' -print0 | xargs -0 rm -f
find build -name '*.os' -print0 | xargs -0 rm -f
find build -name '*.so' -print0 | xargs -0 rm -f 
find build -name '*.a' -print0 | xargs -0 rm -f
find build -name '*.pyc' -print0 | xargs -0 rm -f
find build -name '*.passed' -print0 | xargs -0 rm -f
find build -name '*.skipped' -print0 | xargs -0 rm -f

find . -name 'src' -print0 | xargs -0 rm -rf
rm -rf scons
rm -rf doc/user doc/cookbook 
find doc -name '*.tex' -print0 | xargs -0 rm -f
