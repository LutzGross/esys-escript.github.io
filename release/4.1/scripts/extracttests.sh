#!/bin/bash


# To be run from an esys13 root directory which has had a full build done in it.
#This will make a directory of files which can be shipped elsewhere to test an install


if [ $# -lt 1 ]
then
   echo "Usage: $0 targetdirectory"
   exit 1
fi

if [ -f $1 ]
then
   echo "Target exists and is not a directory"
   exit 2
fi

if [ "$1" == ".." ]
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

if [ ! -d $1 ]
then
   mkdir $1
fi

targetdir=$1

cp -r * $targetdir
cd $targetdir || exit 4

rm -rf esys
rm -rf bin
rm -rf lib
rm -rf utest.sh
rm -rf 
find build -name '*.o' | xargs rm
find build -name '*.os' | xargs rm
find build -name '*.so' | xargs rm
find build -name '*.a' | xargs rm
find build -name '*.pyc' | xargs rm
find . -name 'src' | xargs rm -r
rm -r scons
rm -r doc/user doc/cookbook 
find doc -name '*.tex' | xargs rm
rm -rf debian
