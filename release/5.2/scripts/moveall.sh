#!/bin/bash

# Rewrites dependencies on all libraries in /lib to use /lib versions 
# rather than /build versions

for name in `find lib esys -name '*.dylib' -o -name '*.so'`
do
   echo $name
   scripts/libmover.sh $name `pwd`/lib
done
