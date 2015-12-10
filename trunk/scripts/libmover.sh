

#This script will take a library on macos (.so or .dylib)
#and will replace the relative dependencies with absolute
#deps in targetlibpath
#
# eg:
# for name in `find lib -name '*.so'`;do echo $name;./libmover.sh $name `pwd`/lib;  done


if [ $# -ne 2 ]
then
    echo "Usage: $0 library_file_to_update targetlibpath"
    exit 1
fi

libname=$1
tlibpath=$2

lines=`otool -L $libname | tail +3 | cut -f1 -d\  | cut -f2 | grep build/darwin`
for name in $lines
do
   install_name_tool -change $name $tlibpath/`basename $name` $libname
done 


