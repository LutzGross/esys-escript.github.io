
DEST=../testfiles
mkdir $DEST
if [ $? != 0 ]
then
 echo "$DEST directory exists. Exiting"
 exit 1
fi

cp -r * $DEST
cd $DEST

rm -rf bin
rm -rf debian/
rm -rf lib
rm -rf packaging
rm -f config.log CREDITS.txt log README_LICENSE SConstruct svn_version utest.sh
find . -name '.s*' | xargs rm -rf
find . -name 'src' | xargs rm -rf
find . -name 'py_src' | xargs rm -rf
find . -name 'SConscript' | xargs rm -rf
find . -name '*.c' | xargs rm -rf
find . -name '*.h' | xargs rm -rf
find . -name '*.cpp' | xargs rm -rf
rm -rf tools site_scons scons scripts
rm -rf release
find build -type f | xargs rm
rm -rf esys
find . -name '__pycache__' | xargs rm -rf
find . -name '*.tex' | xargs rm -rf
cd doc
rm -rf epydoc cookbook doxygen install inversion manpage user *.sh *.cls
cd ..
find . -name '*.pyc' | xargs rm

cd ..
tar -czf testfiles.tar.gz testfiles
rm -rf testfiles




