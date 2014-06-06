
#Needs the .deb package, the testing tarball and this script to work

cd /tmp/buildd
dpkg -i *.deb
aptitude install --assume-yes escript python-sympy
tar -zxf testfiles.tar.gz
cd testfiles
./itest.sh `pwd`/build '-t3'
