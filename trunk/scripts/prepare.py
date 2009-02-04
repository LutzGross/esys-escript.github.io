
import shutil, os, datetime, sys


SVNURL="https://shake200.esscc.uq.edu.au/svn/esys13/trunk"
NUMJs=4
TOPDIR=str(datetime.date.today())
ERRMAIL="j.fenwick1@uq.edu.au"

def failure(msg):
    print "Terminating - Error "+str(msg)
    print "Should be sending mail to "+str(ERRMAIL)

def progress(msg):
    print msg

class TestConfiguration(object):
    def __init__(self, name, opts, omp, mpi, binexec, pythonexec):
	self.name=name
	self.opts=opts
	self.omp=omp
	self.mpi=mpi
	self.binexec=binexec
	self.pythonexec=pythonexec

    def getHeader():
	res="#!/bin/bash\n"
	res=res+'MAIL_RECIPIENTS="'+ERRMAIL+'"\n'
	res=res+"function report()\n{\n"
	res=res+"	NOW=`date '+%Y/%m/%d %H:%M'`\n"
	res=res+"cat <<END_MSG | mail -s \"Esys tests `date` $TESTSTATE\" $MAIL_RECIPIENTS\n"
	res=res+"Sucessful configurations:\n"
	res=res+"$SUCCESSFUL\n\n"
	res=res+"Failed on configuration:\n"
	res=res+"$ATTEMPTED\n\n"
	res=res+"Tests ran from $START until $NOW.\n"
	res=res+"Log files can be found in $TOP.\n"
	res=res+"This mail was sent from $SCRIPTNAME, running as $USER on `hostname`.\n"
	res=res+"END_MSG\n"
	res=res+"}\n"
	res=res+"function progress()\n{\n"
	res=res+"  echo $1\n"
	res=res+"  cat $1 >> $LOGFILE\n"
	res=res+"}\n"
	res=res+"function failure()\n{\n  echo $1\n"
	res=res+"  report\n"
	res=res+"  exit 1\n}\nTOP=`pwd`\nLOGFILE=$TOP/Logs\nOLDPYTH=$PYTHONPATH\nOLDLD=$LD_LIBRARY_PATH\n"
	res=res+". /usr/share/modules/init/sh		#So the module command works\n"
	res=res+"module load subversion-1.3.1\nmodule load escript/current\nmodule load pbs\nmodule load mayavi/gcc-4.1.2/mayavi-1.5\n"
	res=res+"module load mplayer/gcc-4.1.2/mplayer-1.0rc2\n\n"
	return res

    def toString(self):
	runcount=1
	ref="cp -r "+self.name+"_src "+self.name+"_test"+str(runcount)+"\n" 
	ref=ref+"cd "+self.name+"_test"+str(runcount)+"\n"
	ref=ref+"TESTROOT=`pwd`\n"
	for o in self.omp:
	    for m in self.mpi:
		cmd="bash utest.sh 'mpiexec -np"+str(m)+"' $TESTROOT/lib/pythonMPI"
		ref=ref+"export OMP_NUM_THREADS="+str(o)+"\n"
		ref=ref+"export PYTHONPATH=`pwd`:$OLDPYTH\n"
		ref=ref+"export LD_LIBRARY_PATH=`pwd`/lib:$OLDLD\n"
		ref=ref+'RUNNAME="'+self.name+' omp='+str(o)+' mpi='+str(m)+'"\n'
		ref=ref+'ATTEMPTING=$RUNNAME'
		ref=ref+'progress "Starting '+cmd+'"'
		ref=ref+cmd+' || failure "'+cmd+'"\n'
		ref=ref+'SUCCESSFUL="$SUCCESSFUL, $RUNNAME"\n'
		ref=ref+'completed "'+cmd+'"'
		ref=ref+'ATTEMPTING=None\n'
		ref=ref+"export OMP_NUM_THREADS=1\n"
		++runcount
	    if len(self.mpi)==0:
		cmd="bash utest.sh '' python"
		ref=ref+"export OMP_NUM_THREADS="+str(o)+"\n"
		ref=ref+"export LD_LIBRARY_PATH=`pwd`/lib:$OLDLD\n"
		ref=ref+"export PYTHONPATH=`pwd`:$OLDPYTH\n"
		ref=ref+'progress "Starting '+cmd+'"' 
		ref=ref+cmd+" || failure \""+cmd+"\" \n" 
		ref=ref+'completed "'+cmd+'"'
		++runcount
	ref=ref+"\ncd $TOP\n\n"
	return ref
 
    getHeader=staticmethod(getHeader)


try:
    os.mkdir(TOPDIR)
    os.chdir(TOPDIR)
except OSError:
	failure("Unable to create top directory "+TOPDIR+" does it exist already?")
	sys.exit(1)

try:
    os.mkdir("Logs")
except OSError:
	failure("Unable to create Logs directory ")
	sys.exit(1)

coresult=os.system("svn export "+SVNURL+" src")
if coresult!=0:
    failure("Unable to export working copy")
    sys.exit(1)

testconfs=[]
testconfs.append(TestConfiguration("OMPNoMPI","",omp=(1,8),mpi=(),binexec="",pythonexec="python"))
testconfs.append(TestConfiguration("MPI","usempi=yes",omp=(1,),mpi=(1,8),binexec="mpiexec -np ",pythonexec="lib/pythonMPI"))


dir=os.getcwd()
for conf in testconfs:
    progress("Creating "+conf.name+"_src")
# Yes I know copytree exists. No I don't trust it (It was having problems doing the copy).
    res=os.system("cp -r src "+conf.name+"_src")
    if res!=0:
	failure("Error copying src to "+conf.name+"_src")
#    shutil.copytree("src",conf.name+"_src")
    os.chdir(conf.name+"_src")
    cmdstr="scons -j"+str(NUMJs)+" "+conf.opts+" install_all build_tests build_py_tests"
    progress(cmdstr)
    res=os.system(cmdstr)
    os.chdir(str(dir))
    if res!=0:
	failure("Error running scons build failed for "+conf.name+"_src")
    
progress("Builds complete")
#print "Removing export copy"
#shutil.rmtree("src",ignore_errors=True)
progress("Building test file")

try:
	testfile=open("dotests.sh","w")
	testfile.write(TestConfiguration.getHeader())
	for c in testconfs:
		testfile.write(c.toString())
	testfile.close()
except IOError:
	failure("Creating testfile")
	
progress("Building test file complete")

print "Should be submitting this test now"
raise "Test not submitted"
#Now we build the script and submit it pbs (that behaviour should be optional?)
