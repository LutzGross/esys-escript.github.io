
import shutil, os, datetime, sys, os.path, time

#This script does not use the python, platform independent path manipulation stuff.
#It probably should

SVNURL="https://shake200.esscc.uq.edu.au/svn/esys13/trunk"
NUMJs=4
TOPDIR=str(datetime.date.today())
ERRMAIL="j.fenwick1@uq.edu.au"
EXECUTELOCATION="/scratch/jfenwick/AUTOTESTS"
OUTSIDEDIR=os.getcwd()
TESTSLEEP=30*60

SRCMSG="This message was sent by prepare.py running as "+str(os.environ['USER'])+" on "+str(os.environ['HOSTNAME']+"\n")

#Settings for actual tests appear below the class declarations


def failure(msg):
    print "Terminating - Error "+str(msg)
    print "Should be sending mail to "+str(ERRMAIL)
    mailcmd="cat << ENDMSG |mail -s 'Esys unit tests failed to execute properly' "+ERRMAIL+"\n"
    mailcmd=mailcmd+"Error preparing for test run:\n"+msg+"\n" 
    mailcmd=mailcmd+SRCMSG
    mailcmd=mailcmd+"ENDMSG\n"
    os.system(mailcmd)
    sys.exit(1)

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
	res=res+"	cat > $LOGDIR/message << END_MSG\n"
	res=res+"Sucessful configurations:\n"
	res=res+"$SUCCESSFUL\n\n"
	res=res+"Failed on configuration:\n"
	res=res+"$ATTEMPTING\n\n"
	res=res+"Tests ran from $START until $NOW.\n"
	res=res+"Log files can be found in $FINALLOGDIR.\n"
	res=res+"END_MSG\n"
	res=res+"}\n"
	res=res+"function progress()\n{\n"
	res=res+"  echo $1\n"
	res=res+"  echo $1 >> $PROGRESSFILE\n"
	res=res+"}\n"
	res=res+"function failure()\n{\n  echo $1\n"
	res=res+"  report\n"
	res=res+"  touch $LOGDIR/Failure\n"
	res=res+"  if [ -f stdout_cpu_0001.out ];then cp std_cpu_* $TESTLOGDIR;fi\n"
	res=res+"  exit 1\n}\n"
	res=res+"cd "+EXECUTELOCATION+"/"+TOPDIR+"\n"
	res=res+"TOP=`pwd`\nLOGDIR=$TOP/Logs\nPROGRESSFILE=$LOGDIR/progress\nOLDPYTH=$PYTHONPATH\nOLDLD=$LD_LIBRARY_PATH\n"
	res=res+". /usr/share/modules/init/sh		#So the module command works\n"
	res=res+"module load subversion-1.3.1\nmodule load escript/current\nmodule load pbs\nmodule load mayavi/gcc-4.1.2/mayavi-1.5\n"
	res=res+"module load mplayer/gcc-4.1.2/mplayer-1.0rc2\n\n"
	res=res+"SCRIPTNAME=$0\n"
	res=res+"START=`date '+%Y/%m/%d %H:%M'`\n"
	res=res+"TESTLOGDIR=$LOGDIR\n"
	res=res+"FINALLOGDIR="+OUTSIDEDIR+"/"+TOPDIR+"_Logs\n"
	return res

    def toString(self):
	runcount=1
	res=""
	print "Processing "+self.name
	for o in self.omp:
	    print "o="+str(o)
	    for m in self.mpi:
		print "   m="+str(m)		
		cmd="bash utest.sh 'mpiexec -np"+str(m)+"' $TESTROOT/lib/pythonMPI  >$TESTLOGDIR/output 2>&1"
		res=res+"cp -r "+self.name+"_src "+self.name+"_test"+str(runcount)+"\n" 
		res=res+"cd "+self.name+"_test"+str(runcount)+"\n"
		res=res+"TESTROOT=`pwd`\n"
		res=res+"TESTLOGDIR=$LOGDIR/"+self.name+"_test"+str(runcount)+"\n"
		res=res+"mkdir $TESTLOGDIR\n"
		res=res+"export OMP_NUM_THREADS="+str(o)+"\n"
		res=res+"export PYTHONPATH=`pwd`:$OLDPYTH\n"
		res=res+"export LD_LIBRARY_PATH=`pwd`/lib:$OLDLD\n"
		res=res+'RUNNAME="'+self.name+' omp='+str(o)+' mpi='+str(m)+'"\n'
		res=res+'ATTEMPTING=$RUNNAME\n'
		res=res+'progress "Starting '+cmd+'"\n'
		res=res+cmd+' || failure "'+cmd+'"\n'
		res=res+"if [ -f stdout_cpu_0001.out ];then cp std_cpu_* $TESTLOGDIR;fi\n"
		res=res+'SUCCESSFUL="$SUCCESSFUL, $RUNNAME"\n'
		res=res+'progress "completed '+cmd+'"\n'
		res=res+'ATTEMPTING=None\n'
		res=res+"export OMP_NUM_THREADS=1\n"
		res=res+"cd $TOP\n"
		res=res+"rm -rf "+self.name+"_test"+str(runcount)+"\n\n"
		runcount=runcount+1
	    if len(self.mpi)==0:
		print "   m=()"
		cmd="bash utest.sh '' python $TESTROOT/lib/pythonMPI  >$TESTLOGDIR/output 2>&1"
		res=res+"cp -r "+self.name+"_src "+self.name+"_test"+str(runcount)+"\n" 
		res=res+"cd "+self.name+"_test"+str(runcount)+"\n"
		res=res+"TESTROOT=`pwd`\n"
		res=res+"TESTLOGDIR=$LOGDIR/"+self.name+"_test"+str(runcount)+"\n"
		res=res+"mkdir $TESTLOGDIR\n"
		res=res+"export OMP_NUM_THREADS="+str(o)+"\n"
		res=res+"export LD_LIBRARY_PATH=`pwd`/lib:$OLDLD\n"
		res=res+"export PYTHONPATH=`pwd`:$OLDPYTH\n"
		res=res+'RUNNAME="'+self.name+' omp='+str(o)+' mpi=n/a"\n'
		res=res+'ATTEMPTING=$RUNNAME\n'
		res=res+'progress "Starting '+cmd+'"\n' 
		res=res+cmd+" || failure \""+cmd+"\" \n"
		res=res+"if [ -f stdout_cpu_0001.out ];then cp std_cpu_* $TESTLOGDIR;fi\n"
		res=res+'ATTEMPTING=None\n'
		res=res+'progress "completed '+cmd+'"\n'
		res=res+"cd $TOP\n"
		res=res+"rm -rf "+self.name+"_test"+str(runcount)+"\n\n"
		runcount=runcount+1
	res=res+"rm -rf "+self.name+"_src\n"
	res=res+"\ncd $TOP\n\n"
	return res
 
    def getFooter():
	res="\ntouch $LOGDIR/Success\n"
	res=res+"report"
	return res

    getHeader=staticmethod(getHeader)
    getFooter=staticmethod(getFooter)

#Test settings
testconfs=[]
testconfs.append(TestConfiguration("OMPNoMPI","",omp=(1,8),mpi=(),binexec="",pythonexec="python"))
testconfs.append(TestConfiguration("MPI","usempi=yes",omp=(1,),mpi=(1,8),binexec="mpiexec -np ",pythonexec="lib/pythonMPI"))

LOGDIR=OUTSIDEDIR+"/"+TOPDIR+"_Logs"

if os.path.exists(LOGDIR):
	failure("Logs directory for "+TOPDIR+" already exists.")
	sys.exit(1)

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




dir=os.getcwd()
for conf in testconfs:
    progress("Creating "+conf.name+"_src")
#    res=os.system("cp -r src "+conf.name+"_src")
#    if res!=0:
#	failure("Error copying src to "+conf.name+"_src")
    try:
    	shutil.copytree("src",conf.name+"_src")
    except Error:
	failure("copying src to "+conf.name+"_src")
    os.chdir(conf.name+"_src")
    cmdstr="scons -j"+str(NUMJs)+" "+conf.opts+" install_all build_tests build_py_tests"
    progress(cmdstr)
    res=os.system(cmdstr)
    os.chdir(str(dir))
    if res!=0:
	failure("running scons build failed for "+conf.name+"_src")
    
progress("Builds complete")
progress("Removing export copy")
shutil.rmtree("src",ignore_errors=True)
progress("Building test file")

try:
	testfile=open("dotests.sh","w")
	testfile.write(TestConfiguration.getHeader())
	for c in testconfs:
		testfile.write(c.toString())
	testfile.write(TestConfiguration.getFooter())
	testfile.close()
	import stat
	os.chmod("dotests.sh",stat.S_IEXEC|stat.S_IREAD)
except IOError:
	failure("Creating testfile")

progress("Building test file complete")
progress("Copying files to exec area")
os.chdir(OUTSIDEDIR)
try:
	shutil.copytree(TOPDIR,EXECUTELOCATION+"/"+TOPDIR)
except Error:
	failure("copying to work area")
progress("Copy to exec area complete")

print "Submitting test"

######### test section

os.chdir(EXECUTELOCATION)
os.chdir(TOPDIR)


try:
#	res=os.system("bash dotests.sh")
	res=os.system("qsub -l select=1:ncpus=8:mem=31gb dotests.sh")
except OSError:
	failure("Submitting tests")

os.chdir(OUTSIDEDIR)

print "Sleeping for "+str(TESTSLEEP)+" seconds to wait for results."
time.sleep(TESTSLEEP)
print "Waking up."

########################

try:
	shutil.copytree(EXECUTELOCATION+"/"+TOPDIR+"/Logs",OUTSIDEDIR+"/"+TOPDIR+"_Logs")
except OSError:
	failure("Log copy failed")

cleanupfailure=False

try:
	shutil.rmtree(EXECUTELOCATION+"/"+TOPDIR)
	shutil.rmtree(OUTSIDEDIR+"/"+TOPDIR)
except OSError:
	cleanupfailure=True

#Now we sum up
#if Logs/Failure or Logs/Success does not exist send message about tests not completing.
if not os.path.exists(LOGDIR+"/Success") and  not os.path.exists(LOGDIR+"/Failure"):
	mailcmd="cat << ENDMSG |mail -s 'Esys unit tests failed to execute properly' "+ERRMAIL+"\n"
	mailcmd=mailcmd+"For some reason no Success or failure is recorded for unit tests in "+LOGDIR+"\n" 
	mailcmd=mailcmd+"\nAlso: the cleanup of work areas failed. "+EXECUTELOCATION+"/"+TOPDIR+" or "+OUTSIDEDIR+"/"+TOPDIR+"\n"
	mailcmd=mailcmd+SRCMSG
	mailcmd=mailcmd+"ENDMSG\n"
	os.system(mailcmd)
	sys.exit(1)


if os.path.exists(LOGDIR+"/Failure"):
	os.system("echo 'Also: the cleanup of work areas failed. "+EXECUTELOCATION+"/"+TOPDIR+" or "+OUTSIDEDIR+"/"+TOPDIR+"' >> "+LOGDIR+"/message\n")
	os.system("echo '"+SRCMSG+"' >> "+LOGDIR+"/message\n")
	mailcmd="cat "+LOGDIR+"/message | mail -s 'Esys unit tests failed' "+ERRMAIL+"\n"
	os.system(mailcmd)
	sys.exit(1)
#so we must have succeeded

if cleanupfailure:
	os.system("echo '"+SRCMSG+"' >> "+LOGDIR+"/message\n")
	mailcmd="cat "+LOGDIR+"/message | mail -s 'Esys unit tests cleanup failed - tests succeeded' "+ERRMAIL+"\n"
	res=os.system(mailcmd)
	sys.exit(res)
else:
	os.system("echo '"+SRCMSG+"' >> "+LOGDIR+"/message\n")	
	mailcmd="cat "+LOGDIR+"/message | mail -s 'Esys unit tests succeeded' "+ERRMAIL+"\n"
	res=os.system(mailcmd)
	sys.exit(res)

