
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#ifndef OVERLORDPATH
  #define OVERLORDPATH ""
#endif

#include <fstream>
#include <string.h>

// added for saveCSV
#include <boost/python.hpp>
#include <boost/scoped_array.hpp>
#include "Data.h"

#include "Utils.h"
#include "DataVector.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ESYS_MPI
#include <mpi.h>
#endif

#ifdef  _WIN32
#include <WinSock2.h>
#else
#include <unistd.h>
#endif

namespace escript {

int getSvnVersion() 
{
#ifdef SVN_VERSION
  return SVN_VERSION;
#else
  return 0;
#endif
}

/* This is probably not very robust, but it works on Savanna today and is useful for performance analysis */
int get_core_id() {
  int processor_num=-1;
#ifdef CORE_ID1
  FILE *fp;
  int i, count_spaces=0;
  char fname[100];
  char buf[1000];

  sprintf(fname, "/proc/%d/stat", getpid());
  fp = fopen(fname, "r");
  if (fp == NULL) return(-1);
  fgets(buf, 1000, fp);
  fclose(fp);

  for (i=strlen(buf)-1; i>=0; i--) {
    if (buf[i] == ' ') count_spaces++;
    if (count_spaces == 4) break;
  }
  processor_num = atoi(&buf[i+1]);
#endif
  return(processor_num);
}


void printParallelThreadCnt() 
{
  int mpi_iam=0, mpi_num=1;
  char hname[64];

#ifdef HAVE_GETHOSTNAME
  gethostname(hname, 64);
  hname[63] = '\0';
#else
  strcpy(hname, "unknown host");
#endif

  #ifdef ESYS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif

  #pragma omp parallel
  {
    int omp_iam=0, omp_num=1;
    #ifdef _OPENMP
    omp_iam = omp_get_thread_num(); /* Call in a parallel region */
    omp_num = omp_get_num_threads();
    #endif
    #pragma omp critical (printthrdcount)
    printf("printParallelThreadCounts: MPI=%03d/%03d OpenMP=%03d/%03d running on %s core %d\n",
      mpi_iam, mpi_num, omp_iam, omp_num, hname, get_core_id());
  }
}

void setNumberOfThreads(const int num_threads) 
{

   #ifdef _OPENMP
   omp_set_num_threads(num_threads);
   #endif

}

int getNumberOfThreads() 
{
   #ifdef _OPENMP
   return omp_get_max_threads();
   #else
   return 1;
   #endif

}

ESCRIPT_DLL_API int getMPISizeWorld() {
  int mpi_num = 1;
  #ifdef ESYS_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif
  return mpi_num;
}

ESCRIPT_DLL_API int getMPIRankWorld() {
  int mpi_iam = 0;
  #ifdef ESYS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  #endif
  return mpi_iam;
}

ESCRIPT_DLL_API int getMPIWorldMax(const int val) {
  #ifdef ESYS_MPI
  int val2 = val;
  int out = val;
  MPI_Allreduce( &val2, &out, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  #else
  int out = val;
  #endif
  return out;
}

ESCRIPT_DLL_API int getMPIWorldSum(const int val) {
  #ifdef ESYS_MPI
  int val2 = val;
  int out = 0;
  MPI_Allreduce( &val2, &out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  #else
  int out = val;
  #endif
  return out;
}

#define CHILD_FAIL 2
#define CHILD_COMPLETE 4

#ifndef _WIN32
#ifdef ESYS_MPI
#include <sys/socket.h>
#include <sys/select.h>
#include <errno.h>
#include <arpa/inet.h>
int prepareSocket(unsigned short *port, int *key) {
    if (getMPIRankWorld() != 0)
        return 0;
    int sfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sfd < 0) {
        perror("socket creation failure");
        return -1;
    }
    int opt = 1;
    if (setsockopt(sfd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(int)) < 0) {
        perror("socket option setting failure");
        close(sfd);
        return -1;
    }

    struct sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(0);
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);

    if (bind(sfd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("bind failure");
        close(sfd);
        return -1;
    }

    if (listen(sfd, SOMAXCONN) < 0) {
        perror("listen failure");
        close(sfd);
        return -1;
    }
    
    struct sockaddr actual;
    unsigned int size = sizeof(actual);
    if (getsockname(sfd, &actual, &size) < 0) {
        perror("failed when determining bound port number");
        close(sfd);
        return -1;
    }
    
    //if size > sizeof(actual), some info was truncated, but certainly not port
    
    *port = ntohs(((struct sockaddr_in *) &actual)->sin_port);

    unsigned int seed = time(NULL) % UINT_MAX;
    *key = rand_r(&seed);
    return sfd;
}

int check_data(unsigned int max, fd_set *all, fd_set *valid, int key, int sfd) {
    int ret = 0;
    for (int i = 0; i <= max; i++) {
        if (i == sfd)
            continue;
        if (FD_ISSET(i, all)) {
            int provided = 0;
            if (recv(i, &provided, sizeof(int), MSG_WAITALL) == sizeof(int)
                    && provided == key) {
                char deadspace[1024];
                while ((provided = recv(i, deadspace, 1024, 0))) {
                    if (provided == -1) {
                        if (errno == EAGAIN || errno == EWOULDBLOCK) {
                            continue;
                        } else {
                            perror("connection failure");
                            return CHILD_FAIL;
                        }
                    }
                }
                return CHILD_COMPLETE;
            } else {
                FD_CLR(i, all);
                close(i);
            }
        }
    }
    return ret;
}

void close_all(unsigned int maxfd, fd_set *all) {
    for (int i = 0; i <= maxfd; i++) {
        if (FD_ISSET(i, all))
            close(i);
    }
}
int waitForCompletion(int sfd, int key) {
    if (getMPIRankWorld() != 0)
        return 0;
    int timeout = 10; //max of 10 seconds with no communication
    
    fd_set all, valid; 
    FD_ZERO(&all);
    FD_SET(sfd, &all);
    time_t last_good_time = time(NULL);
    unsigned int maxfd = sfd;
    while (time(NULL) - last_good_time < timeout) {
        struct timeval timer = {1,0}; //1 sec, 0 usec
        int count = select(maxfd + 1, &all, NULL, NULL, &timer);
        if (count == -1) { //error
            if (errno == EINTR) {
                continue; //just a signal, continue as we were
            } else {
                perror("socket operation error");
                close(sfd);
                return -1;
            }
        } else if (FD_ISSET(sfd, &all)) { //new connection
            int connection = accept(sfd, NULL, NULL);
            if (connection > maxfd)
                maxfd = connection;
            FD_SET(connection, &all);
            FD_CLR(connection, &valid);
            time(&last_good_time);
            count--;
        }
        if (count > 0) { //something to read, either connection key or state
            int res = check_data(maxfd, &all, &valid, key, sfd);
            if (res == CHILD_FAIL) {
                return -1;
            } else if (res == CHILD_COMPLETE) {
                close_all(maxfd, &all);
                return 0;
            }
        }
    }
    close_all(maxfd, &all);
    fprintf(stderr, "Connection to child process timed out\n");
    return -1;
}
#endif //ESYS_MPI
#endif //not _WIN32

ESCRIPT_DLL_API int runMPIProgram(boost::python::list args) {
#ifdef ESYS_MPI
    MPI_Comm intercomm;
    MPI_Info info;
    int errors;
    int nargs = boost::python::extract<int>(args.attr("__len__")());
    std::string cmd = boost::python::extract<std::string>(args[0]);
#ifdef _WIN32
    char** c_args = new char*[nargs];
    char* c_cmd = const_cast<char*>(cmd.c_str());;
    // skip command name in argument list
    for (int i=1; i<nargs; i++) {
	    cpp_args[i-1]=boost::python::extract<std::string>(args[i]);
        c_args[i-1]=const_cast<char*>(cpp_args[i-1].c_str());
    }
    MPI_Info_create(&info);
    MPI_Comm_spawn(c_cmd, c_args, 1, info, 0, MPI_COMM_WORLD, &intercomm, &errors);
    MPI_Info_free(&info);
    delete[] c_args;

    return errors;
#else //#ifdef _WIN32
    char** c_args = new char*[2+nargs];
    std::vector<std::string> cpp_args(2+nargs);//allow for wrapper, port, and key
    char c_cmd[] = OVERLORDPATH"escript-overlord";
    // skip command name in argument list
    for (int i=1; i<nargs; i++) {
	    cpp_args[i+1]=boost::python::extract<std::string>(args[i]);
        c_args[i+1]=const_cast<char*>(cpp_args[i+1].c_str());
    }
    unsigned short port = 0;
    int key = 0;
    int sock = prepareSocket(&port, &key);
    if (getMPIWorldSum(sock) < 0)
        return -1;
    c_args[nargs+1]=NULL;
    char portstr[20] = {'\0'}, keystr[20] = {'\0'};
    sprintf(portstr, "%d", port);
    sprintf(keystr, "%d", key);
    c_args[0] = portstr;
    c_args[1] = keystr;
    c_args[2] = const_cast<char*>(cmd.c_str());
    
    MPI_Info_create(&info);
    //force the gmsh process to run on this host as well for network comm
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int temp = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(hostname, &temp);
    char hoststr[] = "host"; //for warnings
    MPI_Info_set(info, hoststr, hostname);
    MPI_Comm_spawn(c_cmd, c_args, 1, info, 0, MPI_COMM_WORLD, &intercomm, &errors);
    MPI_Info_free(&info);
    delete[] c_args;
    if (errors != MPI_SUCCESS)
        return errors;
    return getMPIWorldMax(waitForCompletion(sock, key));
#endif //#ifdef _WIN32/else
#else //#ifdef ESYS_MPI
    std::string cmd;
    int nargs = boost::python::extract<int>(args.attr("__len__")());
    for (int i=0; i<nargs; i++) {
	    cmd+=boost::python::extract<std::string>(args[i]);
        cmd+=" ";
    }
    return system(cmd.c_str());
#endif //#ifdef ESYS_MPI/else
}
#undef CHILD_COMPLETE
#undef CHILD_FAIL

ESCRIPT_DLL_API double getMachinePrecision() {
   return DBL_EPSILON;
}
ESCRIPT_DLL_API double getMaxFloat() {
   return DBL_MAX;
}
ESCRIPT_DLL_API void MPIBarrierWorld() {
  #ifdef ESYS_MPI
  if (!esysUtils::NoCOMM_WORLD::active())
  {
      MPI_Barrier(MPI_COMM_WORLD );
  }
  else
  {
      throw esysUtils::EsysException("Attempt to use MPI_COMM_WORLD while it is blocked.");
  }
  #endif
}


ESCRIPT_DLL_API
void
saveDataCSV(const std::string& filename, boost::python::dict arg, const std::string& sep, const std::string& csep, 
bool append)
{
    using std::cout;
    using std::endl;
    boost::python::list keys=arg.keys();
    int numdata = boost::python::extract<int>(arg.attr("__len__")());
    bool hasmask=arg.has_key("mask");
    Data mask;
    if (hasmask)
    {
	mask=boost::python::extract<escript::Data>(arg["mask"]);
	keys.remove("mask");
	numdata--;
        if (mask.getDataPointRank()!=0)
	{
		throw DataException("saveDataCSVcpp: masks must be scalar.");
	}
    }
    if (numdata<1)
    {
	throw DataException("saveDataCSVcpp: no data to save specified.");
    }
    std::vector<int> step(numdata);
    std::vector<std::string> names(numdata);
    std::vector<Data> data(numdata);
    std::vector<const DataAbstract::ValueType::value_type*> samples(numdata);
    std::vector<int> offset(numdata);
    std::vector<int> fstypes(numdata);		// FunctionSpace types for each data

    keys.sort();	// to get some predictable order to things

    // We need to interpret the samples correctly even if they are different types
    // for this reason, we should iterate over samples
    for (int i=0;i<numdata;++i)
    {
	names[i]=boost::python::extract<std::string>(keys[i]);
	data[i]=boost::python::extract<escript::Data>(arg[keys[i]]);
	step[i]=(data[i].actsExpanded()?DataTypes::noValues(data[i].getDataPointShape()):0);
	fstypes[i]=data[i].getFunctionSpace().getTypeCode();
	if (i>0) 
	{
	    if (data[i].getDomain()!=data[i-1].getDomain())
	    {
		throw DataException("saveDataCSVcpp: all data must be on the same domain.");
	    }
	}
    }
    int bestfnspace=0;
    if (!data[0].getDomain()->commonFunctionSpace(fstypes, bestfnspace))
    {
	throw DataException("saveDataCSVcpp: FunctionSpaces of data are incompatible");
    }
    // now we interpolate all data to the same type
    FunctionSpace best(data[0].getDomain(),bestfnspace);
    for (int i=0;i<numdata;++i)
    {
	data[i]=data[i].interpolate(best);
    }
    int numsamples=data[0].getNumSamples();		// these must be the same for all data
    int dpps=data[0].getNumDataPointsPerSample();

    
    std::ostringstream os;


    bool first=true;
    
    if (data[0].getDomain()->getMPIRank()==0)
    {
      for (int i=0;i<numdata;++i)
      {
	const DataTypes::ShapeType& s=data[i].getDataPointShape();
        switch (data[i].getDataPointRank())
	{
	case 0: if (!first)
		    {
			os << sep;
		    }
		    else
		    {
			first=false;
		    }
		os << names[i]; break;
	case 1: for (int j=0;j<s[0];++j)
		{
		    if (!first)
		    {
			os << sep;
		    }
		    else
		    {
			first=false;
		    }
		    os << names[i] << csep << j;
		}
		break;
	case 2: for (int j=0;j<s[0];++j)
		{
		    for (int k=0;k<s[1];++k)
		    {
			if (!first)
			{
				os << sep;
			}
			else
			{
				first=false;
			}
		    	os << names[i] << csep << k << csep << j;
		    }
		}
		break;
	case 3: for (int j=0;j<s[0];++j)
		{
		    for (int k=0;k<s[1];++k)
		    {
			for (int l=0;l<s[2];++l)
			{
				if (!first)
				{
					os << sep;
				}
				else
				{
					first=false;
				}
				os << names[i] << csep << k << csep << j << csep << l;
			}
		    }
		}
		break;
	case 4: for (int j=0;j<s[0];++j)
		{
		    for (int k=0;k<s[1];++k)
		    {
			for (int l=0;l<s[2];++l)
			{
			    for (int m=0;m<s[3];++m)
			    {
				if (!first)
				{
					os << sep;
				}
				else
				{
					first=false;
				}
				os << names[i] << csep << k << csep << j << csep << l << csep << m;
			    }
			}
		    }
		}
		break;
	default:
		throw DataException("saveDataCSV: Illegal rank");
	}
      }
      os << endl; 
    }

    const double* masksample=0;

    bool expandedmask=false;		// does the mask act expanded. Are there mask value for each point in the sample
    bool wantrow=true;			// do we output this row?
    if (hasmask)
    {
	if (mask.actsExpanded())
	{
		expandedmask=true;
	}
    }
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os.precision(15);

    // errors prior to this point will occur on all processes anyway
    // so there is no need to explicitly notify other ranks
    int error=0;
    try{
      for (int i=0;i<numsamples;++i)
      {
        if (!best.ownSample(i))
	{
		continue;
	}
	wantrow=true;
	for (int d=0;d<numdata;++d)
	{
	  	samples[d]=data[d].getSampleDataRO(i);
	}
	if (hasmask)
	{
		masksample=mask.getSampleDataRO(i);
		if (!expandedmask)		// mask controls whole sample
		{
			if (masksample[0]<=0)		// masks are scalar
			{
				wantrow=false;
			}
		}
	}
	for (int j=0;j<dpps;++j)
	{
	    // now we need to check if this point is masked off
	    if (expandedmask)
	    {
		wantrow=(masksample[j]>0); // masks are scalar to the relevant value is at [j]
	    }
	    if (wantrow)
	    {
		bool needsep=false; 
		for (int d=0;d<numdata;++d)
		{
			DataTypes::pointToStream(os, samples[d], data[d].getDataPointShape(), offset[d], needsep, sep);
			needsep=true;
			offset[d]+=step[d];
		}
		os << endl;
	    }
	}
	for (int d=0;d<numdata;++d)
	{
	    offset[d]=0;
	}	
      }
    } catch (...)
    {
	error=1;
	if (data[0].getDomain()->getMPISize()==1) {
            throw;
        }
    }
#ifdef ESYS_MPI
    MPI_Comm com=data[0].getDomain()->getMPIComm();
    int rerror=0;
    MPI_Allreduce( &error, &rerror, 1, MPI_INT, MPI_MAX, com );
    error=rerror;
    if (error)
    {
	throw DataException("saveDataCSVcpp: error building output");
    }
#endif

    // at this point os will contain the text to be written
#ifndef ESYS_MPI
    (void) error;

    std::ofstream ofs;
    if (append)
    {
	ofs.open(filename.c_str(), std::ios_base::app);
    }
    else
    {
	ofs.open(filename.c_str());
    }
    if (!ofs.is_open())
    {
	throw DataException("saveDataCSVcpp: unable to open file for writing");
    }
    ofs << os.str();
    ofs.close();

#else
// here we have MPI
    MPI_File mpi_fileHandle_p;
    MPI_Status mpi_status;
    MPI_Info mpi_info = MPI_INFO_NULL;
    boost::scoped_array<char> fname_p(new char[filename.size()+1]);
    strcpy(fname_p.get(), filename.c_str());
     
    int amode = MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_UNIQUE_OPEN;
    if (append)
    {
	amode |= MPI_MODE_APPEND;
    }
    else
    {
	if (data[0].getDomain()->getMPIRank()==0)
	{
	    std::ifstream ifs(fname_p.get());	// if file exists, remove it
	    if (ifs.is_open())
	    {
		ifs.close();
	    	if (remove(fname_p.get()))
		{
		    error=1;
		}
	    }
	}
	data[0].getDomain()->MPIBarrier();
	int rerror=0;
	MPI_Allreduce( &error, &rerror, 1, MPI_INT, MPI_MAX, com );
	if (rerror!=0)
	{
	    std::ostringstream oss;
	    oss << "saveDataCSVcpp: File " << filename << " already exists and could not be removed in preparation for new output.";
	    throw DataException(oss.str());
	}
    }
    int ierr;
    ierr = MPI_File_open(com, fname_p.get(), amode, mpi_info, &mpi_fileHandle_p);
    if (ierr != MPI_SUCCESS) 
    {
	std::ostringstream oss;
	oss << "saveDataCSVcpp: File " << filename << " could not be opened for writing in parallel";
	    // file is not open so we can throw
	throw DataException(oss.str());
    }
    else
    {
            ierr=MPI_File_set_view(mpi_fileHandle_p,MPI_DISPLACEMENT_CURRENT,
                    MPI_CHAR, MPI_CHAR, const_cast<char*>("native"), mpi_info);
// here we are assuming that std::string holds the same type of char as MPI_CHAR
    }

    std::string contents=os.str();
    boost::scoped_array<char> buff(new char[contents.size()+1]);
    strcpy(buff.get(), contents.c_str());
    ierr=MPI_File_write_ordered(mpi_fileHandle_p, buff.get(), contents.size(), MPI_CHAR, &mpi_status);
    if (ierr != MPI_SUCCESS)
    {
	error=1;
    }

    if (MPI_File_close(&mpi_fileHandle_p)!= MPI_SUCCESS)
    {
	error=1;
    }
    data[0].getDomain()->MPIBarrier();
    if (error)	// any errors at this stage are from collective routines
    {		// so there is no need to reduce_all
	throw DataException("saveDataCSVcpp: Error writing and closing file");
    }
    
#endif
}


void 
resolveGroup(boost::python::object obj)
{
	int len=0;
	try
	{
	   len=boost::python::extract<int>(obj.attr("__len__")());
	}
	catch(...)
	{
	   PyErr_Clear();		// tell python the error isn't there anymore
	   throw DataException("Error - resolveGroup expects a sequence object.");
	}
	std::vector<DataLazy*> dats;
	std::vector<Data*> dp;
	for (int i=0;i<len;++i)
	{
		Data* p=0;
		try
		{
			p=boost::python::extract<Data*>(obj[i]);
		}
		catch(...)
		{
			PyErr_Clear();
			throw DataException("Error - resolveGroup only accepts Data objects.");
		}
		if (p->isLazy())
		{
			dats.push_back(dynamic_cast<DataLazy*>(p->borrowData()));
			dp.push_back(p);
		}
	}
	if (!dats.empty())
	{
		dats[0]->resolveGroupWorker(dats);
	}
	// all the data will be identities now but its still lazy
	// convert it to ready
	for (int i=dp.size()-1;i>=0;--i)
	{
		dp[i]->resolve();
	}
}


}  // end of namespace
