
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <fstream>
#include <string.h>

// added for saveCSV
#include <boost/python.hpp>
#include <boost/scoped_ptr.hpp>
#include "Data.h"

#include "Utils.h"
#include "DataVector.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PASO_MPI
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

  #ifdef PASO_MPI
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
  #ifdef PASO_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif
  return mpi_num;
}

ESCRIPT_DLL_API int getMPIRankWorld() {
  int mpi_iam = 0;
  #ifdef PASO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  #endif
  return mpi_iam;
}

ESCRIPT_DLL_API int getMPIWorldMax(const int val) {
  #ifdef PASO_MPI
  int val2 = val;
  int out = val;
  MPI_Allreduce( &val2, &out, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  #else
  int out = val;
  #endif
  return out;
}

ESCRIPT_DLL_API int getMPIWorldSum(const int val) {
  #ifdef PASO_MPI
  int val2 = val;
  int out = 0;
  MPI_Allreduce( &val2, &out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  #else
  int out = val;
  #endif
  return out;
}

ESCRIPT_DLL_API double getMachinePrecision() {
   return DBL_EPSILON;
}
ESCRIPT_DLL_API double getMaxFloat() {
   return DBL_MAX;
}
ESCRIPT_DLL_API void MPIBarrierWorld() {
  #ifdef PASO_MPI
  MPI_Barrier(MPI_COMM_WORLD );
  #endif
}


/*
ESCRIPT_DLL_API
void
saveDataCSV(const std::string& filename, boost::python::dict arg, const std::string& sep, const std::string& csep, 
bool append)
{
    boost::python::list keys=arg.keys();
    int numdata = boost::python::extract<int>(arg.attr("__len__")());
    bool hasmask=arg.has_key("mask");
    Data mask;
    if (hasmask)
    {
	mask=boost::python::extract<escript::Data>(arg["mask"]);
	keys.remove("mask");
        if (mask.getDataPointRank()!=0)
	{
		throw DataException("saveDataCSVcpp: masks must be scalar.");
	}
    }
    if (numdata<1)
    {
	throw DataException("saveDataCSVcpp: no data to save specified.");
    }
    std::vector<std::string> names(numdata);
    std::vector<Data> data(numdata);
    std::vector<int> fstypes(numdata);		// FunctionSpace types for each data
 
    if (hasmask)
    {
	numdata--;
    }

    // We need to interpret the samples correctly even if they are different types
    // for this reason, we should interate over samples
    for (int i=0;i<numdata;++i)
    {
	names[i]=boost::python::extract<std::string>(keys[i]);
	data[i]=boost::python::extract<escript::Data>(arg[keys[i]]);
	fstypes[i]=data[i].getFunctionSpace().getTypeCode();
	if (i>0) 
	{
	    if (data[i].getDomain()!=data[i-1].getDomain())
	    {
		throw DataException("saveDataCSVcpp: all data must be on the same domain.");
	    }
	}
    }
    const_Domain_ptr dom0=data[0].getDomain();
    if (hasmask)
    {
	if (dom0!=mask.getDomain())
	{
	    throw DataException("saveDataCSVcpp: mask be on the same FunctionSpace as data.");
	}
	fstypes[numdata]=mask.getFunctionSpace().getTypeCode();
	names[numdata]="mask";
	data[numdata]=mask;
    }
    int bestfnspace=0;
    if (!dom0->commonFunctionSpace(fstypes, bestfnspace))
    {
	throw DataException("saveDataCSVcpp: FunctionSpaces of data are incompatible");
    }
    // now we interpolate all data to the same type
    FunctionSpace best(dom0,bestfnspace);
    for (int i=0;i<data.size();++i)
    {
	data[i]=data[i].interpolate(best);
    }
    dom0->saveDataCSV(filename, data, names, sep, csep, append, hasmask);
}
*/


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

    // We need to interpret the samples correctly even if they are different types
    // for this reason, we should interate over samples
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
    boost::scoped_ptr<BufferGroup> maskbuffer;	// sample buffer for the mask [if we have one]
    const double* masksample=0;
    int maskoffset=0;
	//the use of shared_ptr here is just to ensure the buffer group is freed
	//I would have used scoped_ptr but they don't work in vectors
    std::vector<boost::shared_ptr<BufferGroup> > bg(numdata);
    for (int d=0;d<numdata;++d)
    {
	bg[d].reset(data[d].allocSampleBuffer());
    }

    bool expandedmask=false;		// does the mask act expanded. Are there mask value for each point in the sample
    bool wantrow=true;			// do we output this row?
    if (hasmask)
    {
	maskbuffer.reset(mask.allocSampleBuffer());
	if (mask.actsExpanded())
	{
		maskoffset=DataTypes::noValues(mask.getDataPointShape());
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
	  	samples[d]=data[d].getSampleDataRO(i,bg[d].get());
	}
	if (hasmask)
	{
		masksample=mask.getSampleDataRO(i, maskbuffer.get());
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
#ifndef PASO_MPI
	throw;
#endif
    }
#ifdef PASO_MPI
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
#ifndef PASO_MPI

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
    const char* mpistr=0;
    MPI_File mpi_fileHandle_p;
    MPI_Status mpi_status;
    MPI_Info mpi_info = MPI_INFO_NULL;
    char* fname_c=new char[filename.size()+1];
    strcpy(fname_c,filename.c_str());
    boost::scoped_ptr<char> fname_p(fname_c);
     
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
                    MPI_CHAR, MPI_CHAR, "native", mpi_info);
// here we are assuming that std::string holds the same type of char as MPI_CHAR
    }

    std::string contents=os.str();
    char* con=new char[contents.size()+1];
    strcpy(con, contents.c_str());
    boost::scoped_ptr<char> buff(con);
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


}  // end of namespace
