
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


#include <string.h>

// added for saveCSV
#include <boost/python.hpp>
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

ESCRIPT_DLL_API
void
saveDataCSV(const std::string& filename, boost::python::dict arg, const std::string& sep, const std::string& csep, 
bool append)
{
    using std::cout;
    using std::endl;
    boost::python::list keys=arg.keys();
    int numdata = boost::python::extract<int>(arg.attr("__len__")());
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
	step[i]=(data[i].actsConstant()?0:DataTypes::noValues(data[i].getDataPointShape()));
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

    
    std::ofstream os;
    if (append)
    {
	os.open(filename.c_str(), std::ios_base::app);
    }
    else
    {
	os.open(filename.c_str());
    }
    if (!os.is_open())
    {
	throw DataException("saveDataCSVcpp: unable to open file for writing");
    }

    bool first=true;
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

	//the use of shared_ptr here is just to ensure the buffer group is freed
	//I would have used scoped_ptr but they don't work in vectors
    std::vector<boost::shared_ptr<BufferGroup> > bg(numdata);
    for (int d=0;d<numdata;++d)
    {
	bg[d].reset(data[d].allocSampleBuffer());
    }

    try{
      for (int i=0;i<numsamples;++i)
      {
	for (int d=0;d<numdata;++d)
	{
	  	samples[d]=data[d].getSampleDataRO(i,bg[d].get());
	}
	for (int j=0;j<dpps;++j)
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
	for (int d=0;d<numdata;++d)
	{
	    offset[d]=0;
	}	
      }
    } catch (...)
    {
        os.close();
	throw;
    }
    os.close();
}

}  // end of namespace
