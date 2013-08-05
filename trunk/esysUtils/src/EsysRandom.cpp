/*****************************************************************************
*
* Copyright (c) 2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "Esys_MPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

namespace {
    
boost::mt19937 base;		// used to seed all the other generators  
vector<boost::mt19937*> gens;
vector<uint32_t> seeds;

void seedGens(long seed)
{
#ifdef _OPENMP
    int numthreads=omp_get_max_threads();
#else
    int numthreads=1;
#endif
    if (gens.size()==0)		// we haven't instantiated the generators yet  
    {
        gens.resize(numthreads);	
        seeds.resize(numthreads);
    }  
    if (seed!=0)
    {
       base.seed((uint32_t)seed);	// without this cast, icc gets confused
       for (int i=0;i<numthreads;++i)
       {
	    uint32_t b=base();
            seeds[i]=b;	// initialise each generator with successive random values      
       }
       #pragma omp parallel for private(i)
       for (int i=0;i<numthreads;++i) 
       {
	   gens[i]=new boost::mt19937(seeds[i]);
       }
    }
}
  
  
}

// Put n random values into array
// Idea here is to create an array of seeds by feeding the original seed into the random generator
// The code at the beginning of the function to compute the seed if one is given is
// just supposed to introduce some variety (and ensure that multiple ranks don't get the same seed).
// I make no claim about how well these initial seeds are distributed
// uses openmp
// don't forget to call CHECK_FOR_EX_WRITE if using this on Data
void randomFillArray(long seed, double* array, size_t n)
{
    static unsigned prevseed=0;	// So if we create a bunch of objects we don't get the same start seed 
    if (seed==0)		// for each one
    {
	if (prevseed==0) 
	{
	    time_t s=time(0);
	    seed=s;
	}
	else
	{
	    seed=prevseed+419;	// these numbers are arbitrary
	    if (seed>3040101)		// I want to avoid overflow on 32bit systems
	    {
		seed=((int)(seed)%0xABCD)+1;
	    }
	}
    }  
    // now we need to consider MPI since we don't want each rank to start with the same seed. Rank in COMM_WORLD will do
#ifdef ESYS_MPI
    Esys_MPI_rank rank;
    int mperr=MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (mperr!=MPI_SUCCESS) {
        rank=0;
    }
    seed+=rank*5;
#endif
    prevseed=seed;  
    
    boost::mt19937::result_type RMAX=base.max();
    seedGens(seed);
    long i;
    
    #pragma omp parallel private(i)
    {
	int tnum=0;
	#ifdef _OPENMP
	tnum=omp_get_thread_num();
	#endif
	boost::mt19937& generator=*(gens[tnum]);
	
    	#pragma omp for schedule(static)
    	for (i=0;i<n;++i)
    	{
#ifdef _WIN32
	    array[i]=((double)generator())/RMAX;
#else
	    array[i]=((double)generator())/RMAX;
#endif
    	}
    }
}
