/*****************************************************************************
*
* Copyright (c) 2013-2014 by University of Queensland
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

#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <cstring>
#include "Esys_MPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

namespace {
    
boost::mt19937 base;		// used to seed all the other generators  
vector<boost::mt19937*> gens;
vector<boost::uint32_t> seeds;

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
       int i;
       base.seed((boost::uint32_t)seed);	// without this cast, icc gets confused
       for (int i=0;i<numthreads;++i)
       {
	    boost::uint32_t b=base();
            seeds[i]=b;	// initialise each generator with successive random values      
       }
       #pragma omp parallel for private(i)
       for (i=0;i<numthreads;++i) 
       {
	   gens[i]=new boost::mt19937(seeds[i]);
       }
    }
}
  
  
}

namespace esysUtils
{

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

// see patternFillArray for details on parameters
void patternFillArray2D(size_t x, size_t y, double* array, size_t spacing, size_t basex, size_t basey, size_t numpoints)
{
      memset(array, 0, x*y*sizeof(double)*numpoints);
      size_t xoff=basex%spacing;
      size_t yoff=basey%spacing;
      for (int r=0;r<y;++r)
      {
	  size_t step=((r+yoff)%spacing)?spacing:1; 
	  for (int c=0;c<x;++c)
	  {
	      if ((c+xoff)%step==0)
	      {
		  for (int p=0;p<numpoints;++p)
		  {
		      array[(c+r*x)*numpoints+p]=1+p;
		  }		
	      }	    
	  }
      } 
}


// fill the array (which we assume is 3D with x by y by z points in it) with a pattern.
// The base? params give the coordinates (in # of elements) of the origin of _this_ rank
//  used to ensure patterns are generated consistantly across multiple ranks
// This is only for internal debug so the patterns (or this function) may disappear 
// without notice
void patternFillArray(int pattern, size_t x, size_t y, size_t z, double* array, size_t spacing, size_t basex, size_t basey, size_t basez, size_t numpoints)
{
    if (pattern==0)	// a cross pattern in the z=0 plane, repeated for each z layer
    {
	memset(array, 0, x*y*sizeof(double)*numpoints);
	size_t xoff=basex%spacing;
	size_t yoff=basey%spacing;
	for (int r=0;r<y;++r)
	{
	    size_t step=((r+yoff)%spacing)?spacing:1;
	    for (int c=0;c<x;++c)
	    {
		if ((c+xoff)%step==0)
		{
		    for (int p=0;p<numpoints;++p)
		    {
			array[(c+r*x)*numpoints+p]=p+1;
		    }
		}
	    }
	}
	for (int l=1;l<z;++l)
	{
	    memcpy(array+(x*y*l*numpoints), array, x*y*sizeof(double)*numpoints);
	}
    }
    else		// pattern 1. A grid in all 3 dimensions 
    {
	if (z<2)
	{
	    patternFillArray(0, x, y, z, array, spacing, basex, basey, basez, numpoints);
	    return;	// this pattern needs a minimum of 2 layers
	}
	size_t xoff=basex%spacing;
	size_t yoff=basey%spacing;
	size_t zoff=basez%spacing;
	
	double* buff1=new double[x*y*numpoints];	// stores the main cross pattern
	double* buff2=new double[x*y*numpoints];	// stores the "verticals"
	memset(buff1, 0, x*y*sizeof(double)*numpoints);
	memset(buff2, 0, x*y*sizeof(double)*numpoints);
	    // fill in buff1
	for (size_t r=0;r<y;++r)
	{
	    size_t step=((r+yoff)%spacing)?spacing:1;
	    for (int c=0;c<x;++c)
	    {
		if ((c+xoff)%step==0)
		{
		    for (int p=0;p<numpoints;++p)
		    {
			buff1[(c+r*x)*numpoints+p]=p+1;
		    }
		}
	    }	    
	}
	
	for (size_t r=(spacing-yoff)%spacing;r<y;r+=spacing)
	{
	    for (size_t c=(spacing-xoff)%spacing;c<x;c+=spacing)
	    {
		for (int p=0;p<numpoints;++p)
		{
		    buff2[(c+r*x)*numpoints+p]=p+1;
		}
	    }
	}	
	for (size_t l=0;l<z;++l)
	{
	    if ((l+zoff)%spacing)
	    {
		memcpy(array+(x*y*l*numpoints), buff2, x*y*sizeof(double)*numpoints);
	    }
	    else
	    {
		memcpy(array+(x*y*l*numpoints), buff1, x*y*sizeof(double)*numpoints);
	    }
	}
	delete[] buff1;
	delete[] buff2;
    }
  
  
  
}

} // end namespace
