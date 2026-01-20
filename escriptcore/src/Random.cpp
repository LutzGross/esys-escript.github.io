/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifdef _WIN32
#include "Random.h"
#endif

#include <escript/EsysMPI.h>

#include <algorithm>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

using namespace std;

namespace {
    
boost::mt19937 base;            // used to seed all the other generators  
vector<boost::mt19937*> gens;
vector<boost::uint32_t> seeds;

void seedGens(long seed)
{
#ifdef _OPENMP
    int numthreads=omp_get_max_threads();
#else
    int numthreads=1;
#endif
    if (gens.size()==0)         // we haven't instantiated the generators yet  
    {
        gens.resize(numthreads);        
        seeds.resize(numthreads);
    }  
    if (seed!=0)
    {
        base.seed((boost::uint32_t)seed); // without this cast, icc gets confused
        for (int i=0;i<numthreads;++i)
        {
            boost::uint32_t b=base();
            seeds[i]=b; // initialise each generator with successive random values      
        }
#pragma omp parallel for
        for (int i=0; i<numthreads; ++i)
        {
            gens[i]=new boost::mt19937(seeds[i]);
        }
    }
}
  
} // anonymous namespace

namespace escript
{

// Put n random values from the interval [0,1] into array.
// Idea here is to create an array of seeds by feeding the original seed into
// the random generator.
// The code at the beginning of the function to compute the seed if one is
// given is just supposed to introduce some variety (and ensure that multiple
// ranks don't get the same seed).
// I make no claim about how well these initial seeds are distributed.
// notes:
// - uses openmp
// - don't forget to call CHECK_FOR_EX_WRITE if using this on Data
void randomFillArray(long seed, double* array, size_t n, JMPI mpiInfo)
{
    // So if we create a bunch of objects we don't get the same start seed
    static unsigned prevseed=0;
    if (seed==0)                // for each one
    {
        if (prevseed==0)
        {
            time_t s=time(0);
            seed=s;
        }
        else
        {
            seed=prevseed+419;  // these numbers are arbitrary
            if (seed>3040101)   // I want to avoid overflow on 32bit systems
            {
                seed=((int)(seed)%0xABCD)+1;
            }
        }
    }
    // now we need to consider MPI since we don't want each rank to start with
    // the same seed. Use rank from the provided communicator
#ifdef ESYS_MPI
    seed += mpiInfo->rank * 5;
#endif
    prevseed=seed;  
    
    boost::mt19937::result_type RMAX=base.max();
    seedGens(seed);
    
#pragma omp parallel
    {
#ifdef _WIN32 // error C3016: 'i': index variable in OpenMP 'for' statement must have signed integral type
#define OMP_LOOP_IDX_T int
#else
#define OMP_LOOP_IDX_T size_t
#endif
        OMP_LOOP_IDX_T i;
#ifdef _OPENMP
        int tnum=omp_get_thread_num();
#else
        int tnum=0;
#endif
        boost::mt19937& generator=*(gens[tnum]);
        
#pragma omp for schedule(static)
        for (i=0;i<n;++i)
        {
            array[i]=((double)generator())/RMAX;
        }
    }
}

// see patternFillArray for details on parameters
void patternFillArray2D(size_t x, size_t y, double* array, size_t spacing,
                        size_t basex, size_t basey, size_t numpoints)
{
    std::fill_n(array, x*y*numpoints, 0.);
    size_t xoff=basex%spacing;
    size_t yoff=basey%spacing;
    for (size_t r=0;r<y;++r)
    {
        size_t step=((r+yoff)%spacing)?spacing:1; 
        for (size_t c=0;c<x;++c)
        {
            if ((c+xoff)%step==0)
            {
                for (size_t p=0;p<numpoints;++p)
                {
                    array[(c+r*x)*numpoints+p]=1+p;
                }             
            }
        }
    }
}


// fill the array (which we assume is 3D with x by y by z points in it) with
// a pattern.
// The base? params give the coordinates (in # of elements) of the origin of
// _this_ rank.
// Used to ensure patterns are generated consistently across multiple ranks.
// This is only for internal debug so the patterns (or this function) may
// disappear without notice
void patternFillArray(int pattern, size_t x, size_t y, size_t z, double* array,
                      size_t spacing, size_t basex, size_t basey, size_t basez,
                      size_t numpoints)
{
    if (pattern==0)     // a cross pattern in the z=0 plane, repeated for each z layer
    {
        std::fill_n(array, x*y*numpoints, 0.);
        size_t xoff=basex%spacing;
        size_t yoff=basey%spacing;
        for (size_t r=0;r<y;++r)
        {
            size_t step=((r+yoff)%spacing)?spacing:1;
            for (size_t c=0;c<x;++c)
            {
                if ((c+xoff)%step==0)
                {
                    for (size_t p=0;p<numpoints;++p)
                    {
                        array[(c+r*x)*numpoints+p]=p+1;
                    }
                }
            }
        }
        for (size_t l=1;l<z;++l)
        {
            std::copy(array, &array[x*y*numpoints], &array[l*x*y*numpoints]);
        }
    }
    else                // pattern 1. A grid in all 3 dimensions 
    {
        if (z<2)
        {
            patternFillArray(0, x, y, z, array, spacing, basex, basey, basez, numpoints);
            return;     // this pattern needs a minimum of 2 layers
        }
        size_t xoff=basex%spacing;
        size_t yoff=basey%spacing;
        size_t zoff=basez%spacing;
        
        double* buff1=new double[x*y*numpoints];        // stores the main cross pattern
        double* buff2=new double[x*y*numpoints];        // stores the "verticals"
        std::fill_n(buff1, x*y*numpoints, 0.);
        std::fill_n(buff2, x*y*numpoints, 0.);
        // fill in buff1
        for (size_t r=0;r<y;++r)
        {
            size_t step=((r+yoff)%spacing)?spacing:1;
            for (size_t c=0;c<x;++c)
            {
                if ((c+xoff)%step==0)
                {
                    for (size_t p=0;p<numpoints;++p)
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
                for (size_t p=0;p<numpoints;++p)
                {
                    buff2[(c+r*x)*numpoints+p]=p+1;
                }
            }
        }       
        for (size_t l=0;l<z;++l)
        {
            if ((l+zoff)%spacing)
            {
                std::copy(buff2, buff2+x*y*numpoints, &array[x*y*l*numpoints]);
            }
            else
            {
                std::copy(buff1, buff1+x*y*numpoints, &array[x*y*l*numpoints]);
            }
        }
        delete[] buff1;
        delete[] buff2;
    }
}

} // end namespace

