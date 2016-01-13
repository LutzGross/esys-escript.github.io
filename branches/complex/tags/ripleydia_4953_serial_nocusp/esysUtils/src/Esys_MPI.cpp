
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "Esys_MPI.h"
#include "index.h"
#include "mem.h"
#include "error.h"
#include "EsysException.h"


#include <iostream>	// temp for debugging

namespace esysUtils
{
  
JMPI makeInfo(MPI_Comm comm, bool owncom)
{
    if (esysUtils::NoCOMM_WORLD::active() && comm==MPI_COMM_WORLD)
    {
	throw esysUtils::EsysException("Attempt to use the MPI_COMM_WORLD communicator when it is blocked.");
    }
    JMPI_* p=new JMPI_(comm, owncom);
    return JMPI(p);
}


JMPI_::JMPI_(MPI_Comm mpicomm, bool owncom)
	: comm(mpicomm), ownscomm(owncom)
{
	msg_tag_counter = 0;
#ifdef ESYS_MPI
	if (MPI_Comm_rank(comm, &rank)!=MPI_SUCCESS || MPI_Comm_size(comm, &size)!=MPI_SUCCESS)
	{
	    Esys_setError( ESYS_MPI_ERROR, "Esys_MPIInfo_alloc : error finding comm rank/size" );
	}
#else
	rank=0;
	size=1;	
#endif	
}

JMPI_::~JMPI_()
{
#ifdef ESYS_MPI
    if (ownscomm)
    {
	MPI_Comm_free(&comm);
    }
#endif
}

dim_t JMPI_::setDistribution(index_t min_id,index_t max_id,index_t* distribution)
{
   int rest=0, p;
   int s=size;
   dim_t N=max_id-min_id+1;
   if (N>0) {
      int local_N=N/s;
      rest=N-local_N*s;
      for (p=0; p<s; ++p) {
         if (p<rest) {
             distribution[p]=min_id+(local_N+1)*p;
         } else {
             distribution[p]=min_id+rest+local_N*p;
         }
      }
      distribution[s]=max_id+1;
      if (rest==0) {
         return local_N;
      } else {
         return local_N+1;
      }
  } else {
      for (p=0; p<s+1; ++p) distribution[p]=min_id;
      return 0;
  }  
  
  
}

void JMPI_::split(dim_t N, dim_t* local_N,index_t* offset) 
{
   int rest=0;
   int s=size;
   int r=rank;
   *local_N=N/s;
   rest=N-(*local_N)*s;
   if (r<rest) {
       (*local_N)++;
       (*offset)=(*local_N)*r;
   } else {
       (*offset)=(*local_N)*r+rest;
   }
}

}



dim_t Esys_MPIInfo_setDistribution(esysUtils::JMPI& mpi_info ,index_t min_id,index_t max_id,index_t* distribution) {
   int rest=0, p;
   int s=mpi_info->size;
   dim_t N=max_id-min_id+1;
   if (N>0) {
      int local_N=N/s;
      rest=N-local_N*s;
      for (p=0; p<s; ++p) {
         if (p<rest) {
             distribution[p]=min_id+(local_N+1)*p;
         } else {
             distribution[p]=min_id+rest+local_N*p;
         }
      }
      distribution[s]=max_id+1;
      if (rest==0) {
         return local_N;
      } else {
         return local_N+1;
      }
  } else {
      for (p=0; p<s+1; ++p) distribution[p]=min_id;
      return 0;
  }
}



/* N = #CPUs, k is a CPU number but out of range or even negative. Return a CPU number in 0...n-1. */
index_t esysUtils::mod_rank(index_t n, index_t k) 
{
    index_t q, out=0;
    if (n>1) {
        q=k/n;
        if (k>0) {
           out=k-n*q;
        } else if (k<0) {
           out=k-n*(q-1);
        }
    }
    return out;
}


/* checks that there is no error across all processes in a communicator */
/* NOTE : does not make guarantee consistency of error string on each process */
bool esysUtils::Esys_MPIInfo_noError( const esysUtils::JMPI& mpi_info )
{
  int errorLocal = Esys_noError() ? 0 : 1;
  int errorGlobal = errorLocal;

#ifdef ESYS_MPI
  if (!checkResult(errorLocal, errorGlobal, mpi_info->comm))
  {
      return false;
  }
  if( (errorLocal==0) && (errorGlobal==1)) 
  {
     Esys_setError( ESYS_MPI_ERROR, "Esys_MPIInfo_noError() : there was an error on another MPI process" );
  }
#endif
  
  return (errorGlobal==0);
}

/* returns the max of inputs on all ranks -- or just sends the input back on nompi */
bool esysUtils::checkResult(int& input, int& output, MPI_Comm& comm)
{
#ifdef ESYS_MPI
    output=0;
    if (MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MAX, comm)!=MPI_SUCCESS)
    {
	return false;
    }
    return true;
#else
    output=input;
    return true;
#endif
}




// ensure that the any ranks with an empty src argument end up with the string from
// one of the other ranks
// with no-mpi, it makes dest point at a copy of src
// Expected use case for this code is to ship error messages between ranks
// as such, it is not written to be speedy
bool esysUtils::shipString(const char* src, char** dest, MPI_Comm& comm)
{
#ifdef ESYS_MPI  
    Esys_MPI_rank rank=0;
    if (MPI_Comm_rank( comm, &rank )!=MPI_SUCCESS)
    {
	return false;	// we have no reason to believe MPI works anymore
    }
    
    int slen=strlen(src);
    // everybody needs to tell everyone if they have a string
    // send your rank if you have a non-empty string else
    // send -1
    int in=(slen?rank:-1);
    int out;
    if (MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MAX, comm)!=MPI_SUCCESS)
    {
	return false;
    }
    if (out==-1)		// should not be called under these conditions, but noone had a string
    {
	*dest=new char[1];
	*dest[0]='\0';
	return true;
    }
    // since we will be using broadcast, we need to tell everyone how big the string is going to be
    // with an additional bcast
    
    if (MPI_Bcast(&slen, 1, MPI_INT, out, comm)!=MPI_SUCCESS)
    {
	return false;
    }
    // now broadcast that string to everyone
    if (rank==out)
    {
	// I could const _ cast src but instead I'll make a copy
	
	*dest=new char[slen+1];
	strcpy(*dest, src);
	
	// this guy should just send the string
	if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm)!=MPI_SUCCESS)
	{
	    return false;
	}
	return true;
    }
    else
    {
	*dest=new char[slen+1];
	if (MPI_Bcast(*dest, slen+1, MPI_CHAR, out, comm)!=MPI_SUCCESS)
	{
	    return false;
	}
	return true;
    }
#else
    *dest=new char[strlen(src)+1];
    strcpy(*dest, src);
    return true;
#endif
  
}

namespace 
{
    // true if a split world call is currently running and MPI_COMM_WORLD should not be allowed by default
    bool nocommworldplease=false;
}

esysUtils::NoCOMM_WORLD::NoCOMM_WORLD()
{
    if (nocommworldplease)
    {
	throw EsysException("NoCOMM_WORLD does not nest.");
    }
    nocommworldplease=true;
}

esysUtils::NoCOMM_WORLD::~NoCOMM_WORLD()
{
    nocommworldplease=false;
}  

bool esysUtils::NoCOMM_WORLD::active()
{
    return nocommworldplease;
}

/**************************************************
                 WRAPPERS 
**************************************************/

int Esys_MPIInfo_initialized( void )
{
  #ifdef ESYS_MPI
     int error=0, initialised=0;
     error = MPI_Initialized( &initialised );
     if( error!=MPI_SUCCESS )
         Esys_setError( ESYS_MPI_ERROR, "mpi_initialised : MPI error" );
     return initialised;
  #else
     return TRUE;
  #endif
}

#ifndef _OPENMP 
int serial_get_max_threads(void) {
   return 1;
}
int serial_get_thread_num(void) {
   return 0;
}
#endif

