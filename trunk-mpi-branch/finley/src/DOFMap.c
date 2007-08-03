/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "DOFMap.h"

Finley_DOFMap* Finley_DOFMap_alloc(dim_t numNodes, index_t* globalID, Paso_Distribution* distribution)
{
  Finley_DOFMap* out=NULL;
  /*  allocate the return value */
  out=MEMALLOC(1,Finley_DOFMap);
  if (Finley_checkPtr(out)) return NULL;

  out->distribution=Paso_Distribution_getReference(distribution);
  out->MPIInfo=Paso_MPIInfo_getReference(distribution->mpi_info);
  out->numNodes=0;
  out->ID=NULL;
  out->myNumDOFs=0; 
  out->numDOFs=0; 
  out->numRemotes=0; 
  out->remoteID=NULL; 
  out->remoteProcessor=NULL; 
  out->offsetInRemoteID=NULL; 
  out->numNeighbours=0;   
  out->neighbours=NULL;  
  out->reference_counter=1;
}

void Finley_DOFMap_free(Finley_DOFMap* in) {
  if ( in ) {
      in->reference_counter--;
      if (in->reference_counter<=0) {
         Paso_Distribution_free(in->distribution);
         Paso_MPIInfo_free(in->MPIInfo);
         MEMFREE(in->ID);
         MEMFREE(in->remoteID); 
         MEMFREE(in->remoteProcessor); 
         MEMFREE(in->offsetInRemoteID); 
         MEMFREE(in->neighbours);  
         MEMFREE(in);  
     }
  }
}
Finley_DOFMap* DOFMap_getReference(Finley_DOFMap *in ) 
{
  if ( in )
    in->reference_counter++;
  return in;
}
