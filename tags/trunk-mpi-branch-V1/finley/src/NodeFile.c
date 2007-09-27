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
/*                                                             */
/*   Finley: Mesh : NodeFile */
/*                                                             */
/*   allocates and deallocates node files                      */
/*                                                             */
/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

/*   allocates a node file to hold nodes */
/*   use Finley_NodeFile_allocTable to allocate the node table (Id,Coordinatess). */
#ifdef PASO_MPI
Finley_NodeFile* Finley_NodeFile_alloc(dim_t numDim, Paso_MPIInfo *MPIInfo){
#else
Finley_NodeFile* Finley_NodeFile_alloc(dim_t numDim){
#endif
  Finley_NodeFile *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_NodeFile);
  if (Finley_checkPtr(out)) return NULL;
  out->isPrepared=FINLEY_UNKNOWN;
  out->numNodes=0;
  out->numDegreesOfFreedom=0;
  out->reducedNumDegreesOfFreedom=0;
  out->reducedNumNodes=0;
  out->numDim=numDim;
  out->Id=NULL;
  out->Tag=NULL;
  out->Coordinates=NULL;
  out->degreeOfFreedom=NULL;
  out->degreeOfFreedomId=NULL;
  out->reducedDegreeOfFreedom=NULL;
  out->reducedDegreeOfFreedomId=NULL;
  out->toReduced=NULL;
  out->status=FINLEY_INITIAL_STATUS;
#ifdef PASO_MPI
  out->Dom=NULL;
  out->MPIInfo = Paso_MPIInfo_getReference( MPIInfo );
  out->degreeOfFreedomDistribution = Finley_NodeDistribution_alloc( MPIInfo );
  out->reducedDegreeOfFreedomDistribution = Finley_NodeDistribution_alloc( MPIInfo );
  out->CommBuffer = Paso_CommBuffer_alloc( MPIInfo, __g_nodeTag++ );
  out->reducedCommBuffer = Paso_CommBuffer_alloc( MPIInfo, __g_nodeTag++ );
#endif
  return out;
}

/*  deallocates a node file: */

void Finley_NodeFile_dealloc(Finley_NodeFile* in) {
  if (in!=NULL) {
     #ifdef Finley_TRACE
     printf("node file is deallocated.\n");
     #endif
     Finley_NodeFile_deallocTable(in);
#ifdef PASO_MPI
     Paso_MPIInfo_dealloc( in->MPIInfo );
     Finley_NodeDistribution_dealloc( in->degreeOfFreedomDistribution ); 
     Finley_NodeDistribution_dealloc( in->reducedDegreeOfFreedomDistribution );
     Paso_CommBuffer_dealloc( in->CommBuffer );
     Paso_CommBuffer_dealloc( in->reducedCommBuffer );
#endif
     MEMFREE(in);      
  }
}
