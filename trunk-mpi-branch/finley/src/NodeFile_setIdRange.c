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

/*   Finley: Mesh: NodeFile */

/*   returns the maximum and minimum node id number of nodes: */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"
#include "Util.h"

/**************************************************************/


void Finley_NodeFile_setGlobalIdRange(index_t* min_id,index_t* max_id,Finley_NodeFile* in) {
   index_t min_id_local, max_id_local;
   #ifdef PASO_MPI
   index_t global_id_range[2], id_range[2];
   #endif

   Finley_NodeFile_setIdRange(&min_id_local, &max_id_local,in);

   #ifdef PASO_MPI
   id_range[0]=-min_id_local;
   id_range[1]=max_id_local;
   MPI_Allreduce( id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm );
   *min_id=-global_id_range[0];
   *max_id=global_id_range[1];
   #else
   *min_id=min_id_local;
   *max_id=max_id_local;
   #endif
}

void Finley_NodeFile_setIdRange(index_t* min_id,index_t* max_id,Finley_NodeFile* in) {
   *min_id=Finley_Util_getMinInt(1,in->numNodes,in->Id);
   *max_id=Finley_Util_getMaxInt(1,in->numNodes,in->Id);
}
void Finley_NodeFile_setGlobalDOFRange(index_t* min_id,index_t* max_id,Finley_NodeFile* in) {
   index_t min_id_local, max_id_local;
   #ifdef PASO_MPI
   index_t global_id_range[2], id_range[2];
   #endif

   Finley_NodeFile_setDOFRange(&min_id_local, &max_id_local,in);

   #ifdef PASO_MPI
   id_range[0]=-min_id_local;
   id_range[1]=max_id_local;
   MPI_Allreduce( id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm );
   *min_id=-global_id_range[0];
   *max_id=global_id_range[1];
   #else
   *min_id=min_id_local;
   *max_id=max_id_local;
   #endif
}

void Finley_NodeFile_setDOFRange(index_t* min_id,index_t* max_id,Finley_NodeFile* in) {
   *min_id=Finley_Util_getMinInt(1,in->numNodes,in->globalDegreesOfFreedom);
   *max_id=Finley_Util_getMaxInt(1,in->numNodes,in->globalDegreesOfFreedom);
}

