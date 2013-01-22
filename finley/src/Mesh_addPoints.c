/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

/************************************************************************************

   Finley: adds points to Points table of the mesh

************************************************************************************/

#include "Mesh.h"

#ifdef ESYS_MPI    
void Finley_Mesh_MPI_minimizeDistance( void *, void*, int *, MPI_Datatype * );
 
void Finley_Mesh_MPI_minimizeDistance(void *invec_p, void *inoutvec_p, int *len, MPI_Datatype *dtype)
{
    const dim_t numPoints = (*len)/2;
    double *invec = (double*) invec_p;
    double *inoutvec = (double*)  inoutvec_p;
    int i;
    for ( i=0; i<numPoints; i++ ) {
        if ( invec[2*i] < inoutvec[2*i] ) { 
	     inoutvec[2*i]=invec[2*i];
	     inoutvec[2*i+1]=invec[2*i+1];
	} 
    }
}
#endif

void Finley_Mesh_addPoints(Finley_Mesh* mesh, const dim_t numPoints, const double *points_ptr, const index_t *tags_ptr)
{
  dim_t i,n,k, numOldPoints, numNewPoints;
  double *dist_p=NULL;
  index_t *node_id_p=NULL, *point_index_p=NULL;
  index_t order=mesh->integrationOrder;
  index_t reduced_order=mesh->reducedIntegrationOrder;
  const dim_t numDim=mesh->Nodes->numDim;
  Esys_MPIInfo * mpi_info= Esys_MPIInfo_getReference(mesh->MPIInfo);
  Finley_ElementFile *oldPoints=mesh->Points, *newPoints=NULL;
  Finley_ReferenceElementSet *refPoints=NULL;
  const index_t firstDOF=Paso_Distribution_getFirstComponent(mesh->Nodes->degreesOfFreedomDistribution);
  const index_t lastDOF=Paso_Distribution_getLastComponent(mesh->Nodes->degreesOfFreedomDistribution);
  
  
  if ( oldPoints == NULL) {
    refPoints=Finley_ReferenceElementSet_alloc(Finley_Point1, order, reduced_order);
    numOldPoints=0;
  } else {
    refPoints=Finley_ReferenceElementSet_reference(oldPoints->referenceElementSet);
    numOldPoints=mesh->Points->numElements;
  }  
  newPoints=Finley_ElementFile_alloc(refPoints, mpi_info);
  /* first we find the node which is the closest on this processor: */
  dist_p=TMPMEMALLOC(numPoints, double);
  point_index_p=TMPMEMALLOC(numPoints,index_t);
  node_id_p=TMPMEMALLOC(numPoints,index_t);
  
  for (i=0; i< numPoints; ++i) {
     dist_p[i]=LARGE_POSITIVE_FLOAT;
     node_id_p[i]=-1;
     node_id_p[i]=-1;
  }
  if ( numDim == 3 ) {
      #pragma omp parallel private(i)
      {
	  for (i=0; i< numPoints; ++i) {
	    register const double X0=points_ptr[INDEX2(0,i,numDim)]; 
	    register const double X1=points_ptr[INDEX2(1,i,numDim)];
	    register const double X2=points_ptr[INDEX2(2,i,numDim)];
	    register double dist_local=LARGE_POSITIVE_FLOAT;
	    register index_t node_id_local=-1;
	    #pragma omp for private(n)
	    for (n=0; n< mesh-> Nodes->numNodes; n++) {
		register const double D0=mesh-> Nodes->Coordinates[INDEX2(0,n,numDim)] - X0;
		register const double D1=mesh-> Nodes->Coordinates[INDEX2(1,n,numDim)] - X1;
		register const double D2=mesh-> Nodes->Coordinates[INDEX2(2,n,numDim)] - X2;
		register const double d= D0 * D0 + D1 * D1 + D2 * D2;
		if ( d < dist_local) {
		    dist_local = d;
		    node_id_local=n;
		}
	    }
	    #pragma omp critical
	    {
		if ( dist_local < dist_p[i] ) {
		    dist_p[i] = dist_local;
		    node_id_p[i] = node_id_local;
		}
	    }
	  }
      }
  } else if ( numDim == 2 ) {
	#pragma omp parallel private(i)
	{
	    for (i=0; i< numPoints; ++i) {
	      register const double X0=points_ptr[INDEX2(0,i,numDim)]; 
	      register const double X1=points_ptr[INDEX2(1,i,numDim)];
	      register double dist_local=LARGE_POSITIVE_FLOAT;
	      register index_t node_id_local=-1;
	      #pragma omp for private(n)
	      for (n=0; n< mesh-> Nodes->numNodes; n++) {
		  register const double D0=mesh-> Nodes->Coordinates[INDEX2(0,n,numDim)] - X0;
		  register const double D1=mesh-> Nodes->Coordinates[INDEX2(1,n,numDim)] - X1;
		  register const double d= D0 * D0 + D1 * D1;
		  if ( d < dist_local) {
		      dist_local = d;
		      node_id_local=n;
		  }
	      }
	      #pragma omp critical
	      {
		  if ( dist_local < dist_p[i] ) {
		      dist_p[i] = dist_local;
		      node_id_p[i] = node_id_local;
		  }
	      }
	    }
	}
  }  else {
	#pragma omp parallel private(i)
	{
	    for (i=0; i< numPoints; ++i) {
	      register const double X0=points_ptr[INDEX2(0,i,numDim)]; 
	      register double dist_local=LARGE_POSITIVE_FLOAT;
	      register index_t node_id_local=-1;
	      #pragma omp for private(n)
	      for (n=0; n< mesh-> Nodes->numNodes; n++) {
		  register const double D0=mesh-> Nodes->Coordinates[INDEX2(0,n,numDim)] - X0;
		  register const double d= D0 * D0;
		  if ( d < dist_local) {
		      dist_local = d;
		      node_id_local=n;
		  }
	      }
	      #pragma omp critical
	      {
		  if ( dist_local < dist_p[i] ) {
		      dist_p[i] = dist_local;
		      node_id_p[i] = node_id_local;
		  }
	      }
	    }
	}
  }

  /* now we need to reduce this across all processors */
  #ifdef ESYS_MPI
    { 
         MPI_Op op;
         int count = 2* numPoints;
         double *sendbuf =NULL, *recvbuf=NULL;
	 sendbuf=TMPMEMALLOC(count,double);
	 recvbuf=TMPMEMALLOC(count,double);
	 
	 for (i=0; i< numPoints; ++i) {
              sendbuf[2*i  ]= dist_p[i];
	      sendbuf[2*i+1]= (double) (mesh-> Nodes->Id[node_id_p[i]]);
         }
	 MPI_Op_create(Finley_Mesh_MPI_minimizeDistance, TRUE,&op);
         MPI_Allreduce ( sendbuf,  recvbuf, count, MPI_DOUBLE, op, mpi_info->comm);
	 MPI_Op_free(&op);
	 /* if the node id has changed we found another node which is closer elsewhere */
         for (i=0; i< numPoints; ++i) {
	      register const int best_fit_Id = (int) (recvbuf[2*i+1] + 0.5);
	      if ( best_fit_Id != mesh-> Nodes->Id[node_id_p[i]] ) {
		    node_id_p[i] =-1;
	      }
         }
         TMPMEMFREE(sendbuf);
	 TMPMEMFREE(recvbuf);
    }
  #endif
  /*  we pick the points to be used on this processor */
  numNewPoints=0;
  for (i=0; i< numPoints; ++i) {
     if (node_id_p[i]>-1) { /* this processor uses a node which is identical to point i */

        if  ( mesh-> Nodes->globalReducedDOFIndex[node_id_p[i]]  >-1) { /* the point is also used in the reduced mesh */

	  register const index_t global_id=mesh->Nodes->globalDegreesOfFreedom[node_id_p[i]];
	  if ( (firstDOF<= global_id) && ( global_id <lastDOF) ) {  /* is this point actually relevant */


	     /* is this point already in the Point table? */
	       bool_t notAnOldPoint=TRUE;

	       if (numOldPoints >0) {
		    for (k=0; k<numOldPoints; ++k) {
		        if ( global_id == oldPoints->Nodes[k]) {
			  notAnOldPoint=FALSE;
			  break;
			}
		    }
	       }
	       if (notAnOldPoint) {
		    /* is this point unique in the new list of points? */
		    bool_t notANewPoint=TRUE;
		    for (k=0; k<numNewPoints; ++k) {
		        if ( global_id == mesh->Nodes->globalDegreesOfFreedom[node_id_p[point_index_p[k]]]) {
			  notANewPoint=FALSE;
			  break;
			}
		    }
		    if (notANewPoint) {
		      point_index_p[numNewPoints]=i;
		      numNewPoints++;
		   }
	       } 
	   }
	}
     }
   }
   /* now we are ready to create the new Point table */
   Finley_ElementFile_allocTable(newPoints,numOldPoints+numNewPoints);
   if (numOldPoints > 0) {
       #pragma omp parallel for private(n) schedule(static)
       for(n=0;n<numOldPoints;n++) {
             newPoints->Owner[n]=oldPoints->Owner[n];
             newPoints->Id[n]   =oldPoints->Id[n];
             newPoints->Tag[n]  =oldPoints->Tag[n];
             newPoints->Nodes[n]=oldPoints->Nodes[n];
	     newPoints->Color[n]=0;
       }
   }
   #pragma omp parallel for private(n) schedule(static)
   for(n=0;n<numNewPoints;n++) {
             register const index_t idx = point_index_p[n];
             newPoints->Owner[numOldPoints+n]=mpi_info->rank;
             newPoints->Id[numOldPoints+n]   =0;
             newPoints->Tag[numOldPoints+n]  =tags_ptr[idx];
             newPoints->Nodes[numOldPoints+n]=node_id_p[idx];
	     newPoints->Color[numOldPoints+n]=0;
   }
   newPoints->minColor=0;
   newPoints->maxColor=0;
  
  /* all done */
  TMPMEMFREE(dist_p);
  TMPMEMFREE(node_id_p);
  TMPMEMFREE(point_index_p);
  Finley_ReferenceElementSet_dealloc(refPoints);
  Esys_MPIInfo_free(mpi_info);
  if (Finley_noError()) {
      Finley_ElementFile_free(oldPoints);
      mesh->Points=newPoints;
  } else {
     Finley_ElementFile_free(newPoints);
  }
}
