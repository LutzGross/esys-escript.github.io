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

/*   Finley: Mesh */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*   allocates a Mesh with name name for elements of type id using an integration order. If order is negative, */
/*   the most appropriate order is selected indepently. */

extern Finley_RefElementInfo Finley_RefElement_InfoList[];

#ifndef PASO_MPI
Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order) 
#else
Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order, Paso_MPIInfo *mpi_info) 
#endif
{
  Finley_Mesh *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_Mesh);
  if (Finley_checkPtr(out)) return NULL;
  out->Name=NULL;  
  out->Nodes=NULL; 
  out->Elements=NULL;   
  out->FaceElements=NULL; 
  out->Points=NULL;      
  out->ContactElements=NULL;      
  out->TagMap=NULL;      
  out->reference_counter=0;

  out->FullFullPattern=NULL;
  out->FullReducedPattern=NULL;
  out->ReducedFullPattern=NULL;
  out->ReducedReducedPattern=NULL;

#ifdef PASO_MPI 
  out->MPIInfo = NULL;
 
  /* get MPI info */
  out->MPIInfo = Paso_MPIInfo_getReference( mpi_info );
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
#endif
  /*   copy name: */
  
  out->Name=MEMALLOC(strlen(name)+1,char);
  if (Finley_checkPtr(out->Name)) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  strcpy(out->Name,name);
  
  /*   allocate node table: */
#ifdef PASO_MPI
  out->Nodes=Finley_NodeFile_alloc( numDim, mpi_info );
#else
  out->Nodes=Finley_NodeFile_alloc(numDim);
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  out->order=order;
  out->Elements=NULL;
  out->FaceElements=NULL;
  out->Points=NULL;
  out->ContactElements=NULL;
  out->reference_counter++;
  return out;
}

/* returns a reference to Finley_Mesh in */

Finley_Mesh* Finley_Mesh_reference(Finley_Mesh* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

/*   deallocates a mesh: */

void Finley_Mesh_dealloc(Finley_Mesh* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<1) {
       #ifdef Finley_TRACE
       if (in->Name!=NULL) {
           printf("Finley_Mesh_dealloc: mesh %s is deallocated.\n",in->Name);
       } else {
           printf("Finley_Mesh_dealloc\n");
       }
       #endif
       MEMFREE(in->Name);
       Finley_NodeFile_dealloc(in->Nodes);
       Finley_ElementFile_dealloc(in->Elements);   
       Finley_ElementFile_dealloc(in->FaceElements);
       Finley_ElementFile_dealloc(in->ContactElements);
       Finley_ElementFile_dealloc(in->Points);
       Finley_TagMap_free(in->TagMap);
       Paso_SystemMatrixPattern_dealloc(in->FullFullPattern);
       Paso_SystemMatrixPattern_dealloc(in->FullReducedPattern);
       Paso_SystemMatrixPattern_dealloc(in->ReducedFullPattern);
       Paso_SystemMatrixPattern_dealloc(in->ReducedReducedPattern);
#ifdef PASO_MPI
       Paso_MPIInfo_dealloc( in->MPIInfo );
#endif
       MEMFREE(in);      
     }
  }
}

/**************************************************************/

/*  returns the spatial dimension of the mesh: */

dim_t Finley_Mesh_getDim(Finley_Mesh *in) {
  return in->Nodes->numDim;
}

/**************************************************************/

/*  returns the number of nodes in the mesh: */

dim_t Finley_Mesh_getNumNodes(Finley_Mesh *in) {
  return in->Nodes->numNodes;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

dim_t Finley_Mesh_getNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->numDegreesOfFreedom;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

dim_t Finley_Mesh_getReducedNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->reducedNumDegreesOfFreedom;
}

#ifdef PASO_MPI
void print_mesh_statistics( Finley_Mesh *out, bool_t reduced )
{
  index_t i, r, j, N, M, dim;
	int *ref;  
	dim = out->Nodes->numDim;

	if( !reduced ){
		printf( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\nMESH STATISTICS\n\nFULL MESH\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
		printf( "\nNodes\n========\n\n" );
		printf( "\t%d internal DOF\n\t%d boundary DOF\n\t%d local DOF\n\t%d external DOF\n", out->Nodes->degreeOfFreedomDistribution->numInternal, out->Nodes->degreeOfFreedomDistribution->numBoundary, out->Nodes->degreeOfFreedomDistribution->numLocal, out->Nodes->degreeOfFreedomDistribution->numExternal);
		for( i=0; i<out->Nodes->numNodes; i++ ){
			printf( "node %d\t: id %d   \tDOF %d   \t: tag %d  \t: Dom %d  \t: coordinates [%3g", i, out->Nodes->Id[i], out->Nodes->degreeOfFreedom[i], out->Nodes->Tag[i], out->Nodes->Dom[i], out->Nodes->Coordinates[INDEX2(0,i,dim)] );
			for( j=1; j<dim; j++ )
				printf( ", %3g",  out->Nodes->Coordinates[INDEX2(j,i,dim)]);
			printf( " ]\n" );
		}

		printf( "Elements\n=========\n\n" );
		printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->Elements->elementDistribution->numInternal, out->Elements->elementDistribution->numBoundary, out->Elements->elementDistribution->numLocal );
		N = out->Elements->ReferenceElement->Type->numNodes;
		for( i=0; i<out->Elements->numElements; i++ ){
			printf( "element %d    \t: id %d  \t: dom %d  \t: nodes [ %2d", i, out->Elements->Id[i], out->Elements->Dom[i], out->Elements->Nodes[INDEX2(0,i,N)] );
			for( j=1; j<N; j++ )
				printf( ", %2d", out->Elements->Nodes[INDEX2(j,i,N)]);	
			printf( " ] -> " );	
			if( N>8 )
				printf( "\n\t\t\t\t\t\t" );
			printf( ": DOF   [ %2d", out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,i,N)]] );	
			for( j=1; j<N; j++ )
				printf( ", %2d", out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(j,i,N)]]);	
			printf( " ]\n" );	
		}

		printf( "\nFace Elements\n==========\n\n" );
		printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->FaceElements->elementDistribution->numInternal, out->FaceElements->elementDistribution->numBoundary, out->FaceElements->elementDistribution->numLocal );
		N = out->FaceElements->ReferenceElement->Type->numNodes;
		for( i=0; i<out->FaceElements->numElements; i++ ){
			printf( "face element %d \t: id %d  \t: dom %d  \t: nodes [ %2d", i, out->FaceElements->Id[i], out->FaceElements->Dom[i], out->FaceElements->Nodes[INDEX2(0,i,N)] );
			for( j=1; j<N; j++ )
				printf( ", %2d", out->FaceElements->Nodes[INDEX2(j,i,N)]  );	
			printf( " ] -> " );	
			if( N>8 )
				printf( "\n\t\t\t\t\t\t" );
			printf( ": DOF   [ %2d", out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,i,N)]]);	
			for( j=1; j<N; j++ )
				printf( ", %2d", out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(j,i,N)]]);	
			printf( " ]\n" );	
		}
		printf( "\nDistribution Data\n==========\n\n" );
		Finley_NodeDistribution_print( out->Nodes->degreeOfFreedomDistribution, stdout );
	}
	else{
		printf( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\nMESH STATISTICS\n\nREDUCED MESH\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
		printf( "\nNodes\n========\n\n" );
		printf( "\t%d internal DOF\n\t%d boundary DOF\n\t%d local DOF\n\t%d external DOF\n", out->Nodes->reducedDegreeOfFreedomDistribution->numInternal, out->Nodes->reducedDegreeOfFreedomDistribution->numBoundary, out->Nodes->reducedDegreeOfFreedomDistribution->numLocal, out->Nodes->reducedDegreeOfFreedomDistribution->numExternal);
		for( i=0, r=0; i<out->Nodes->numNodes; i++ ){
			if( out->Nodes->toReduced[i]>=0 ){
				printf( "node %d   \t: id %d   \tDOF %d   \t: tag %d  \t: Dom %d  \t: coordinates [%3g", r, out->Nodes->Id[i], out->Nodes->reducedDegreeOfFreedom[i], out->Nodes->Tag[i], out->Nodes->Dom[i], out->Nodes->Coordinates[INDEX2(0,i,dim)] );
				for( j=1; j<dim; j++ )
					printf( ", %3g",  out->Nodes->Coordinates[INDEX2(j,i,dim)]);
				printf( " ]\n" );
				r++;
			}
		}

		printf( "Elements\n=========\n\n" );
		printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->Elements->elementDistribution->numInternal, out->Elements->elementDistribution->numBoundary, out->Elements->elementDistribution->numLocal );
		N = out->Elements->LinearReferenceElement->Type->numNodes;
		M = out->Elements->ReferenceElement->Type->numNodes;
		ref = out->Elements->ReferenceElement->Type->linearNodes;
		for( i=0; i<out->Elements->numElements; i++ ){
			printf( "element %d    \t: id %d  \t: dom %d  \t: nodes [ %3d", i, out->Elements->Id[i], out->Elements->Dom[i], out->Nodes->toReduced[out->Elements->Nodes[INDEX2(ref[0],i,M)]] );
			for( j=1; j<N; j++ ){
				printf( ", %3d", out->Nodes->toReduced[out->Elements->Nodes[INDEX2(ref[j],i,M)]] );
			}
			printf( " ] DOF [ %3d", out->Nodes->reducedDegreeOfFreedom[out->Elements->Nodes[INDEX2(ref[0],i,M)]] );	
			for( j=1; j<N; j++ )
				printf( ", %3d", out->Nodes->reducedDegreeOfFreedom[out->Elements->Nodes[INDEX2(ref[j],i,M)]]);	
			printf( " ]\n" );	
		}

		printf( "\nFace Elements\n=================================================\n\n" );
		printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->FaceElements->elementDistribution->numInternal, out->FaceElements->elementDistribution->numBoundary, out->FaceElements->elementDistribution->numLocal );
		N = out->FaceElements->LinearReferenceElement->Type->numNodes;
		M = out->FaceElements->ReferenceElement->Type->numNodes;
		ref = out->FaceElements->ReferenceElement->Type->linearNodes;
		for( i=0; i<out->FaceElements->numElements; i++ ){
			printf( "face element %d \t: id %d  \t: dom %d  \t: nodes [ %3d", i, out->FaceElements->Id[i], out->FaceElements->Dom[i], out->Nodes->toReduced[out->FaceElements->Nodes[INDEX2(ref[0],i,M)]] );
			for( j=1; j<N; j++ )
				printf( ", %3d", out->Nodes->toReduced[out->FaceElements->Nodes[INDEX2(j,i,M)]]  );	
			printf( " ] DOF [ %3d", out->Nodes->reducedDegreeOfFreedom[out->FaceElements->Nodes[INDEX2(ref[0],i,M)]]);	
			for( j=1; j<N; j++ )
				printf( ", %3d", out->Nodes->reducedDegreeOfFreedom[out->FaceElements->Nodes[INDEX2(j,i,M)]]);	
			printf( " ]\n" );	
		}

		printf( "\nDistribution Data\n==================\n\n" );
		Finley_NodeDistribution_print( out->Nodes->reducedDegreeOfFreedomDistribution, stdout );
	}
}
#endif
/* 
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:51  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:32  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.3  2005/06/29 02:34:51  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:13  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/
