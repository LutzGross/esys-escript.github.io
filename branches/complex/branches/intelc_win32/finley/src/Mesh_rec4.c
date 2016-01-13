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

/*   Finley: generates rectangular meshes  */

/*   Generates a numElements[0] x numElements[1] mesh with first order elements (Rec4) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */


/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "RectangularMesh.h"

/**************************************************************/

#ifdef PASO_MPI
static void Finley_Mesh_printDistributed( Finley_Mesh *in )
{
  dim_t k, i0, NUMNODES;
  Finley_Mesh *out=in;

  NUMNODES = in->FaceElements->ReferenceElement->Type->numNodes;

  printf( "\n============NODES=============\n" );
  for( k=0; k<in->Nodes->numNodes; k++ )
    printf( "\tId %d\tDOF %d\tcoord [%3g%3g]\n", out->Nodes->Id[k], out->Nodes->degreeOfFreedom[k] , out->Nodes->Coordinates[INDEX2(0,k,2)], out->Nodes->Coordinates[INDEX2(1,k,2)] );
    for( k=0; k<out->Nodes->degreeOfFreedomDistribution->numNeighbours; k++ )
    {
      if( out->Nodes->degreeOfFreedomDistribution->neighbours[k]>=0 )
      {
        printf( "\t%d boundary DOF { ", out->Nodes->degreeOfFreedomDistribution->edges[k]->numForward ); 
        for( i0=0; i0<out->Nodes->degreeOfFreedomDistribution->edges[k]->numForward; i0++ )
          printf( "%d ", out->Nodes->degreeOfFreedomDistribution->edges[k]->indexForward[i0] );
        printf("} to %d\n", out->Nodes->degreeOfFreedomDistribution->neighbours[k] );
        printf( "\t%d external DOF { ", out->Nodes->degreeOfFreedomDistribution->edges[k]->numBackward ); 
        for( i0=0; i0<out->Nodes->degreeOfFreedomDistribution->edges[k]->numBackward; i0++ )
          printf( "%d ", out->Nodes->degreeOfFreedomDistribution->edges[k]->indexBackward[i0] );
        printf("} from %d\n", out->Nodes->degreeOfFreedomDistribution->neighbours[k] );
      }
    } 

  printf( "\n============ELEMENTS=============\n");
  for( k=0; k<in->Elements->elementDistribution->numInternal; k++ )
  {
    printf( "I\tId %d : nodes [%d %d %d %d]->DOF [%d %d %d %d]\n", out->Elements->Id[k], out->Elements->Nodes[INDEX2(0,k,4)], out->Elements->Nodes[INDEX2(1,k,4)], out->Elements->Nodes[INDEX2(2,k,4)], out->Elements->Nodes[INDEX2(3,k,4)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(2,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(3,k,4)]] );
  }

  for( k=in->Elements->elementDistribution->numInternal; k<in->Elements->elementDistribution->numInternal+in->Elements->elementDistribution->numBoundary; k++ )
  {
    printf( "B\tId %d : nodes [%d %d %d %d]->DOF [%d %d %d %d]\n", out->Elements->Id[k], out->Elements->Nodes[INDEX2(0,k,4)], out->Elements->Nodes[INDEX2(1,k,4)], out->Elements->Nodes[INDEX2(2,k,4)], out->Elements->Nodes[INDEX2(3,k,4)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(2,k,4)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(3,k,4)]] );
  }

  for( k=in->FaceElements->elementDistribution->numInternal; k<in->FaceElements->elementDistribution->numInternal+in->FaceElements->elementDistribution->numBoundary; k++ )
  {
    if( NUMNODES==4 )
      printf( "F\tId %d : nodes [%d %d %d %d]->DOF [%d %d %d %d]\n", out->FaceElements->Id[k], out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)], out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)], out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)], out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]] );
    else
      printf( "F\tId %d : nodes [%d %d]->DOF [%d %d]\n", out->FaceElements->Id[k], out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)], out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]] );
  }

}
#endif



/* get the number of nodes/elements for domain with rank=rank, of size processors
   where n is the total number of nodes/elements in the global domain */
static index_t domain_MODdim( index_t rank, index_t size, index_t n )
{ 
  rank = size-rank-1;

  if( rank < n%size )
    return (index_t)floor(n/size)+1;
  return (index_t)floor(n/size);
}


/* Determines the number of nodes/elements etc along an axis which is numElementsGlobal long for domain rank */
/* A bit messy, but it only has to be done once... */
static void domain_calculateDimension( index_t rank, dim_t size, dim_t numElementsGlobal, bool_t periodic, dim_t *numNodesLocal, dim_t *numDOFLocal, dim_t *numElementsLocal, dim_t *numElementsInternal, dim_t *firstNode, dim_t *nodesExternal, dim_t *DOFExternal, dim_t *numNodesExternal, bool_t *periodicLocal )
{
  index_t i0;
  dim_t numNodesGlobal = numElementsGlobal+1;

  (*numNodesLocal) = domain_MODdim( rank, size, numNodesGlobal );
  
  numElementsLocal[0] = numNodesLocal[0]+1;
  periodicLocal[0] = periodicLocal[1] = FALSE;
  nodesExternal[0] = nodesExternal[1] = 1;
  if( periodic ) 
  {
    if( size==1 )
    {
      numElementsLocal[0] = numElementsGlobal;
      nodesExternal[0] = nodesExternal[1] = 0;
      periodicLocal[0] = periodicLocal[1] = TRUE;
    }
    else
    {
      if( rank==0 )
      {
        periodicLocal[0] = TRUE;
        numNodesLocal[0]++;
      }
      if( rank==(size-1) )
      {      
        periodicLocal[1] = TRUE;
        numNodesLocal[0]--;
        numElementsLocal[0]--;
      }
    }
  }
  else if( !periodic )
  {
    if( rank==0 ){
      nodesExternal[0]--;
      numElementsLocal[0]--;
    }
    if( rank==(size-1) )
    {
      nodesExternal[1]--;
      numElementsLocal[0]--;
    }
  }
  numNodesExternal[0] = nodesExternal[0]+nodesExternal[1];
  
  numElementsInternal[0] = numElementsLocal[0];
  if( (rank==0) && (rank==size-1) );
    
  else if( !periodic && ( (rank==0) ^ (rank==size-1) ) )
      numElementsInternal[0] -= 1;
  else
    numElementsInternal[0] -= 2;

  firstNode[0] = 0;
  for( i0=0; i0<rank; i0++ )
    firstNode[0] += domain_MODdim( i0, size, numNodesGlobal );

  numDOFLocal[0] = numNodesLocal[0];
  if( periodicLocal[0] )
  {
    numDOFLocal[0]--;
  }
  DOFExternal[0] = nodesExternal[0];
  DOFExternal[1] = nodesExternal[1];

  /* some debugging printf statements */
  //printf( "rank/size = %d/%d\nNodes : %d Local, %d External[%d %d], First = %d\nElements : %d Local\nDOF : %d Local, External [%d %d]\nperiodicLocal [%d %d]\n\n", rank, size, *numNodesLocal, *numNodesExternal, nodesExternal[0], nodesExternal[1], *firstNode, *numElementsLocal, *numDOFLocal, DOFExternal[0], DOFExternal[1], periodicLocal[0], periodicLocal[1] );
}


Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) 
#ifndef PASO_MPI
{
  dim_t N0,N1,NE0,NE1,i0,i1,totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NUMNODES,M0,M1;
  index_t k,node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=NE0+1;
  N1=NE1+1;

  if (N0<=N1) {
     M0=1;
     M1=N0;
  } else {
     M0=N1;
     M1=1;
  }

  NFaceElements=0;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements+=2*NE1;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*NE0;
  } else {
      NDOF1=N1-1;
  }

  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d mesh",N0,N1);
  out=Finley_Mesh_alloc(name,2,order);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Rec4,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec4Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec4Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line2,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Line2_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0*N1);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  /*  set nodes: */
                                                                                                                                                                                                     
  #pragma omp parallel for private(i0,i1,k)
  for (i1=0;i1<N1;i1++) {
    for (i0=0;i0<N0;i0++) {
      k=M0*i0+M1*i1;
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(i0)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=i0+N0*i1;
      out->Nodes->Tag[k]=0;
      out->Nodes->degreeOfFreedom[k]=M0*(i0%NDOF0) +M1*(i1%NDOF1);
    }
  }
  /* tags for the faces: */
  if (!periodic[1]) {
    for (i0=0;i0<N0;i0++) {
      out->Nodes->Tag[M0*i0+M1*0]+=10;
      out->Nodes->Tag[M0*i0+M1*(N1-1)]+=20;
    }
  }
  if (!periodic[0]) {
    for (i1=0;i1<N1;i1++) {
      out->Nodes->Tag[M0*0+M1*i1]+=1;
      out->Nodes->Tag[M0*(N0-1)+M1*i1]+=2;
    }
  }

  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,k,node0) 
  for (i1=0;i1<NE1;i1++) {
    for (i0=0;i0<NE0;i0++) {
      k=i0+NE0*i1;
      node0=i0+i1*N0;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);

      out->Elements->Nodes[INDEX2(0,k,4)]=node0;
      out->Elements->Nodes[INDEX2(1,k,4)]=node0+1;
      out->Elements->Nodes[INDEX2(2,k,4)]=node0+N0+1;
      out->Elements->Nodes[INDEX2(3,k,4)]=node0+N0;

    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0);
  
  /*   face elements: */
  
  if (useElementsOnFace) {
     NUMNODES=4;
  } else {
     NUMNODES=2;
  }
  totalNECount=NE0*NE1;
  faceNECount=0;
  
  /* elements on boundary 010 (x2=0): */
  if (!periodic[1]) {
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
          k=i0+faceNECount;
          node0=i0;
   
          out->FaceElements->Id[k]=i0+totalNECount;
          out->FaceElements->Tag[k]=10;
          out->FaceElements->Color[k]=i0%2;
   
          if (useElementsOnFace) {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
              out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
              out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
          } else {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
          }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
     
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
         k=i0+faceNECount;
         node0=i0+(NE1-1)*N0;
   
         out->FaceElements->Id[k]=i0+totalNECount;
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Color[k]=i0%2+2;
   
         if (useElementsOnFace) {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
             out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
             out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
         } else {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
         }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
  }
  if (!periodic[0]) {
     /* **  elements on boundary 001 (x1=0): */
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
          k=i1+faceNECount;
          node0=i1*N0;

          out->FaceElements->Id[k]=i1+totalNECount;
          out->FaceElements->Tag[k]=1;
          out->FaceElements->Color[k]=i1%2+4;

          if (useElementsOnFace) {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
              out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
          } else {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
   
          }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  
     /* **  elements on boundary 002 (x1=1): */
  
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
          k=i1+faceNECount;
          node0=(NE0-1)+i1*N0;

          out->FaceElements->Id[k]=i1+totalNECount;
          out->FaceElements->Tag[k]=2;
          out->FaceElements->Color[k]=i1%2+6;
   
          if (useElementsOnFace) {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0;
             out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
          } else {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
          }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=7;

  
  /*  face elements done: */
  
  /*   condense the nodes: */
  
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculatuions:*/

  Finley_Mesh_prepare(out) ;

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  return out;
}

/*----------------------------------------------------------------------------
         MPI VERSION 


  ASSUMPTIONS
  ===========

  - the domain dimension is large enough in the NE0 direction for each domain
    to be 2 internal nodes wide in that direction. Face element calculation gets
    buggered up otherwise. 
  - if dividing into two domains with periodic[0]=TRUE , then the global domain 
    must be at least 4 elements along the NE0 direction, otherwise redundant 
    nodes are generated.

  These problems can be avoided by ensuring that numElements[0]/mpiSize >= 2

  if domains are small enough to cause these problems, you probably don't 
  need to solve them in parallel. If you really do, rotate the domain so that the
  long axis is aligned with NE0.
------------------------------------------------------------------------------*/


#else
{
  dim_t N0,N1,NE0,NE1,i0,i1, NE0_local, NDOF0,NDOF1, NFaceElements, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, faceNECount, totalNECount,M0,M1;
  index_t innerLoop=0;
  index_t NUMNODES,k,firstNode=0, DOFcount=0, forwardDOF[2], backwardDOF[2], node0, node1;
  index_t targetDomain=-1;
  index_t *indexForward=NULL;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=NE0+1;
  N1=NE1+1;
  Paso_MPIInfo *mpi_info = NULL;

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError())
        return NULL;

  if( mpi_info->rank==0 )
    domLeft = TRUE;
  if( mpi_info->rank==mpi_info->size-1 )
    domRight = TRUE;
  if( mpi_info->rank>0 && mpi_info->rank<mpi_info->size-1 )
    domInternal = TRUE;

  /* dimensions of the local subdomain */
  domain_calculateDimension( mpi_info->rank, mpi_info->size, NE0, periodic[0], &numNodesLocal, &numDOFLocal, &numElementsLocal, &numElementsInternal, &firstNode, nodesExternal, DOFExternal, &numNodesExternal, periodicLocal );  

  if (N0<=N1) {
     M0=1;
     M1=N0;
  } else {
     M0=N1;
     M1=1;
  }

  NFaceElements=0;
  if (!periodic[0]) 
  {
      if( domLeft )
        NFaceElements+=NE1;
      if( domRight )
        NFaceElements+=NE1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*NE0;
  } else {
      NDOF1=N1-1;
  }

  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d mesh",N0,N1);
  out=Finley_Mesh_alloc( name, 2, order, mpi_info );

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc( Rec4, out->order, mpi_info );
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec4Face,out->order, mpi_info );
     out->ContactElements=Finley_ElementFile_alloc(Rec4Face_Contact,out->order, mpi_info );
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line2,out->order, mpi_info );
     out->ContactElements=Finley_ElementFile_alloc(Line2_Contact,out->order, mpi_info );
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order, mpi_info );
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable( out->Nodes, (numNodesLocal + numNodesExternal)*N1 );
  Finley_NodeDistibution_allocTable( out->Nodes->degreeOfFreedomDistribution, numNodesLocal*NDOF1, numNodesExternal*NDOF1, 0 );
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  /*  set nodes: */
                                                                                                                                                                                                     
  //#pragma omp parallel for private(i0,i1,k)
  if( mpi_info->size==1 )
    innerLoop = numNodesLocal;
  else
    innerLoop = numNodesLocal-periodicLocal[0];
  
  for (k=0,i1=0;i1<N1;i1++) {
    for ( i0=0;i0<innerLoop;i0++,k++) {
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(i0+firstNode)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=i0+innerLoop*i1;
      out->Nodes->Tag[k]=0;
      out->Nodes->degreeOfFreedom[k]=i0%numDOFLocal + numDOFLocal*(i1%NDOF1);
    }
  }

  /* tag 'em */
  for( i0=0; i0<numNodesLocal; i0++ )
  {
    out->Nodes->Tag[i0] += 10;                        // bottom
    out->Nodes->Tag[i0+numNodesLocal*(N1-1)] += 20;   // top
  }
  if( domLeft && !periodicLocal[0] )
    for( i1=0; i1<N1; i1++ ) 
      out->Nodes->Tag[i1*numNodesLocal] += 1;         // left

  if( domRight && !periodicLocal[0] )
    for( i1=0; i1<N1; i1++ )
      out->Nodes->Tag[(i1+1)*numNodesLocal-1] += 2;   // right

  /* external nodes */
  k = innerLoop*N1;
  DOFcount = out->Nodes->degreeOfFreedom[k-1]+1;
  if( mpi_info->size>1 )
  {
    /* add forward communication information for boundary nodes */
    indexForward = TMPMEMALLOC( NDOF1, index_t );
    if( domInternal || domRight || periodicLocal[0] )
    {
        for( int i=0; i<NDOF1; i++ )
          indexForward[i] = i*numDOFLocal;
        targetDomain = mpi_info->rank-1>=0 ? mpi_info->rank-1 : mpi_info->size-1;
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexForward );
    }
    if( domInternal || domLeft || periodicLocal[1] )
    {
      for( int i=0; i<NDOF1; i++ )
          indexForward[i] = (i+1)*numDOFLocal - 1;
      targetDomain = (mpi_info->rank+1) % mpi_info->size;
      Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexForward );
    }
    TMPMEMFREE( indexForward );

    /* LHS */
    if( periodicLocal[0] )
    { 
      for (i1=0; i1<N1; i1++, k++) 
      {
        out->Nodes->Coordinates[INDEX2(0,k,2)]=Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k] = (i1 % NDOF1)*numDOFLocal;
      }
      out->Nodes->Tag[k-N1] = 10;
      out->Nodes->Tag[k-1]  = 20;
      for (i1=0; i1<N1; i1++, k++) 
      {
        out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(NE0-1)/DBLE(NE0)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k] = DOFcount + i1%NDOF1;
      }
      out->Nodes->Tag[k-N1] = 10;
      out->Nodes->Tag[k-1]  = 20;
      DOFcount = out->Nodes->degreeOfFreedom[k-1]+1;

      targetDomain = mpi_info->size-1;
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, out->Nodes->degreeOfFreedom+(numNodesLocal*N1) );
    }    
    if( mpi_info->rank>0 )
    {
      for (i1=0; i1<N1; i1++, k++) 
      {
        out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(firstNode-1)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k] = DOFcount + i1%NDOF1;
      }
      out->Nodes->Tag[k-N1] = 10;
      out->Nodes->Tag[k-1]  = 20;
      DOFcount = out->Nodes->degreeOfFreedom[k-1]+1;

      targetDomain = mpi_info->rank-1;
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, out->Nodes->degreeOfFreedom+(k-N1) );
    }
    
    /* RHS */
    if( periodicLocal[1] || (mpi_info->rank < mpi_info->size-1) )
    {
      for (i1=0; i1<N1; i1++, k++) 
      {
        out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(firstNode+numNodesLocal-periodicLocal[0])/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k] = DOFcount + i1%NDOF1;
      }
      out->Nodes->Tag[k-N1] = 10;
      out->Nodes->Tag[k-1]  = 20;
      DOFcount = out->Nodes->degreeOfFreedom[k-1]+1;
      
      targetDomain = (mpi_info->rank+1) % mpi_info->size;
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, out->Nodes->degreeOfFreedom+(k-N1) );
    }
  }

  /*   set the elements: */
  
 /* internal */
 //#pragma omp parallel for private(i0,i1,k,node0)  
 for ( i1=0, k=0; i1<NE1; i1++) 
 {  
    for ( i0=0; i0<numElementsInternal; i0++, k++) 
    {
      node0=i0+i1*innerLoop;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

      out->Elements->Nodes[INDEX2(0,k,4)]=node0;
      out->Elements->Nodes[INDEX2(1,k,4)]=node0+1;
      out->Elements->Nodes[INDEX2(2,k,4)]=node0+innerLoop+1;
      out->Elements->Nodes[INDEX2(3,k,4)]=node0+innerLoop;

    }
  }
  /* boundary */
  if( mpi_info->size>1 )
  {
    if( domInternal )
    {
      /* left */
      for( i1=0; i1<NE1; i1++, k++ )
      {
        node1=numNodesLocal*N1 + i1;
        node0=i1*numNodesLocal;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

        out->Elements->Nodes[INDEX2(0,k,4)]=node1;
        out->Elements->Nodes[INDEX2(1,k,4)]=node0;
        out->Elements->Nodes[INDEX2(2,k,4)]=node0+numNodesLocal;
        out->Elements->Nodes[INDEX2(3,k,4)]=node1+1;
      }

      /* right */
      for( i1=0; i1<NE1; i1++, k++ )
      {
        node0 = (i1+1)*numNodesLocal-1;
        node1 = (numNodesLocal+1)*N1 + i1;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

        out->Elements->Nodes[INDEX2(0,k,4)]=node0;
        out->Elements->Nodes[INDEX2(1,k,4)]=node1;
        out->Elements->Nodes[INDEX2(2,k,4)]=node1+1;
        out->Elements->Nodes[INDEX2(3,k,4)]=node0+numNodesLocal;
      }
    }
    else if( domLeft )
    {
      /* left */
      if( periodicLocal[0] )
      {
        node1 = numNodesLocal*N1;
        node0 = innerLoop*N1;
        for( i1=0; i1<NE1; i1++, k++, node1++, node0++ )
        {
          out->Elements->Id[k]=k;
          out->Elements->Tag[k]=0;
          out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

          out->Elements->Nodes[INDEX2(0,k,4)]=node1;
          out->Elements->Nodes[INDEX2(1,k,4)]=node0;
          out->Elements->Nodes[INDEX2(2,k,4)]=node0+1;
          out->Elements->Nodes[INDEX2(3,k,4)]=node1+1;
        }
      }
      /* right */
      for( i1=0; i1<NE1; i1++, k++ )
      {
        node0 = (i1+1)*innerLoop-1;
        node1 = (numNodesLocal+periodicLocal[0])*N1 + i1;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

        out->Elements->Nodes[INDEX2(0,k,4)]=node0;
        out->Elements->Nodes[INDEX2(1,k,4)]=node1;
        out->Elements->Nodes[INDEX2(2,k,4)]=node1+1;
        out->Elements->Nodes[INDEX2(3,k,4)]=node0+numNodesLocal-periodicLocal[0];
      }
    }
    else if( domRight )
    {
      /* left */
      for( i1=0; i1<NE1; i1++, k++ )
      {
        node1=numNodesLocal*N1 + i1;
        node0=i1*numNodesLocal;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

        out->Elements->Nodes[INDEX2(0,k,4)]=node1;
        out->Elements->Nodes[INDEX2(1,k,4)]=node0;
        out->Elements->Nodes[INDEX2(2,k,4)]=node0+numNodesLocal;
        out->Elements->Nodes[INDEX2(3,k,4)]=node1+1;
      }
      /* right */
      if( periodicLocal[1] )
      {
        for( i1=0; i1<NE1; i1++, k++ )
        {
          node0=(i1+1)*numNodesLocal-1;
          node1 = (numNodesLocal+1)*N1 + i1;

          out->Elements->Id[k]=k;
          out->Elements->Tag[k]=0;
          out->Elements->Color[k]=0; //COLOR_MOD(i0)+3*COLOR_MOD(i1);

          out->Elements->Nodes[INDEX2(0,k,4)]=node0;
          out->Elements->Nodes[INDEX2(1,k,4)]=node1;
          out->Elements->Nodes[INDEX2(2,k,4)]=node1+1;
          out->Elements->Nodes[INDEX2(3,k,4)]=node0+numNodesLocal;
        }
      }
    }  
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=0; //COLOR_MOD(0)+3*COLOR_MOD(0);

  if( domInternal )
  {
    out->Elements->elementDistribution->numInternal = NE1*(numElementsLocal-2);
    out->Elements->elementDistribution->numBoundary = NE1*2;
  }
  else if( mpi_info->size==1 )
  {
    out->Elements->elementDistribution->numInternal = NE1*numElementsLocal;
  }
  else
  {
      out->Elements->elementDistribution->numInternal += NE1*(numElementsLocal-1-periodic[0]);
      out->Elements->elementDistribution->numBoundary += NE1*(1 + periodic[0]);
  }
  
  out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal + out->FaceElements->elementDistribution->numBoundary;

  /*   face elements: */
  
  if (useElementsOnFace) {
     NUMNODES=4;
  } else {
     NUMNODES=2;
  }
  totalNECount=numElementsLocal*NE1;
  faceNECount=0;
  
  /* do the internal elements first */
  if (!periodic[1]) 
  {
     /* elements on boundary 010 (x1=0): */
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<numElementsLocal-numNodesExternal;i0++) {
          k=i0+faceNECount;
          node0=i0;
   
          out->FaceElements->Id[k]=i0+totalNECount;
          out->FaceElements->Tag[k]   = 10;
          out->FaceElements->Color[k] = 0; // i0%2;
   
          if (useElementsOnFace) {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
              out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal+1;
              out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numNodesLocal;
          } else {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
          }
     }
     totalNECount += numElementsLocal-numNodesExternal;
     faceNECount  += numElementsLocal-numNodesExternal;
     
     /* **  elements on boundary 020 (x1=1): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<numElementsLocal-numNodesExternal;i0++) {
         k=i0+faceNECount;
         node0=i0+(NE1-1)*numNodesLocal;
   
         out->FaceElements->Id[k]=i0+totalNECount;
         out->FaceElements->Tag[k]   = 20;
         out->FaceElements->Color[k] = 0; // i0%2+2;
   
         if (useElementsOnFace) {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numNodesLocal+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal;
             out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
             out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
         } else {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numNodesLocal+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal;
         }
     }
     totalNECount += numElementsLocal-numNodesExternal;
     faceNECount  += numElementsLocal-numNodesExternal;
  }
  if ( !periodic[0]) 
  {
     /* **  elements on boundary 001 (x0=0): */
     if( domLeft )
     {
       #pragma omp parallel for private(i1,k,node0) 
       for (i1=0;i1<NE1;i1++) {
            k=i1+faceNECount;
            node0=i1*numNodesLocal;

            out->FaceElements->Id[k]=i1+totalNECount;
            out->FaceElements->Tag[k]   = 1;
            out->FaceElements->Color[k] = 0; //i1%2+4;

            if (useElementsOnFace) {
                out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numNodesLocal;
                out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
                out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
                out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numNodesLocal+1;
            } else {
                out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numNodesLocal;
                out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
     
            }
       }
       totalNECount+=NE1;
       faceNECount+=NE1;
     }
     /* **  elements on boundary 002 (x0=1): */
     if( domRight )
     {
        #pragma omp parallel for private(i1,k,node0) 
        for (i1=0;i1<NE1;i1++) {
            k=i1+faceNECount;
            node0=(numNodesLocal-2)+i1*numNodesLocal;

            out->FaceElements->Id[k]=i1+totalNECount;
            out->FaceElements->Tag[k]   = 2;
            out->FaceElements->Color[k] = 0; //i1%2+6;
      
            if (useElementsOnFace) {
                out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
                out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal+1;
                out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal;
                out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
            } else {
                out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
                out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal+1;
            }
        }
        totalNECount+=NE1;
        faceNECount+=NE1;      
     }
  }

  if( mpi_info->size>1 )
{
  if( !periodic[1] && (domInternal || domRight) )
  {
    // bottom left
    k     = faceNECount;
    node0 = numNodesLocal*N1;

    out->FaceElements->Id[k]=totalNECount;
    out->FaceElements->Tag[k]   = 10;
    out->FaceElements->Color[k] = 0; //i1%2+6;

    if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=0;
        out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=numNodesLocal;
        out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
    } else {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=0;
    }
    totalNECount++;
    faceNECount++;

    // top left
    k     = faceNECount;
    node0 = (numNodesLocal+1)*N1 - 2;

    out->FaceElements->Id[k]=totalNECount;
    out->FaceElements->Tag[k]   = 20;
    out->FaceElements->Color[k] = 0; //i1%2+6;

    if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=numNodesLocal*(N1-1);
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
        out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
        out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=numNodesLocal*(N1-2);
    } else {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=numNodesLocal*(N1-1);
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
    }
    totalNECount++;
    faceNECount++;
  }
  if( !periodic[1] && (domInternal || domLeft) )
  { 
    // bottom right
    k     = faceNECount;
    node0 = (numNodesLocal+nodesExternal[0])*N1;

    out->FaceElements->Id[k]=totalNECount;
    out->FaceElements->Tag[k]   = 10;
    out->FaceElements->Color[k] = 0; //i1%2+6;

    if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=numNodesLocal-1;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
        out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
        out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=2*numNodesLocal-1;
    } else {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=numNodesLocal-1;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
    }
    totalNECount++;
    faceNECount++;

    // top right
    k     = faceNECount;
    node0 = (numNodesLocal+1+nodesExternal[0])*N1-2;

    out->FaceElements->Id[k]=totalNECount;
    out->FaceElements->Tag[k]   = 20;
    out->FaceElements->Color[k] = 0; //i1%2+6;

    if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=numNodesLocal*N1-1;
        out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=numNodesLocal*(N1-1)-1;
        out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
    } else {
        out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
        out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=numNodesLocal*N1-1;
    }
    totalNECount++;
    faceNECount++;
  }
}
  out->FaceElements->minColor = 0;
  out->FaceElements->maxColor = 0; //7;

  if( domInternal )
  {
    if( !periodic[1] )
    {
      out->FaceElements->elementDistribution->numInternal = 2*(numElementsLocal-2);
      out->FaceElements->elementDistribution->numBoundary = 4;
    }
  }
  else if( mpi_info->size==1 )
  {
    if( !periodic[0] )
      out->FaceElements->elementDistribution->numInternal += NE1*2;
    if( !periodic[1] )
      out->FaceElements->elementDistribution->numInternal += numElementsLocal*2;
  }
  else
  {
    if( !periodic[0] )
      out->FaceElements->elementDistribution->numInternal += NE1;
    if( !periodic[1] )
    {
      out->FaceElements->elementDistribution->numInternal += 2*(numElementsLocal-1-periodic[0]);
      out->FaceElements->elementDistribution->numBoundary += 2*(1 + periodic[0]);
    }
  }
  
  out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal + out->FaceElements->elementDistribution->numBoundary;

  if( out->FaceElements->elementDistribution->numLocal!=faceNECount )
    printf( "ERROR : face element numbers generated %d + %d = %d != %d\n",out->FaceElements->elementDistribution->numInternal, out->FaceElements->elementDistribution->numBoundary, out->FaceElements->elementDistribution->numLocal, faceNECount );

  

  /*  face elements done: */
  
  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds( out );

  /* prepare mesh for further calculatuions:*/

  //Finley_Mesh_prepare(out) ;

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  if ( !Finley_MPI_noError( mpi_info )) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  Paso_MPIInfo_dealloc( mpi_info );

  return out;
}
#endif



