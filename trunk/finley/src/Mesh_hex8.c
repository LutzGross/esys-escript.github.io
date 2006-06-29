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

/*   Finley: generates rectangular meshes */

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with first order elements (Hex8) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "RectangularMesh.h"

/**************************************************************/

#ifdef PASO_MPI
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
}

void print_mesh_statistics( Finley_Mesh *out )
{
  index_t i, j, N;
  
  printf( "\nNodes\n=====\n\n" );
	printf( "\t%d internal DOF\n\t%d boundary DOF\n\t%d local DOF\n\t%d external DOF\n", out->Nodes->degreeOfFreedomDistribution->numInternal, out->Nodes->degreeOfFreedomDistribution->numBoundary, out->Nodes->degreeOfFreedomDistribution->numLocal, out->Nodes->degreeOfFreedomDistribution->numExternal);
  for( i=0; i<out->Nodes->numNodes; i++ )
    printf( "node %d\t: id %d   \tDOF %d   \t: tag %d  \t: coordinates [%3g, %3g, %3g]\n", i, out->Nodes->Id[i], out->Nodes->degreeOfFreedom[i], out->Nodes->Tag[i], out->Nodes->Coordinates[INDEX2(0,i,3)], out->Nodes->Coordinates[INDEX2(1,i,3)], out->Nodes->Coordinates[INDEX2(2,i,3)] );

  printf( "Elements\n========\n\n" );
	printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->Elements->elementDistribution->numInternal, out->Elements->elementDistribution->numBoundary, out->Elements->elementDistribution->numLocal );
	N = out->Elements->ReferenceElement->Type->numNodes;
  for( i=0; i<out->Elements->numElements; i++ ){
    printf( "element %d    \t: id %d  \t: nodes [ %3d, %3d, %3d, %3d, %3d, %3d, %3d, %3d ]", i, out->Elements->Id[i], out->Elements->Nodes[INDEX2(0,i,8)], out->Elements->Nodes[INDEX2(1,i,8)], out->Elements->Nodes[INDEX2(2,i,8)], out->Elements->Nodes[INDEX2(3,i,8)], out->Elements->Nodes[INDEX2(4,i,8)], out->Elements->Nodes[INDEX2(5,i,8)], out->Elements->Nodes[INDEX2(6,i,8)], out->Elements->Nodes[INDEX2(7,i,8)] );
		printf( " DOF [ %3d", out->Nodes->degreeOfFreedom[out->Nodes->Id[out->Elements->Nodes[INDEX2(0,i,N)]]] );	
		for( j=1; j<N; j++ )
			printf( ", %3d", out->Nodes->degreeOfFreedom[out->Nodes->Id[out->Elements->Nodes[INDEX2(j,i,N)]]]  );	
		printf( " ]\n" );	
  }

	printf( "\nFace Elements\n==============\n\n" );
	printf( "\t%d internal\n\t%d boundary\n\t%d local\n", out->FaceElements->elementDistribution->numInternal, out->FaceElements->elementDistribution->numBoundary, out->FaceElements->elementDistribution->numLocal );
	N = out->FaceElements->ReferenceElement->Type->numNodes;
  for( i=0; i<out->FaceElements->numElements; i++ ){
		printf( "face element %d \t: id %d  \t: nodes [ %3d", i, out->FaceElements->Id[i], out->FaceElements->Nodes[INDEX2(0,i,N)] );
		for( j=1; j<N; j++ )
			printf( ", %3d", out->FaceElements->Nodes[INDEX2(j,i,N)]  );	
		printf( " ] DOF [ %3d", out->Nodes->degreeOfFreedom[out->Nodes->Id[out->FaceElements->Nodes[INDEX2(0,i,N)]]] );	
		for( j=1; j<N; j++ )
			printf( ", %3d", out->Nodes->degreeOfFreedom[out->Nodes->Id[out->FaceElements->Nodes[INDEX2(j,i,N)]]]  );	
		printf( " ]\n" );	
  }
}

#endif

#ifndef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace)
#else
Finley_Mesh* Finley_RectangularMesh_Hex8_singleCPU(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace)
#endif
{
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements,NUMNODES,M0,M1,M2;
  index_t node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=NE0+1;
  N1=NE1+1;
  N2=NE2+1;
#ifdef PASO_MPI
  Paso_MPIInfo *mpi_info = NULL;

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError())
        return NULL;
#endif

  if (N0<=MIN(N1,N2)) {
     if (N1 <= N2) {
        M0=1;
        M1=N0;
        M2=N0*N1;
     } else {
        M0=1;
        M2=N0;
        M1=N0*N2;
     }
  } else if (N1<=MIN(N2,N0)) {
     if (N2 <= N0) {
        M1=1;
        M2=N1;
        M0=N2*N1;
     } else {
        M1=1;
        M0=N1;
        M2=N1*N0;
     }
  } else {
     if (N0 <= N1) {
        M2=1;
        M0=N2;
        M1=N2*N0;
     } else {
        M2=1;
        M1=N2;
        M0=N1*N2;
     }
  }


  NFaceElements=0;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements+=2*NE1*NE2;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*NE0*NE2;
  } else {
      NDOF1=N1-1;
  }
  if (!periodic[2]) {
      NDOF2=N2;
      NFaceElements+=2*NE0*NE1;
  } else {
      NDOF2=N2-1;
  }
  
  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d x %d mesh",N0,N1,N2);

#ifndef PASO_MPI
  out=Finley_Mesh_alloc(name,3,order);
#else
  out=Finley_Mesh_alloc(name,3,order,mpi_info);
#endif
  if (! Finley_noError()) return NULL;

#ifdef PASO_MPI
  out->Elements=Finley_ElementFile_alloc(Hex8,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex8Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Hex8Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec4,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec4_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
#else
  out->Elements=Finley_ElementFile_alloc(Hex8,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex8Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Hex8Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec4,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec4_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */ 
  Finley_NodeFile_allocTable(out->Nodes,N0*N1*N2);
#ifdef PASO_MPI
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0*NDOF1*NDOF2, 0, 0 );
#endif
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */

  #pragma omp parallel for private(i0,i1,i2,k)
  for (i2=0;i2<N2;i2++) {
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0;i0++) {
        k=M0*i0+M1*i1+M2*i2;
        out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE(i0)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=i0+N0*i1+N0*N1*i2;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=M0*(i0%NDOF0) +M1*(i1%NDOF1) +M2*(i2%NDOF2);
      }
    }
  }
  /* tags for the faces: */
  /* tags for the faces: */
  if (!periodic[2]) {
     for (i1=0;i1<N1;i1++) {
       for (i0=0;i0<N0;i0++) {
         out->Nodes->Tag[M0*i0+M1*i1+M2*0]+=100;
         out->Nodes->Tag[M0*i0+M1*i1+M2*(N2-1)]+=200;
       }
     }
  }
  if (!periodic[1]) {
    for (i2=0;i2<N2;i2++) {
      for (i0=0;i0<N0;i0++) {
         out->Nodes->Tag[M0*i0+M1*0+M2*i2]+=10;
         out->Nodes->Tag[M0*i0+M1*(N1-1)+M2*i2]+=20;
      }
    }
  }
  if (!periodic[0]) {
    for (i2=0;i2<N2;i2++) {
      for (i1=0;i1<N1;i1++) {
        out->Nodes->Tag[M0*0+M1*i1+M2*i2]+=1;
        out->Nodes->Tag[M0*(N0-1)+M1*i1+M2*i2]+=2;
      }
    }
  }

  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,i2,k,node0) 
  for (i2=0;i2<NE2;i2++) {
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<NE0;i0++) {
        k=i0+NE0*i1+NE0*NE1*i2;
        node0=i0+i1*N0+N0*N1*i2;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);;

        out->Elements->Nodes[INDEX2(0,k,8)]=node0;
        out->Elements->Nodes[INDEX2(1,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(2,k,8)]=node0+N0+1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node0+N0;
        out->Elements->Nodes[INDEX2(4,k,8)]=node0+N0*N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0*N1+1;
        out->Elements->Nodes[INDEX2(6,k,8)]=node0+N0*N1+N0+1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0*N1+N0;

      }
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);
  
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=4;
  }
  totalNECount=NE0*NE1*NE2;
  faceNECount=0;
  
  /*   these are the quadrilateral elements on boundary 1 (x3=0): */
  if (!periodic[2]) {
     /* **  elements on boundary 100 (x3=0): */
  
     #pragma omp parallel for private(i0,i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i1+faceNECount;
         node0=i0+i1*N0;
   
         out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
         out->FaceElements->Tag[k]=100;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2);

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
         }
       }
     }
     totalNECount+=NE1*NE0;
     faceNECount+=NE1*NE0;
  
     /* **  elements on boundary 200 (x3=1) */
  
     #pragma omp parallel for private(i0,i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i1+faceNECount;
         node0=i0+i1*N0+N0*N1*(NE2-1);
   
         out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
         out->FaceElements->Tag[k]=200;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2)+4;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ N0 * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ N0 * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ N0 * N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ N0 * N1+N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ N0 * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ N0 * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ N0 * N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ N0 * N1+N0;
         }

   
       }
     }
     totalNECount+=NE1*NE0;
     faceNECount+=NE1*NE0;
  }
  if (!periodic[0]) {
     /* **  elements on boundary 001 (x1=0): */
  
     #pragma omp parallel for private(i1,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i1=0;i1<NE1;i1++) {
         k=i1+NE1*i2+faceNECount;
         node0=i1*N0+N0*N1*i2;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=1;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+8;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
         }
       }
     }
     totalNECount+=NE1*NE2;
     faceNECount+=NE1*NE2;
     
     /* **  elements on boundary 002 (x1=1): */
  
     #pragma omp parallel for private(i1,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i1=0;i1<NE1;i1++) {
         k=i1+NE1*i2+faceNECount;
         node0=(NE0-1)+i1*N0+N0*N1*i2 ;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=2;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+12;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0*N1+1;
         }
       }
     }
     totalNECount+=NE1*NE2;
     faceNECount+=NE1*NE2;
  }
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i2+faceNECount;
         node0=i0+N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+16;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1*N0+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*N0+N0;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*N0;
         }
       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i2+faceNECount;
         node0=i0+(NE1-1)*N0+N0*N1*i2;
   
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+20;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
         }

       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=23;
   
  /*  face elements done: */
  
#ifdef PASO_MPI
  /* make sure that the trivial distribution data is correct */
	out->FaceElements->elementDistribution->numBoundary = 0;
  out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal = faceNECount; 
	out->Elements->elementDistribution->numBoundary = 0;
  out->Elements->elementDistribution->numLocal = out->Elements->elementDistribution->numInternal = out->Elements->numElements; 
  out->ContactElements->elementDistribution->numLocal = out->ContactElements->elementDistribution->numInternal = out->ContactElements->elementDistribution->numInternal = 0;
  out->Points->elementDistribution->numLocal = out->Points->elementDistribution->numInternal = out->Points->elementDistribution->numInternal = 0;
	
	out->Nodes->degreeOfFreedomDistribution->numInternal = out->Nodes->degreeOfFreedomDistribution->numLocal;
  out->Nodes->degreeOfFreedomDistribution->numBoundary = 0;
#endif
  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds(out);

#ifdef PASO_MPI
  /* setup the CommBuffer */
  Finley_NodeDistribution_formCommBuffer( out->Nodes->degreeOfFreedomDistribution, out->Nodes->CommBuffer );
  if ( !Finley_MPI_noError( mpi_info )) {
    if( Finley_noError() )
      Finley_setError( PASO_MPI_ERROR, "Error on another MPI process" );
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  }
#endif

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

#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) 
{
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements,NUMNODES;//,M0,M1,M2;
  dim_t idCount, NE0_local, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal;
  bool_t dom_left, dom_right, dom_internal;

  index_t N0dom;
  index_t firstNode=0, DOFcount=0, forwardDOF[2], backwardDOF[2], node0, node1, node2;
  index_t targetDomain=-1;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE;
  index_t *indexForward=NULL;
  Finley_Mesh* out;

  char name[50];
  Paso_MPIInfo *mpi_info = NULL;
  double time0=Finley_timer();

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=NE0+1;
  N1=NE1+1;
  N2=NE2+1;


  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError())
        return NULL;

  /* use the serial version to generate the mesh for the 1-CPU case */
  if( mpi_info->size==1 )
  {
    Paso_MPIInfo_dealloc( mpi_info );
    out =  Finley_RectangularMesh_Hex8_singleCPU( numElements, Length, periodic, order, useElementsOnFace );
  	//print_mesh_statistics( out );
		return out;
  }    

  if( mpi_info->rank==0 )
    domLeft = TRUE;
  if( mpi_info->rank==mpi_info->size-1 )
    domRight = TRUE;
  if( mpi_info->rank>0 && mpi_info->rank<mpi_info->size-1 )
    domInternal = TRUE;

  /* dimensions of the local subdomain */
  domain_calculateDimension( mpi_info->rank, mpi_info->size, NE0, periodic[0], &numNodesLocal, &numDOFLocal, &numElementsLocal, &numElementsInternal, &firstNode, nodesExternal, DOFExternal, &numNodesExternal, periodicLocal );  

  /* count Degrees of Freedom along each axis */
  NDOF0 = N0 - periodic[0];
  NDOF1 = N1 - periodic[1];
  NDOF2 = N2 - periodic[2];

  /* count face elements */
  /* internal face elements */
  NFaceElements = 0;
  if( !periodic[0] )
    NFaceElements += (domLeft+domRight)*NE1*NE2;
  if( !periodic[1] )
    NFaceElements += 2*numElementsInternal*NE2;
  if( !periodic[2] )
    NFaceElements += 2*numElementsInternal*NE1;
  /* boundary face elements */
	/* this is looks nasty, but it beats a bunch of nested if/then/else carry-on */
  NFaceElements += 2*( 2 - (domLeft + domRight)*(!periodic[0]) )*( (!periodic[1])*NE2 + (!periodic[2])*NE1 );

  
  /*  allocate mesh: */
  sprintf(name,"Rectangular %d x %d x %d mesh",N0,N1,N2);

  out=Finley_Mesh_alloc(name,3,order,mpi_info);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex8,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex8Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Hex8Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec4,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec4_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,(numNodesLocal+2-!periodic[0]*(domLeft+domRight))*N1*N2);
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, numDOFLocal*NDOF1*NDOF2, (DOFExternal[0]+DOFExternal[1])*NDOF1*NDOF2, 0 );
  Finley_ElementFile_allocTable(out->Elements,(numElementsLocal)*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */
  /* INTERNAL/BOUNDARY NODES */
	k=0;
  #pragma omp parallel for private(i0,i1,i2,k)
  for (i2=0;i2<N2;i2++) { 
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<numNodesLocal-domLeft*periodic[0];i0++,k++) {         
        out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE(i0+firstNode)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=k;
      }
    }
  }
  if( domLeft && periodic[0] ) {
    for (i2=0;i2<N2;i2++)  {
      for (i1=0;i1<N1;i1++, k++) {
        out->Nodes->Coordinates[INDEX2(0,k,3)]=Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=i1*numNodesLocal + i2*numNodesLocal*N1;
      }
    }
    /* tag the faces for this special case */
    if( !periodic[1] )
    {
      for( i2=0; i2<N2; i2++ ){
        out->Nodes->Tag[k + (i2-N2)*N1     ] += 10;
        out->Nodes->Tag[k + (i2+1-N2)*N1 -1] += 20; 
      }
    }
    if( !periodic[2] )
    {
      for( i1=0; i1<N1; i1++ ){
        out->Nodes->Tag[k -N1*N2 +i1] += 100;
        out->Nodes->Tag[k -N1    +i1] += 200;       
      }
    }
  }
  /* tags for the faces: */
  N0dom = (numNodesLocal-periodicLocal[0]);
  if (!periodic[2]) {    
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0dom;i0++) {   
         out->Nodes->Tag[i0 + N0dom*i1]+=100;
         out->Nodes->Tag[i0 + N0dom*i1 + N0dom*N1*(N2-1)]+=200;
       }
     }
  }
  if (!periodic[1]) {
    for (i2=0;i2<N2;i2++) {
      for (i0=0;i0<N0dom;i0++) {
         out->Nodes->Tag[i0 + i2*N1*N0dom]+=10;
         out->Nodes->Tag[i0 + (i2+1)*N1*N0dom-N0dom]+=20;
      }
    }
  }
  if (!periodic[0] && !domInternal ) {
    for (i2=0;i2<N2;i2++) {
      for (i1=0;i1<N1;i1++) {
        if( domLeft )
          out->Nodes->Tag[i1*N0dom + i2*N0dom*N1]+=1;
        if( domRight )
          out->Nodes->Tag[(i1+1)*N0dom-1 + i2*N0dom*N1]+=2;
      }
    }
  }
	/* setup the forward communication data for the boundary nodes that we have just defined */
	/* the case where there are 2 subdomains and periodic[0]=true has to be treated
		 as a special case to because the two domains have two interface boundaries to one-another */
  indexForward = TMPMEMALLOC( NDOF1*NDOF2, index_t );
	if( mpi_info->size>2 || !periodic[0] ){
		if( domInternal || domRight || periodicLocal[0] )
		{
				for( int i=0; i<NDOF2; i++ )
					for( int j=0; j<NDOF1; j++ )
						indexForward[j+i*NDOF1] = numDOFLocal*j+NDOF1*numDOFLocal*i;
				targetDomain = mpi_info->rank-1>=0 ? mpi_info->rank-1 : mpi_info->size-1;
				Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );
		}
		if( domInternal || domLeft || periodicLocal[1] )
		{
			for( int i=0; i<NDOF2; i++ )
				for( int j=0; j<NDOF1; j++ )
					indexForward[j+i*NDOF1] = numDOFLocal*(j+1)-1+NDOF1*numDOFLocal*i;
			targetDomain = (mpi_info->rank+1) % mpi_info->size;
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );
		}
	}
	else {
			for( int i=0; i<NDOF2; i++ )
				for( int j=0; j<NDOF1; j++ )
					indexForward[j+i*NDOF1] = numDOFLocal*(j+1)-1+NDOF1*numDOFLocal*i;
			targetDomain = (mpi_info->rank+1) % mpi_info->size;
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );

			for( int i=0; i<NDOF2; i++ )
				for( int j=0; j<NDOF1; j++ )
					indexForward[j+i*NDOF1] = numDOFLocal*j+NDOF1*numDOFLocal*i;
			targetDomain = mpi_info->rank-1>=0 ? mpi_info->rank-1 : mpi_info->size-1;
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );
	}
	
  /* EXTERNAL NODES */
  /* left hand boundary */
  DOFcount = NDOF1*NDOF2*numDOFLocal;
  if( (domLeft && periodic[0]) || !domLeft ) {
    if( (domLeft && periodic[0]) )
      for (i2=0;i2<N2;i2++)  {
        for (i1=0;i1<N1;i1++, k++) {
          out->Nodes->Coordinates[INDEX2(0,k,3)]=(1.-DBLE(1)/DBLE(N0-1))*Length[0];
          out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
          out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
          out->Nodes->Id[k]=k;
          out->Nodes->Tag[k]=0;
					out->Nodes->degreeOfFreedom[k]=DOFcount+i1%NDOF1+(i2%NDOF2)*NDOF1;
        }
      }
    else
      for (i2=0;i2<N2;i2++)  {
        for (i1=0;i1<N1;i1++, k++) {
          out->Nodes->Coordinates[INDEX2(0,k,3)]=(DBLE(firstNode-1)/DBLE(N0-1))*Length[0];
          out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
          out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
          out->Nodes->Id[k]=k;
          out->Nodes->Tag[k]=0;
					out->Nodes->degreeOfFreedom[k]=DOFcount+i1%NDOF1+(i2%NDOF2)*NDOF1;
        }
      }      
		DOFcount += NDOF1*NDOF2;
		targetDomain = mpi_info->rank-1>=0 ? mpi_info->rank-1 : mpi_info->size-1;
		if( !periodic[1] ){
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, out->Nodes->degreeOfFreedom + (k-NDOF1*NDOF2) );
		}
		else {
			for( int i=0; i<NDOF2; i++ )
				for( int j=0; j<NDOF1; j++ )
					indexForward[j+i*NDOF1] = DOFcount - NDOF1*NDOF2 + j + i*NDOF1; 
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );
		}
				
    /* tag the faces for this special case */
    if( !periodic[1] )
    {
      for( i1=0; i1<N1; i1++ ){
        out->Nodes->Tag[k -N1*N2 +i1] += 10;
        out->Nodes->Tag[k -N1    +i1] += 20; 
      }
    }
    if( periodic[2] )
    {
      for( i2=0; i2<N2; i2++ ){
        out->Nodes->Tag[k +(i2-N2)*N1     ] += 100;
        out->Nodes->Tag[k +(i2-N2+1)*N1 -1] += 200;       
      }
    }
  }
  if( (domRight && periodic[0]) || !domRight )
  {
    if( domRight && periodic[0] )
      for (i2=0;i2<N2;i2++)  {
        for (i1=0;i1<N1;i1++, k++) {
          out->Nodes->Coordinates[INDEX2(0,k,3)]=Length[0];
          out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
          out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
          out->Nodes->Id[k]=k;
          out->Nodes->Tag[k]=0;
					out->Nodes->degreeOfFreedom[k]=DOFcount+i1%NDOF1+(i2%NDOF2)*NDOF1;
        }
      }
    else
      for (i2=0;i2<N2;i2++)  {
        for (i1=0;i1<N1;i1++, k++) {
          out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE(firstNode+numNodesLocal-periodicLocal[0])/DBLE(N0-1)*Length[0];
          out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
          out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
          out->Nodes->Id[k]=k;
          out->Nodes->Tag[k]=0;
					out->Nodes->degreeOfFreedom[k]=DOFcount+i1%NDOF1+(i2%NDOF2)*NDOF1;
        }
      }
		DOFcount += NDOF1*NDOF2;

		targetDomain = mpi_info->rank+1 < mpi_info->size? mpi_info->rank+1 : 0;
		if( !periodic[1] ){
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, out->Nodes->degreeOfFreedom + (k-NDOF1*NDOF2) );
		}
		else {
			for( int i=0; i<NDOF2; i++ )
				for( int j=0; j<NDOF1; j++ )
					indexForward[j+i*NDOF1] = DOFcount - NDOF1*NDOF2 + j + i*NDOF1; 
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, indexForward );
		}
			
    /* tag the faces for this special case */
    if( !periodic[1] )
    {
      for( i1=0; i1<N1; i1++ ){
        out->Nodes->Tag[k -N1*N2 +i1] += 10;
        out->Nodes->Tag[k -N1    +i1] += 20; 
      }
    }
    if( !periodic[2] )
    {
      for( i2=0; i2<N2; i2++ ){
        out->Nodes->Tag[k +(i2-N2)*N1     ] += 100;
        out->Nodes->Tag[k +(i2-N2+1)*N1 -1] += 200;       
      }
    }
  }
	out->Nodes->degreeOfFreedomDistribution->numInternal = NDOF1*NDOF2*(numDOFLocal - 2 + domRight*(!periodic[0]) + domLeft*(!periodic[0]));
	out->Nodes->degreeOfFreedomDistribution->numBoundary = out->Nodes->degreeOfFreedomDistribution->numLocal - out->Nodes->degreeOfFreedomDistribution->numInternal; 
  
	TMPMEMFREE( indexForward );
  /*   set the elements: */

  /* INTERNAL elements */
  N0dom = (numNodesLocal-periodicLocal[0]);
  k = 0;
  #pragma omp parallel for private(i0,i1,i2,k,node0) 
  for (i2=0;i2<NE2;i2++) {
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<numElementsInternal;i0++,k++) {
        node0=i0+i1*N0dom+N0dom*N1*i2;

        out->Elements->Id[k]=k;
				
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0;//COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);;

        out->Elements->Nodes[INDEX2(0,k,8)]=node0;
        out->Elements->Nodes[INDEX2(1,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(2,k,8)]=node0+N0dom+1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node0+N0dom;
        out->Elements->Nodes[INDEX2(4,k,8)]=node0+N0dom*N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0dom*N1+1;
        out->Elements->Nodes[INDEX2(6,k,8)]=node0+N0dom*N1+N0dom+1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0dom*N1+N0dom;

      }
    }
  }
	out->Elements->elementDistribution->numInternal = NE1*NE2*numElementsInternal;
	out->Elements->elementDistribution->numBoundary = 0;

  /* BOUNDARY Elements */
  /* left boundary */
  if( !domLeft )
  {
    for (i2=0;i2<NE2;i2++) {
      node0 = numNodesLocal*N1*N2 + i2*N1;
      for (i1=0;i1<NE1;i1++,node0++,k++) {
        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0;//COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);

        out->Elements->Nodes[INDEX2(0,k,8)]=node0;
        out->Elements->Nodes[INDEX2(1,k,8)]=i1*numNodesLocal + i2*numNodesLocal*N1;
        out->Elements->Nodes[INDEX2(2,k,8)]=(i1+1)*numNodesLocal + i2*numNodesLocal*N1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(4,k,8)]=node0+N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=i1*numNodesLocal + (i2+1)*numNodesLocal*N1;
        out->Elements->Nodes[INDEX2(6,k,8)]=(i1+1)*numNodesLocal + (i2+1)*numNodesLocal*N1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node0+N1+1;    
      }
    }
		out->Elements->elementDistribution->numBoundary += NE1*NE2;
  }
  /* the left periodic boundary is done a little differently to a left internal boundary */
  else if( (domLeft && periodic[0]) )
  {
    for (i2=0;i2<NE2;i2++) {
      node0 = numDOFLocal*N1*N2 + i2*N1;
      node1 = node0 + N1*N2;
      for (i1=0;i1<NE1;i1++,k++,node0++,node1++) {

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0;//COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);

        out->Elements->Nodes[INDEX2(0,k,8)]=node1;
        out->Elements->Nodes[INDEX2(1,k,8)]=node0;
        out->Elements->Nodes[INDEX2(2,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node1+1;
        out->Elements->Nodes[INDEX2(4,k,8)]=node1+N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=node0+N1;
        out->Elements->Nodes[INDEX2(6,k,8)]=node0+N1+1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node1+N1+1;    
      }
    }  
		out->Elements->elementDistribution->numBoundary += NE1*NE2;
  }
  /* right boundary */
  if( !domRight || (domRight && periodic[0]) ){
    for (i2=0;i2<NE2;i2++) {
      for (i1=0;i1<NE1;i1++,node0++,node1+=numDOFLocal,k++) {
				node1 = numDOFLocal -1 + numDOFLocal*i1 + N1*numDOFLocal*i2;
				node0 = (numNodesLocal + domInternal + periodicLocal[0])*N1*N2 + i2*N1 + i1;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=0;//COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);;

        out->Elements->Nodes[INDEX2(0,k,8)]=node1;
        out->Elements->Nodes[INDEX2(1,k,8)]=node0;
        out->Elements->Nodes[INDEX2(2,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node1+N0dom;
        out->Elements->Nodes[INDEX2(4,k,8)]=node1+N0dom*N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=node0+N1;
        out->Elements->Nodes[INDEX2(6,k,8)]=node0+N1+1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node1+N0dom*N1+N0dom;   
      }
    }
		out->Elements->elementDistribution->numBoundary += NE1*NE2;
	}

	out->Elements->minColor=0;
  out->Elements->maxColor=0;//COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);
	out->Elements->elementDistribution->numLocal = out->Elements->elementDistribution->numInternal + out->Elements->elementDistribution->numBoundary; 
	
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=4;
  }
  totalNECount=k;
  faceNECount=0;
	idCount = totalNECount;
  
	/*   Do internal face elements for each boundary face first */
  /*   these are the quadrilateral elements on boundary 1 (x3=0): */
	numElementsInternal = numElementsLocal-nodesExternal[0]-nodesExternal[1];
  if (!periodic[2]) {
     /*  elements on boundary 100 (x3=0): */
  
     #pragma omp parallel for private(i0,i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0; i0<numElementsInternal; i0++) {
         k=i0+numElementsInternal*i1+faceNECount;
         node0=i0+i1*numDOFLocal;
   
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=100;
         out->FaceElements->Color[k]=0;//(i0%2)+2*(i1%2);

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+numDOFLocal*N1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numDOFLocal*N1+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
         }
       }
     }
     totalNECount+=NE1*numElementsInternal;
     faceNECount+=NE1*numElementsInternal;
 		 
     /* **  elements on boundary 200 (x3=1) */
  
     #pragma omp parallel for private(i0,i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0;i0<numElementsInternal;i0++) {
         k=i0+numElementsInternal*i1+faceNECount;
         node0=i0+i1*numDOFLocal+numDOFLocal*N1*(NE2-1);
   
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=200;
         out->FaceElements->Color[k]=0;//(i0%2)+2*(i1%2)+4;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ numDOFLocal * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ numDOFLocal * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ numDOFLocal * N1+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ numDOFLocal * N1+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numDOFLocal;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ numDOFLocal * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ numDOFLocal * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ numDOFLocal * N1+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ numDOFLocal * N1+numDOFLocal;
         }
       }
     }
     totalNECount+=NE1*numElementsInternal;
     faceNECount+=NE1*numElementsInternal;
  }
  if (!periodic[0] && !domInternal) {
     /* **  elements on boundary 001 (x1=0): */
  	 if( domLeft ){
			 #pragma omp parallel for private(i1,i2,k,node0) 
			 for (i2=0;i2<NE2;i2++) {
				 for (i1=0;i1<NE1;i1++) {
					 k=i1+NE1*i2+faceNECount;
					 node0=i1*numDOFLocal+numDOFLocal*N1*i2;
		 
					 out->FaceElements->Id[k]=idCount++;
					 out->FaceElements->Tag[k]=1;
					 out->FaceElements->Color[k]=0; //(i2%2)+2*(i1%2)+8;

					 if  (useElementsOnFace) {
							out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
							out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal*N1;
							out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
							out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal;
							out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
							out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal*N1+1;
							out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
							out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numDOFLocal+1;
					 } else {
							out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
							out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal*N1;
							out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
							out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal;
					 }
				 }
			 }
			 totalNECount+=NE1*NE2;
			 faceNECount+=NE1*NE2;
     }
     /* **  elements on boundary 002 (x1=1): */
 		 if( domRight ) { 
			 #pragma omp parallel for private(i1,i2,k,node0) 
			 for (i2=0;i2<NE2;i2++) {
				 for (i1=0;i1<NE1;i1++) {
					 k=i1+NE1*i2+faceNECount;
					 node0=(numDOFLocal-2)+i1*numDOFLocal+numDOFLocal*N1*i2 ;
		 
					 out->FaceElements->Id[k]=idCount++;
					 out->FaceElements->Tag[k]=2;
					 out->FaceElements->Color[k]=0;//(i2%2)+2*(i1%2)+12;

					 if  (useElementsOnFace) {
							out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
							out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal+1;
							out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
							out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal*N1+1;
							out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
							out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal;
							out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
							out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numDOFLocal*N1;
					 } else {
							out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
							out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal+1;
							out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
							out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal*N1+1;
					 }
				 }
		 	 }
       totalNECount+=NE1*NE2;
       faceNECount+=NE1*NE2;
     }
  }
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<numElementsInternal;i0++) {
         k=i0+numElementsInternal*i2+faceNECount;
         node0=i0+numDOFLocal*N1*i2;
   
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Color[k]=0;//(i0%2)+2*(i2%2)+16;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*numDOFLocal;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1*numDOFLocal+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*numDOFLocal+numDOFLocal;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*numDOFLocal;
         }
       }
     }
     totalNECount+=numElementsInternal*NE2;
     faceNECount+=numElementsInternal*NE2;
  
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<numElementsInternal;i0++) {
         k=i0+numElementsInternal*i2+faceNECount;
         node0=i0+(NE1-1)*numDOFLocal+numDOFLocal*N1*i2;
   
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+20;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal*N1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numDOFLocal*N1+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal+1;
         }

       }
     }
     totalNECount+=numElementsInternal*NE2;
     faceNECount+=numElementsInternal*NE2;
  }
	out->FaceElements->elementDistribution->numInternal = faceNECount;
	
	/* now do the boundary face elements */
	/* LHS */
	if( !domLeft  )
	{
		if( !periodic[2] ) {

			/* x3=0 */
			for( i1=0; i1<NE1; i1++ )
			{
				k = i1+faceNECount;
				node0 = i1*numNodesLocal;
				node1 = numNodesLocal*N1*N2 + i1;

				out->FaceElements->Tag[k]=200;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numNodesLocal*N1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
				}
			}	
			faceNECount  += NE1;
			totalNECount += NE1;

			/* x3=1 */
			for( i1=0; i1<NE1; i1++ )
			{
				k = i1+faceNECount;
				node0 = numNodesLocal*N1*(NE2-1) + i1*numNodesLocal;
				node1 = numNodesLocal*N1*N2 + i1 + (NE2-1)*N1;

				out->FaceElements->Tag[k]=200;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal*N1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1+1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numNodesLocal*N1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1+1;
				}
			}	
			faceNECount  += NE1;
			totalNECount += NE1;
		}

		if( !periodic[1] ) {
			/* x2=0 */
			for( i2=0; i2<NE2; i2++ )
			{
				k = i2+faceNECount;
				node0 = i2*numNodesLocal*N1;
				node1 = numNodesLocal*N1*N2 + i2*N1;

				out->FaceElements->Tag[k]=20;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1+N1+1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1;
				}
			}	
			faceNECount  += NE2;
			totalNECount += NE2;

			/* x2=1 */
			for( i2=0; i2<NE2; i2++ )
			{
				k = i2+faceNECount;
				node0 = i2*numNodesLocal*N1 + numNodesLocal*(NE1-1); 
				node1 = numNodesLocal*N1*N2 + i2*N1 + (NE1-1);

				out->FaceElements->Tag[k]=20;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+numNodesLocal*N1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1+1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+numNodesLocal*N1+numNodesLocal;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1+1;
				}
			}	
			faceNECount  += NE2;
			totalNECount += NE2;
		}
	}
  
	/* RHS */
	if( !domRight || periodicLocal[1] )
	{
		/* the case of left hand boundary domain and periodic boundary condition on its left hand boundary */
		if( domLeft && periodic[0] ){
			if( !periodic[2] ) {

				/* x3=0 */
				for( i1=0; i1<NE1; i1++ )
				{
					k = i1+faceNECount;
					node0 = numDOFLocal*N1*N2 + i1; 
					node1 = numNodesLocal*N1*N2 + i1;

					out->FaceElements->Tag[k]=200;
					out->FaceElements->Id[k]=idCount++;
					out->FaceElements->Color[k]=0;

					if( useElementsOnFace )
					{
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1+N1;
						out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1+N1+1;
						out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1;
					} else {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
					}
				}	
				faceNECount  += NE1;
				totalNECount += NE1;

				/* x3=1 */
				for( i1=0; i1<NE1; i1++ )
				{
					k = i1+faceNECount;
					node0 = numDOFLocal*N1*N2 + i1 + (NE2-1)*N1;
					node1 = numNodesLocal*N1*N2 + i1 + (NE2-1)*N1;

					out->FaceElements->Tag[k]=200;
					out->FaceElements->Id[k]=idCount++;
					out->FaceElements->Color[k]=0;
					
					if( useElementsOnFace )
					{
						out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+N1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N1;
					} else {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+N1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N1;
					}
				}	
				faceNECount  += NE1;
				totalNECount += NE1;
			}

			if( !periodic[1] ) {
				/* x2=0 */
				for( i2=0; i2<NE2; i2++ )
				{
					k = i2+faceNECount;
					node0 = numDOFLocal*N1*N2 + i2*N1; 
					node1 = numNodesLocal*N1*N2 + i2*N1;

					out->FaceElements->Tag[k]=20;
					out->FaceElements->Id[k]=idCount++;
					out->FaceElements->Color[k]=0;

					if( useElementsOnFace )
					{
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1;
						out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1+N1+1;
					} else {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+N1;
					}
				}	
				faceNECount  += NE2;
				totalNECount += NE2;

				/* x2=1 */
				for( i2=0; i2<NE2; i2++ )
				{
					k = i2+faceNECount;
					node0 = numDOFLocal*N1*N2 + i2*N1 + N1-2; 
					node1 = numNodesLocal*N1*N2 + i2*N1 + N1-2;

					out->FaceElements->Tag[k]=20;
					out->FaceElements->Id[k]=idCount++;
					out->FaceElements->Color[k]=0;

					if( useElementsOnFace )
					{
						out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node1;
						out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1;
						out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1+N1;
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1+1;
					} else {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node1+1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1+1;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1+1;
					}
				}	
				faceNECount  += NE2;
				totalNECount += NE2;
			}
			
		}
		if( !periodic[2] ) {
			/* x3=0 */
			for( i1=0; i1<NE1; i1++ )
			{
				k = i1+faceNECount;
				node0 = numDOFLocal*(i1+1) - 1;
				node1 = (numNodesLocal+periodicLocal[0]+periodicLocal[1])*N1*N2 + i1 + domInternal*N1*N2;

				out->FaceElements->Tag[k]=200;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+numDOFLocal*N1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1+N1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1;
				}
			}	
			faceNECount  += NE1;
			totalNECount += NE1;

			/* x3=1 */
			for( i1=0; i1<NE1; i1++ )
			{
				k = i1+faceNECount;
				node0 = numDOFLocal*N1*(NE2-1) + (i1+1)*numDOFLocal - 1;
				node1 = numNodesLocal*N1*N2 + i1 + (NE2-1)*N1 + (domInternal+periodicLocal[1]+periodicLocal[0])*N1*N2;

				out->FaceElements->Tag[k]=200;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;
				
				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal*N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal*N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+numDOFLocal*N1+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1+N1;
				}
			}	
			faceNECount  += NE1;
			totalNECount += NE1;
		}
		if( !periodic[1] ) {
			/* x2=0 */
			for( i2=0; i2<NE2; i2++ )
			{
				k = i2+faceNECount;
				node0 = N1*numDOFLocal*i2 + numDOFLocal - 1;
				node1 = numNodesLocal*N1*N2 + i2*N1 + (domInternal+periodicLocal[0]+periodicLocal[1])*N1*N2;

				out->FaceElements->Tag[k]=20;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;

				if( useElementsOnFace )
				{
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*numDOFLocal;
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*numDOFLocal+numDOFLocal;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*numDOFLocal;
				}
			}	
			faceNECount  += NE2;
			totalNECount += NE2;

			/* x2=1 */
			for( i2=0; i2<NE2; i2++ )
			{
				k = i2+faceNECount;
				node0 = numDOFLocal*N1*i2 + NE1*numDOFLocal - 1;
				node1 = numNodesLocal*N1*N2 + i2*N1 + (NE1-1) + (domInternal+periodicLocal[0]+periodicLocal[1])*N1*N2;

				out->FaceElements->Tag[k]=20;
				out->FaceElements->Id[k]=idCount++;
				out->FaceElements->Color[k]=0;
				
				if( useElementsOnFace ){
					out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node1;
					out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node1+N1;
					out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N1*numDOFLocal;
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N1*numDOFLocal+numDOFLocal;
				} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+numDOFLocal;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node1+1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node1+N1+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N1*numDOFLocal+numDOFLocal;
				}
			}
			faceNECount  += NE2;
			totalNECount += NE2;
		}
	}
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=0;//23;

	out->FaceElements->elementDistribution->numBoundary = faceNECount - out->FaceElements->elementDistribution->numInternal;
  out->FaceElements->elementDistribution->numLocal = faceNECount; 


  /* setup distribution info for other elements */
  out->ContactElements->elementDistribution->numLocal = out->ContactElements->elementDistribution->numInternal = out->ContactElements->elementDistribution->numInternal = 0;
  out->Points->elementDistribution->numLocal = out->Points->elementDistribution->numInternal = out->Points->elementDistribution->numInternal = 0;
	
  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds( out );

  /* setup the CommBuffer */
  Finley_NodeDistribution_formCommBuffer( out->Nodes->degreeOfFreedomDistribution, out->Nodes->CommBuffer );
  if ( !Finley_MPI_noError( mpi_info )) {
    if( Finley_noError() )
      Finley_setError( PASO_MPI_ERROR, "Error on another MPI process" );
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  }

  Finley_NodeDistribution_calculateIndexExternal( out->Nodes->degreeOfFreedomDistribution, out->Nodes->CommBuffer );

  /* prepare mesh for further calculatuions:*/
  Finley_Mesh_prepare(out) ;

//  print_mesh_statistics( out );

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif
	
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  return out;
}
#endif

