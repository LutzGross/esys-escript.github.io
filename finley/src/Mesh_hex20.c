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

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with second order elements (Hex20) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */


/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$

/**************************************************************/

#include "RectangularMesh.h"

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
static void domain_calculateDimension( index_t rank, dim_t size, dim_t numElementsGlobal, bool_t periodic, dim_t *numNodesLocal, dim_t *numDOFLocal, dim_t *numElementsLocal, dim_t *numElementsInternal, dim_t *firstNode, dim_t *nodesExternal, dim_t *DOFExternal, dim_t *numNodesExternal, bool_t *periodicLocal, dim_t *DOFBoundary )
{
  index_t i0;
  dim_t numNodesGlobal = numElementsGlobal+1;

  (*numNodesLocal) = 2*domain_MODdim( rank, size, numNodesGlobal ) - 1;
  if( rank<size-1 ) // add on node for right hand boundary
    (*numNodesLocal) += 1;

  numElementsLocal[0] = domain_MODdim( rank, size, numNodesGlobal )+1;
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
  nodesExternal[0]*=2;
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
  firstNode[0] *= 2;
  
  numDOFLocal[0] = numNodesLocal[0];
  if( periodicLocal[0] )
    numDOFLocal[0]--;
  
  DOFExternal[0] = nodesExternal[0];
  DOFExternal[1] = nodesExternal[1];

  if( size>1 )
  {
    DOFBoundary[0] = periodicLocal[0]*1 + (rank>0)*1;
    DOFBoundary[1] = periodicLocal[1]*2 + (rank<(size-1))*2;
  }
  else
  {
    DOFBoundary[0] = DOFBoundary[1] = 0;
  }
}
#endif
/**************************************************************/
#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Hex20_singleCPU(dim_t* numElements,double* Length,bool_t* periodic,index_t order,bool_t useElementsOnFace, Paso_MPIInfo *mpi_info) {
#else
Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t* numElements,double* Length,bool_t* periodic,index_t order,bool_t useElementsOnFace) {
#endif
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements,NUMNODES,M0,M1,M2;
  index_t node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=2*NE0+1;
  N1=2*NE1+1;
  N2=2*NE2+1;

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
#ifdef PASO_MPI
  out=Finley_Mesh_alloc(name,3,order,mpi_info);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex20,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex20Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Hex20Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec8,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec8_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0*N1*N2);
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0*NDOF1*NDOF2, 0, 0 );
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
#else
  out=Finley_Mesh_alloc(name,3,order);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex20,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex20Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Hex20Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec8,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec8_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0*N1*N2);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

	/* create nodes */
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
#ifdef PASO_MPI
        out->Nodes->Dom[k]=NODE_INTERNAL;
#endif
      }
    }
  }

  
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
        node0=2*i0+2*i1*N0+2*N0*N1*i2;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);
#ifdef PASO_MPI
        out->Elements->Dom[k]=ELEMENT_INTERNAL;
#endif

        out->Elements->Nodes[INDEX2(0,k,20)]=node0;
        out->Elements->Nodes[INDEX2(1,k,20)]=node0+2;
        out->Elements->Nodes[INDEX2(2,k,20)]=node0+2*N0+2;
        out->Elements->Nodes[INDEX2(3,k,20)]=node0+2*N0;
        out->Elements->Nodes[INDEX2(4,k,20)]=node0+2*N0*N1;
        out->Elements->Nodes[INDEX2(5,k,20)]=node0+2*N0*N1+2;
        out->Elements->Nodes[INDEX2(6,k,20)]=node0+2*N0*N1+2*N0+2;
        out->Elements->Nodes[INDEX2(7,k,20)]=node0+2*N0*N1+2*N0;
        out->Elements->Nodes[INDEX2(8,k,20)]=node0+1;
        out->Elements->Nodes[INDEX2(9,k,20)]=node0+N0+2;
        out->Elements->Nodes[INDEX2(10,k,20)]=node0+2*N0+1;
        out->Elements->Nodes[INDEX2(11,k,20)]=node0+N0;
        out->Elements->Nodes[INDEX2(12,k,20)]=node0+N0*N1;
        out->Elements->Nodes[INDEX2(13,k,20)]=node0+N0*N1+2;
        out->Elements->Nodes[INDEX2(14,k,20)]=node0+N0*N1+2*N0+2;
        out->Elements->Nodes[INDEX2(15,k,20)]=node0+N0*N1+2*N0;
        out->Elements->Nodes[INDEX2(16,k,20)]=node0+2*N0*N1+1;
        out->Elements->Nodes[INDEX2(17,k,20)]=node0+2*N0*N1+N0+2;
        out->Elements->Nodes[INDEX2(18,k,20)]=node0+2*N0*N1+2*N0+1;
        out->Elements->Nodes[INDEX2(19,k,20)]=node0+2*N0*N1+N0;
      }
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);
  
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=20;
  } else {
     NUMNODES=8;
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
        node0=2*i0+2*i1*N0;
  
        out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
        out->FaceElements->Tag[k]=100;
        out->FaceElements->Color[k]=(i0%2)+2*(i1%2);
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
  
        if  (useElementsOnFace) {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0;
           out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+1;
           out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0*N1;
           out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0*N1+2;
           out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+2*N0*N1+N0;
           out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+2*N0*N1+1;
        } else {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
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
        node0=2*i0+2*i1*N0+2*N0*N1*(NE2-1);
  
        out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
        out->FaceElements->Tag[k]=200;
        out->FaceElements->Color[k]=(i0%2)+2*(i1%2)+4;
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

        if  (useElementsOnFace) {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2*N0;

           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0;

           out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+2*N0*N1+1;
           out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+2*N0*N1+N0;

           out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0*N1;
           out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0*N1+2;
           out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0*N1+2*N0;

           out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+1;
           out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0;

        } else {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0*N1+1;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1+N0;
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
         node0=2*i1*N0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=1;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+8;
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+2;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N0;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0*N1+1;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N0+1;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0+2;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
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
         node0=2*(NE0-1)+2*i1*N0+2*N0*N1*i2 ;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=2;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+12;
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N0*N1+2;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N0*N1+1;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0*N1;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1+2;
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
         node0=2*i0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Color[k]=(i2%2)+2*(i0%2)+16;
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N1*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N1*N0;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N1*N0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N1*N0+2*N0;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N1*N0;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N1*N0+N0+2;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N1*N0+N0;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N1*N0+2*N0;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N1*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N1*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*N0;
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
         node0=2*i0+2*(NE1-1)*N0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Color[k]=(i2%2)+2*(i0%2)+20;
#ifdef PASO_MPI
        out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N1*N0+2*N0;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+2*N0+1;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0+2;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N1*N0;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N1*N0+2*N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+1;
         }
       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=24;

#ifdef PASO_MPI
	Finley_ElementFile_setDomainFlags( out->Elements );
	Finley_ElementFile_setDomainFlags( out->FaceElements );
	Finley_ElementFile_setDomainFlags( out->ContactElements );
	Finley_ElementFile_setDomainFlags( out->Points );

	/* reorder the degrees of freedom */
	Finley_Mesh_resolveDegreeOfFreedomOrder( out, TRUE );
#endif
	
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
#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t* numElements,double* Length,bool_t* periodic,index_t order,bool_t useElementsOnFace) {
  dim_t N0,N0t,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF0t,NDOF1,NDOF2,NFaceElements,NUMNODES,M0,M1,M2;
  dim_t kk,iI, NE0_local, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, DOFBoundary[2];
	bool_t dom_left, dom_right, dom_internal;
  index_t firstNode=0, DOFcount=0, node0, node1, node2, idCount;
  index_t targetDomain=-1, firstNodeConstruct, j;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE, boundaryLeft=FALSE, boundaryRight=FALSE;
	index_t *indexBackward=NULL, *indexForward=NULL,*facePerm=NULL, *forwardDOF=NULL, *backwardDOF=NULL;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  Paso_MPIInfo *mpi_info = NULL;

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=2*NE0+1;
  N1=2*NE1+1;
  N2=2*NE2+1;

	index_t face0[] = {0,4,7,3,1,5,6,2,12,19,15,11,8,16,18,10,13,17,14,9}; 
	index_t face1[] = {1,2,6,5,0,3,7,4,9,14,17,13,8,10,18,16,11,15,19,12}; 
	index_t face2[] = {0,1,5,4,3,2,6,7,8,13,16,12,11,9,17,19,14,18,15,10}; 
	index_t face3[] = {3,7,6,2,0,4,5,1,15,18,14,10,11,19,17,9,12,16,13,8}; 
	index_t face4[] = {0,3,2,1,4,7,6,5,11,10,9,8,15,14,13,12,19,18,17,16};
	index_t face5[] = {4,5,6,7,0,1,2,3,16,17,18,19,12,13,14,15,8,9,10,11}; 
	
	index_t face0R[] = {0,4,7,3,12,19,15,11}; 
	index_t face1R[] = {1,2,6,5,9,14,17,13}; 
	index_t face2R[] = {0,1,5,4,8,13,16,12}; 
	index_t face3R[] = {3,7,6,2,15,18,14,10}; 
	index_t face4R[] = {0,3,2,1,11,10,9,8};
	index_t face5R[] = {4,5,6,7,16,17,18,19}; 
	
  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError())
        return NULL;

  /* use the serial version to generate the mesh for the 1-CPU case */
  if( mpi_info->size==1 )
  {
    out =  Finley_RectangularMesh_Hex20_singleCPU( numElements, Length, periodic, order, useElementsOnFace, mpi_info);
		return out;
  }    

  if( mpi_info->rank==0 )
    domLeft = TRUE;
  if( mpi_info->rank==mpi_info->size-1 )
    domRight = TRUE;
  if( mpi_info->rank>0 && mpi_info->rank<mpi_info->size-1 )
    domInternal = TRUE;

  /* dimensions of the local subdomain */
  domain_calculateDimension( mpi_info->rank, mpi_info->size, NE0, periodic[0], &numNodesLocal, &numDOFLocal, &numElementsLocal, &numElementsInternal, &firstNode, nodesExternal, DOFExternal, &numNodesExternal, periodicLocal, DOFBoundary );  

  NFaceElements=0;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements+=(domRight+domLeft)*NE1*NE2;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*numElementsLocal*NE2;
  } else {
      NDOF1=N1-1;
  }
  if (!periodic[2]) {
      NDOF2=N2;
      NFaceElements+=2*numElementsLocal*NE1;
  } else {
      NDOF2=N2-1;
  }
  
	boundaryLeft = !domLeft || periodicLocal[0];
	boundaryRight = !domRight || periodicLocal[1];
	N0t = numNodesLocal + boundaryRight + boundaryLeft*2;
	NDOF0t = numDOFLocal + boundaryRight + boundaryLeft*2;
	firstNodeConstruct = firstNode - 2*boundaryLeft;
	firstNodeConstruct = firstNodeConstruct<0 ? N0+firstNodeConstruct-1 : firstNodeConstruct;

  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d x %d mesh",N0,N1,N2);
  out=Finley_Mesh_alloc(name,3,order,mpi_info);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex20,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex20Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Hex20Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec8,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec8_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0t*N1*N2);
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, numDOFLocal*NDOF1*NDOF2, NDOF1*NDOF2*3, 0 );
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  #pragma omp parallel for private(i0,i1,i2,k)
  for (k=0,i2=0;i2<N2;i2++) {
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0t;i0++,k++) {         
        out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE((i0+firstNodeConstruct) % N0)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=i0 + (i1%NDOF1)*N0t + (i2%NDOF2)*N0t*N1;
		 		out->Nodes->Dom[k]=NODE_INTERNAL;
      }
    }
  }

	/* mark the nodes that reference external and boundary DOF as such */
	if( boundaryLeft ){
		for( i1=0; i1<N1; i1++ )
			for( i2=0; i2<N2; i2++ ) {
				out->Nodes->Dom[N1*N0t*i2+N0t*i1] = NODE_EXTERNAL;
				out->Nodes->Dom[N1*N0t*i2+N0t*i1+1] = NODE_EXTERNAL;
				out->Nodes->Dom[N1*N0t*i2+N0t*i1+2] = NODE_BOUNDARY; 
			}
	}
	if( boundaryRight ){
		for( i1=0; i1<N1; i1++ )
			for( i2=0; i2<N2; i2++ ) {
				out->Nodes->Dom[N1*N0t*i2+N0t*(i1+1)-1] = NODE_EXTERNAL;
				out->Nodes->Dom[N1*N0t*i2+N0t*(i1+1)-2] = NODE_BOUNDARY; 
				out->Nodes->Dom[N1*N0t*i2+N0t*(i1+1)-3] = NODE_BOUNDARY; 
			}
	}
	if( periodicLocal[0] ){
		for( i1=0; i1<N1; i1++ )
			for( i2=0; i2<N2; i2++ ) {
				out->Nodes->degreeOfFreedom[N1*N0t*i2+i1*N0t+3] = out->Nodes->degreeOfFreedom[i1*N0t+2];
				out->Nodes->Dom[N1*N0t*i2+N0t*i1+3] = NODE_BOUNDARY; 
			}
	}
		
  /* tag Nodes that are referenced by face elements */
  if (!periodic[2]) {    
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0t;i0++) {   
         out->Nodes->Tag[i0 + N0t*i1]+=100;
         out->Nodes->Tag[i0 + N0t*i1 + N0t*N1*(N2-1)]+=200;
       }
     }
  }
  if (!periodic[1]) {
    for (i2=0;i2<N2;i2++) {
      for (i0=0;i0<N0t;i0++) {
         out->Nodes->Tag[i0 + i2*N1*N0t]+=10;
         out->Nodes->Tag[i0 + (i2+1)*N1*N0t-N0t]+=20;
      }
    }
  }
  if (!periodic[0] && !domInternal ) {
    for (i2=0;i2<N2;i2++) {
      for (i1=0;i1<N1;i1++) {
        if( domLeft )
          out->Nodes->Tag[i1*N0t + i2*N0t*N1]+=1;
        if( domRight )
          out->Nodes->Tag[(i1+1)*N0t-1 + i2*N0t*N1]+=2;
      }
    }
  }

	/* form the boudary communication information */
	forwardDOF  = MEMALLOC(NDOF1*NDOF2*2,index_t);
	backwardDOF = MEMALLOC(NDOF1*NDOF2*2,index_t);
	if( !(mpi_info->size==2 && periodicLocal[0])){
		if( boundaryLeft  ) {
			targetDomain = mpi_info->rank-1 < 0 ? mpi_info->size-1 : mpi_info->rank-1;
			for( iI=0, i2=0; i2<NDOF2; i2++ ){
				for( i1=0; i1<NDOF1; i1++, iI+=2 ){
					forwardDOF[i1+i2*NDOF1]  = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t+2];
					backwardDOF[iI] = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t];
					backwardDOF[iI+1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t+1];
				}
			} 
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2*2, backwardDOF );
		}
		if( boundaryRight ) {
			targetDomain = mpi_info->rank+1 > mpi_info->size-1 ? 0 : mpi_info->rank+1;
			for( iI=0,i2=0; i2<NDOF2; i2++ ){
				for( i1=0; i1<NDOF1; i1++, iI+=2 ){
					forwardDOF[iI] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-3];
					forwardDOF[iI+1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-2];
					backwardDOF[i1+i2*NDOF1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-1];
				}
			} 
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2*2, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, backwardDOF );
		}
	} else{
		/* periodic boundary conditions with 2 domains, need to change the order in which domain 0 passes boundary data */
		targetDomain = 1;
		
		for( iI=0,i2=0; i2<NDOF2; i2++ ){
			for( i1=0; i1<NDOF1; i1++, iI+=2 ){
				forwardDOF[iI] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-3];
				forwardDOF[iI+1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-2];
				backwardDOF[i1+i2*NDOF1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+(i1+1)*N0t-1];
			}
		} 
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2*2, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, backwardDOF );

		for( iI=0, i2=0; i2<NDOF2; i2++ ){
			for( i1=0; i1<NDOF1; i1++, iI+=2 ){
				forwardDOF[i1+i2*NDOF1]  = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t+2];
				backwardDOF[iI] = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t];
				backwardDOF[iI+1] = out->Nodes->degreeOfFreedom[i2*N0t*N1+i1*N0t+1];
			}
		} 
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*NDOF2*2, backwardDOF );
	}
	MEMFREE( forwardDOF );
	MEMFREE( backwardDOF );
  /*   set the elements: */
  
	k=0;
  #pragma omp parallel for private(i0,i1,i2,k,node0) 
  for (i2=0;i2<NE2;i2++) {
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<numElementsLocal;i0++,k++) {
				node0 = (periodicLocal[0] && !i0) ? 2*(i1*N0t + i2*N1*N0t) : 2*(i1*N0t + i2*N1*N0t + i0) + periodicLocal[0];

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);
        out->Elements->Dom[k]=ELEMENT_INTERNAL;

        out->Elements->Nodes[INDEX2(0,k,20)]=node0;
        out->Elements->Nodes[INDEX2(1,k,20)]=node0+2;
        out->Elements->Nodes[INDEX2(2,k,20)]=node0+2*N0t+2;
        out->Elements->Nodes[INDEX2(3,k,20)]=node0+2*N0t;
        out->Elements->Nodes[INDEX2(4,k,20)]=node0+2*N0t*N1;
        out->Elements->Nodes[INDEX2(5,k,20)]=node0+2*N0t*N1+2;
        out->Elements->Nodes[INDEX2(6,k,20)]=node0+2*N0t*N1+2*N0t+2;
        out->Elements->Nodes[INDEX2(7,k,20)]=node0+2*N0t*N1+2*N0t;
        out->Elements->Nodes[INDEX2(8,k,20)]=node0+1;
        out->Elements->Nodes[INDEX2(9,k,20)]=node0+N0t+2;
        out->Elements->Nodes[INDEX2(10,k,20)]=node0+2*N0t+1;
        out->Elements->Nodes[INDEX2(11,k,20)]=node0+N0t;
        out->Elements->Nodes[INDEX2(12,k,20)]=node0+N0t*N1;
        out->Elements->Nodes[INDEX2(13,k,20)]=node0+N0t*N1+2;
        out->Elements->Nodes[INDEX2(14,k,20)]=node0+N0t*N1+2*N0t+2;
        out->Elements->Nodes[INDEX2(15,k,20)]=node0+N0t*N1+2*N0t;
        out->Elements->Nodes[INDEX2(16,k,20)]=node0+2*N0t*N1+1;
        out->Elements->Nodes[INDEX2(17,k,20)]=node0+2*N0t*N1+N0t+2;
        out->Elements->Nodes[INDEX2(18,k,20)]=node0+2*N0t*N1+2*N0t+1;
        out->Elements->Nodes[INDEX2(19,k,20)]=node0+2*N0t*N1+N0t;
      }
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);

	if( boundaryLeft ) 
		for( i2=0; i2<NE2; i2++ )
			for( i1=0; i1<NE1; i1++ )
				out->Elements->Dom[i2*NE1*numElementsLocal+i1*numElementsLocal]=ELEMENT_BOUNDARY;
	if( boundaryRight ) 
		for( i2=0; i2<NE2; i2++ )
			for( i1=0; i1<NE1; i1++ )
				out->Elements->Dom[i2*NE1*numElementsLocal+(i1+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;

	out->Elements->numElements = numElementsLocal*NE1*NE2;
	Finley_ElementFile_setDomainFlags( out->Elements );
  
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=20;
  } else {
     NUMNODES=8;
  }
  totalNECount=out->Elements->numElements;
  faceNECount=0;
	idCount = totalNECount;
  
  /*   these are the quadrilateral elements on boundary 1 (x3=0): */
	numElementsInternal = numElementsLocal-nodesExternal[0]-nodesExternal[1];
  if (!periodic[0] && !domInternal) {
     /* **  elements on boundary 001 (x1=0): */
  	 if( domLeft ){
			 #pragma omp parallel for private(i1,i2,k) 
			 for (i2=0;i2<NE2;i2++) {
				 for (i1=0;i1<NE1;i1++) {
					 k=i1+NE1*i2+faceNECount;
					 kk=i1*numElementsLocal + i2*numElementsLocal*NE1;
					 facePerm =!useElementsOnFace ? face0R :  face0;
		 
					 out->FaceElements->Id[k]=idCount++;
					 out->FaceElements->Tag[k]=1;
					 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
					 out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+8;

					 for( j=0; j<NUMNODES; j++ )
						out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
				 }
			 }
			 totalNECount+=NE1*NE2;
			 faceNECount+=NE1*NE2;
     }
     /* **  elements on boundary 002 (x1=1): */
 		 if( domRight ) { 
			 #pragma omp parallel for private(i1,i2,k) 
			 for (i2=0;i2<NE2;i2++) {
				 for (i1=0;i1<NE1;i1++) {
					 k=i1+NE1*i2+faceNECount;
					 kk=(i1+1)*numElementsLocal + i2*numElementsLocal*NE1 - 1;
					 facePerm =!useElementsOnFace ? face1R :  face1;
		 
					 out->FaceElements->Id[k]=idCount++;
					 out->FaceElements->Tag[k]=2;
         	 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
					 out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+12;

					 for( j=0; j<NUMNODES; j++ )
						out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
				 }
		 	 }
       totalNECount+=NE1*NE2;
       faceNECount+=NE1*NE2;
     }
  }
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,i2,k) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<numElementsLocal;i0++) {
         k=i0+numElementsLocal*i2+faceNECount;
				 kk=i0+numElementsLocal*NE1*i2;
				 facePerm =!useElementsOnFace ? face2R :  face2;
   
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+16;

				 for( j=0; j<NUMNODES; j++ )
					out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
       }
     }
 	   if( boundaryLeft ){
		 	for( i2=0; i2<NE2; i2++ )
		 		out->FaceElements->Dom[faceNECount+i2*numElementsLocal]=ELEMENT_BOUNDARY;
		 	if( periodicLocal[0] )
				for( i2=0; i2<NE2; i2++ )
					out->FaceElements->Dom[faceNECount+i2*numElementsLocal+1]=ELEMENT_BOUNDARY;
 		 }
 	   if( boundaryRight )
			for( i2=0; i2<NE2; i2++ )
				out->FaceElements->Dom[faceNECount+(i2+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;
     totalNECount+=numElementsLocal*NE2;
     faceNECount+=numElementsLocal*NE2;
  
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,i2,k) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<numElementsLocal;i0++) {
         k=i0+numElementsLocal*i2+faceNECount;
				 kk=i0+numElementsLocal*NE1*(i2+1)-numElementsLocal;
				 facePerm =!useElementsOnFace ? face3R :  face3;
   
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+20;

				 for( j=0; j<NUMNODES; j++ )
					out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
       }
     }
 	   if( boundaryLeft ){
		 	for( i2=0; i2<NE2; i2++ )
		 		out->FaceElements->Dom[faceNECount+i2*numElementsLocal]=ELEMENT_BOUNDARY;
		 	if( periodicLocal[0] )
				for( i2=0; i2<NE2; i2++ )
					out->FaceElements->Dom[faceNECount+i2*numElementsLocal+1]=ELEMENT_BOUNDARY;
 		 }
 	   if( boundaryRight )
			for( i2=0; i2<NE2; i2++ )
				out->FaceElements->Dom[faceNECount+(i2+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;
     totalNECount+=numElementsLocal*NE2;
     faceNECount+=numElementsLocal*NE2;
  }
  if (!periodic[2]) {
     /*  elements on boundary 100 (x3=0): */
  
     #pragma omp parallel for private(i0,i1,k) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0; i0<numElementsLocal; i0++) {
         k=i0+numElementsLocal*i1+faceNECount;
				 kk=i0 + i1*numElementsLocal;
				 facePerm =!useElementsOnFace ? face4R :  face4;
   
         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=100;
         out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2);

				 for( j=0; j<NUMNODES; j++ )
					out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
       }
     }
 	   if( boundaryLeft ){
		 	for( i1=0; i1<NE1; i1++ )
		 		out->FaceElements->Dom[faceNECount+i1*numElementsLocal]=ELEMENT_BOUNDARY;
		 	if( periodicLocal[0] )
				for( i1=0; i1<NE1; i1++ )
					out->FaceElements->Dom[faceNECount+i1*numElementsLocal+1]=ELEMENT_BOUNDARY;
 		 }
 	   if( boundaryRight )
			for( i1=0; i1<NE1; i1++ )
				out->FaceElements->Dom[faceNECount+(i1+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;
     totalNECount+=NE1*numElementsLocal;
     faceNECount+=NE1*numElementsLocal;
 		 
     /* **  elements on boundary 200 (x3=1) */
  
     #pragma omp parallel for private(i0,i1,k) 
     for (i1=0;i1<NE1;i1++) {
       for (i0=0;i0<numElementsLocal;i0++) {
         k=i0+numElementsLocal*i1+faceNECount;
				 kk=i0+i1*numElementsLocal+numElementsLocal*NE1*(NE2-1);
				 facePerm = !useElementsOnFace ? face5R : face5;

         out->FaceElements->Id[k]=idCount++;
         out->FaceElements->Tag[k]=200;
         out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2)+4;

				 for( j=0; j<NUMNODES; j++ )
					out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,20)];
       }
     }
 	   if( boundaryLeft ){
		 	for( i1=0; i1<NE1; i1++ )
		 		out->FaceElements->Dom[faceNECount+i1*numElementsLocal]=ELEMENT_BOUNDARY;
		 	if( periodicLocal[0] )
				for( i1=0; i1<NE1; i1++ )
					out->FaceElements->Dom[faceNECount+i1*numElementsLocal+1]=ELEMENT_BOUNDARY;
 		 }
 	   if( boundaryRight )
			for( i1=0; i1<NE1; i1++ )
				out->FaceElements->Dom[faceNECount+(i1+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;
     totalNECount+=NE1*numElementsLocal;
     faceNECount+=NE1*numElementsLocal;
  }
	out->FaceElements->elementDistribution->numInternal = faceNECount;
	
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=23;
	out->FaceElements->numElements=faceNECount;
	Finley_ElementFile_setDomainFlags( out->FaceElements );

  /* setup distribution info for other elements */
	Finley_ElementFile_setDomainFlags( out->ContactElements );
	Finley_ElementFile_setDomainFlags( out->Points );

	/* reorder the degrees of freedom */
	Finley_Mesh_resolveDegreeOfFreedomOrder( out, TRUE );

  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds(out);
  if( !Finley_MPI_noError(mpi_info) )
  {
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  } 

  /* prepare mesh for further calculatuions:*/
  Finley_Mesh_prepare(out);
  if( !Finley_MPI_noError(mpi_info) )
  {
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  } 

  /* free up memory */
  Paso_MPIInfo_dealloc( mpi_info );

	//print_mesh_statistics( out, FALSE );

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  return out; 
}
#endif

