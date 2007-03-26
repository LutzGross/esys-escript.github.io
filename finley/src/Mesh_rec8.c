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

/*   Generates a numElements[0] x numElements[1] mesh with second order elements (Rec8) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

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
Finley_Mesh* Finley_RectangularMesh_Rec8_singleCPU(int* numElements,double* Length,int* periodic,int order, index_t reduced_order, int useElementsOnFace, Paso_MPIInfo *mpi_info) 
#else
Finley_Mesh* Finley_RectangularMesh_Rec8(int* numElements,double* Length,int* periodic,int order, index_t reduced_order, int useElementsOnFace) 
#endif
{
  dim_t N0,N1,NE0,NE1,i0,i1,totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NUMNODES,M0,M1;
  index_t k,node0;
  Finley_Mesh* out=NULL;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=2*NE0+1;
  N1=2*NE1+1;

#ifndef PASO_MPI
  if (N0<=N1) {
     M0=1;
     M1=N0;
  } else {
     M0=N1;
     M1=1;
  }
#else
     M0=1;
     M1=N0;
#endif

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
#ifdef PASO_MPI
  out=Finley_Mesh_alloc(name,2,order, reduced_order,mpi_info);
#else
  out=Finley_Mesh_alloc(name,2,order, reduced_order);
#endif
  if (! Finley_noError()) return NULL;

#ifdef PASO_MPI
  out->Elements=Finley_ElementFile_alloc(Rec8,out->order, out->reduced_order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec8Face,out->order, out->reduced_order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec8Face_Contact,out->order, out->reduced_order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line3,out->order, out->reduced_order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line3_Contact,out->order, out->reduced_order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order, out->reduced_order,mpi_info);
#else
  out->Elements=Finley_ElementFile_alloc(Rec8,out->order, out->reduced_order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec8Face,out->order, out->reduced_order);
     out->ContactElements=Finley_ElementFile_alloc(Rec8Face_Contact,out->order, out->reduced_order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line3,out->order, out->reduced_order);
     out->ContactElements=Finley_ElementFile_alloc(Line3_Contact,out->order, out->reduced_order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order, out->reduced_order);
#endif

  if (! Finley_noError()) {
       Finley_Mesh_dealloc(out);
       return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0*N1);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
#ifdef PASO_MPI
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0*NDOF1, 0, 0 );
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, NE0*NE1, NE0*NE1);
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, NFaceElements );
#endif
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
#ifdef PASO_MPI
      out->Nodes->Dom[k]=NODE_INTERNAL;
#endif
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
      node0=2*i0+2*i1*N0;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);
#ifdef PASO_MPI
      out->Elements->Dom[k]=ELEMENT_INTERNAL;
#endif

      out->Elements->Nodes[INDEX2(0,k,8)]=node0;
      out->Elements->Nodes[INDEX2(1,k,8)]=node0+2;
      out->Elements->Nodes[INDEX2(2,k,8)]=node0+2*N0+2;
      out->Elements->Nodes[INDEX2(3,k,8)]=node0+2*N0;
      out->Elements->Nodes[INDEX2(4,k,8)]=node0+1;
      out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0+2;
      out->Elements->Nodes[INDEX2(6,k,8)]=node0+2*N0+1;
      out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0;

    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0);
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=3;
  }
  
  totalNECount=NE0*NE1;
  faceNECount=0;
  
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
       k=i0+faceNECount;
       node0=2*i0;

       out->FaceElements->Id[k]=i0+totalNECount;
       out->FaceElements->Tag[k]=10;
       out->FaceElements->Color[k]=i0%2;
#ifdef PASO_MPI
	     out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
       }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
  
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
       k=i0+faceNECount;
       node0=2*i0+2*(NE1-1)*N0;

       out->FaceElements->Id[k]=i0+totalNECount;
       out->FaceElements->Tag[k]=20;
       out->FaceElements->Color[k]=i0%2+2;
#ifdef PASO_MPI
	     out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0+2;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+1;
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
       node0=2*i1*N0;

       out->FaceElements->Id[k]=i1+totalNECount;
       out->FaceElements->Tag[k]=1;
       out->FaceElements->Color[k]=(i1%2)+4;
#ifdef PASO_MPI
	     out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+1;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0;
       }
   
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  
     /* **  elements on boundary 002 (x1=1): */
     
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       k=i1+faceNECount;
       node0=2*(NE0-1)+2*i1*N0;
   
       out->FaceElements->Id[k]=i1+totalNECount;
       out->FaceElements->Tag[k]=2;
       out->FaceElements->Color[k]=(i1%2)+6;
#ifdef PASO_MPI
	     out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
   
       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+2;
       }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;

  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=7;

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
	
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  /* prepare mesh for further calculatuions:*/
  Finley_Mesh_prepare(out) ;


  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  if (Finley_noError()) {
       if ( ! Finley_Mesh_isPrepared(out)) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  } else {
      Finley_Mesh_dealloc(out);
  }
  return out;
}

#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Rec8(int* numElements,double* Length,int* periodic,int order, index_t reduced_order ,int useElementsOnFace) {
  dim_t N0,N1,NE0,NE1,i0,i1,totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NUMNODES, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, DOFBoundary[2];
	dim_t N0t, NDOF0t;
  index_t *numForward=NULL, *numBackward=NULL;
  index_t kk, k,i,j,node0,firstNodeConstruct,firstNode=0, DOFcount=0, forwardDOF=NULL, backwardDOF=NULL, targetDomain=-1,faceNEBoundary;
  Finley_Mesh* out=NULL;
  char name[50];
	index_t *indexForward=NULL, *indexBackward=NULL;
  double time0=Finley_timer();
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE, boundaryLeft=FALSE, boundaryRight=FALSE;
  Paso_MPIInfo *mpi_info = NULL;
	index_t faceStart, numFaceLeft;

	/* these are used in constructing the face elements, and give the permutation
	   of the node order of a normal element for each face */
	index_t *facePerm;	 
	index_t face0[]={0, 1, 2, 3, 4, 5, 6, 7};
	index_t face1[]={2, 3, 0, 1, 6, 7, 4, 5};
	index_t face2[]={3, 0, 1, 2, 7, 4, 5, 6};
	index_t face3[]={1, 2, 3, 0, 5, 6, 7, 4};
	index_t face0nc[]={0, 1, 4};
	index_t face1nc[]={2, 3, 6};
	index_t face2nc[]={3, 0, 7};
	index_t face3nc[]={1, 2, 5};

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=2*NE0+1;
  N1=2*NE1+1;

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError()) {
        Finley_Mesh_dealloc(out);
        return NULL;
  }

	// use the serial method for the single CPU case
	if( mpi_info->size==1 ){
		out = Finley_RectangularMesh_Rec8_singleCPU( numElements,Length,periodic,order, reduced_order,useElementsOnFace, mpi_info);
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
      NFaceElements+=(domLeft + domRight)*NE1;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*numElementsLocal;
  } else {
      NDOF1=N1-1;
  }

	boundaryLeft = !domLeft || periodicLocal[0];
	boundaryRight = !domRight || periodicLocal[1];
	N0t = numNodesLocal + boundaryRight + boundaryLeft*2;
	NDOF0t = numDOFLocal + boundaryRight + boundaryLeft*2;
	firstNodeConstruct = firstNode - 2*boundaryLeft;
	firstNodeConstruct = firstNodeConstruct<0 ? N0+firstNodeConstruct-1 : firstNodeConstruct;

  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d mesh",N0,N1);
  out=Finley_Mesh_alloc(name,2,order,reduced_order,mpi_info);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Rec8,out->order, out->reduced_order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec8Face,out->order, out->reduced_order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec8Face_Contact,out->order, out->reduced_order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line3,out->order, out->reduced_order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line3_Contact,out->order, out->reduced_order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order, out->reduced_order,mpi_info);

  if (! Finley_noError()) {
       Finley_Mesh_dealloc(out);
       return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0t*N1);
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
	
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0t*NDOF1, NDOF1*3, 0 );
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, numElementsLocal*NE1, NE1*(numElementsLocal-boundaryRight*(!periodic[1])) );
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, NFaceElements-2*boundaryRight*(!periodic[1]) );
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /*  set nodes: */
  #pragma omp parallel for private(i0,i1,k)
  for (i1=0;i1<N1;i1++) {
    for (i0=0;i0<N0t;i0++, k++) {
			k=i0+i1*N0t;
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE((firstNodeConstruct + i0) % N0)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=i0 + i1*N0t;
      out->Nodes->Tag[k]=0;
      out->Nodes->Dom[k]=NODE_INTERNAL;
      out->Nodes->degreeOfFreedom[k]=i0 + (i1%NDOF1)*N0t;
    }
		/* take care to make the DOF reflect the internal/external position of the DOF */
		if( boundaryLeft ){
			out->Nodes->Dom[i1*N0t] = NODE_EXTERNAL;	
			out->Nodes->Dom[i1*N0t+1] = NODE_EXTERNAL;	
			out->Nodes->Dom[i1*N0t+2] = NODE_BOUNDARY;	
			if( periodicLocal[0] ){
				out->Nodes->degreeOfFreedom[i1*N0t+3] = out->Nodes->degreeOfFreedom[i1*N0t+2];
				out->Nodes->Dom[i1*N0t+3] = NODE_BOUNDARY;	
			}
		}
		if( boundaryRight ){
			out->Nodes->Dom[(i1+1)*N0t-1] = NODE_EXTERNAL;	
			out->Nodes->Dom[(i1+1)*N0t-2] = NODE_BOUNDARY;	
			out->Nodes->Dom[(i1+1)*N0t-3] = NODE_BOUNDARY;	
		}
  }
  /* tags for the faces: */
  if (!periodic[1]) {
    for (i0=0;i0<N0t;i0++) {
      out->Nodes->Tag[i0]+=10;
      out->Nodes->Tag[i0+N0t*(N1-1)]+=20;
    }
  }
  if (domLeft && !periodic[0]) {
    for (i1=0;i1<N1;i1++) {
      out->Nodes->Tag[i1*N0t]+=1;
    }
  }
	if (domRight && !periodic[0]){
    for (i1=0;i1<N1;i1++) {
      out->Nodes->Tag[(i1+1)*N0t-1]+=2;
		}
	}
  
	/* setup the receive information */
	indexBackward = TMPMEMALLOC( NDOF1*2, index_t );
	indexForward  = TMPMEMALLOC( NDOF1*2, index_t );
	if( !(mpi_info->size==2 && periodicLocal[0])){
		/* Backward (receive) DOF */
		if( boundaryLeft ){
			targetDomain = mpi_info->rank-1<0 ? mpi_info->size-1 : mpi_info->rank-1;

			for( i=0, i1=0; i<NDOF1; i++, i1+=2 ){
				indexBackward[i1] = out->Nodes->degreeOfFreedom[i*N0t];
				indexBackward[i1+1] = out->Nodes->degreeOfFreedom[i*N0t+1];
				indexForward[i] = out->Nodes->degreeOfFreedom[i*N0t+2];
			}

			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexForward );
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*2, indexBackward );
		}
		
		if( boundaryRight ){
			targetDomain = mpi_info->rank+1 > mpi_info->size-1 ? 0 : mpi_info->rank+1;

			for( i=0, i1=0; i<NDOF1; i++, i1+=2 ){
				indexBackward[i] = out->Nodes->degreeOfFreedom[(i+1)*N0t-1];
				indexForward[i1]   = out->Nodes->degreeOfFreedom[(i+1)*N0t-3];
				indexForward[i1+1] = out->Nodes->degreeOfFreedom[(i+1)*N0t-2];
			}

			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*2, indexForward );
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexBackward );
		}
	}
	else{
		targetDomain = 1;

		for( i=0, i1=0; i<NDOF1; i++, i1+=2 ){
				indexBackward[i] = out->Nodes->degreeOfFreedom[(i+1)*N0t-1];
				indexForward[i1]   = out->Nodes->degreeOfFreedom[(i+1)*N0t-3];
				indexForward[i1+1] = out->Nodes->degreeOfFreedom[(i+1)*N0t-2];
		}

		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*2, indexForward );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexBackward );

		for( i=0, i1=0; i<NDOF1; i++, i1+=2 ){
				indexBackward[i1] = out->Nodes->degreeOfFreedom[i*N0t];
				indexBackward[i1+1] = out->Nodes->degreeOfFreedom[i*N0t+1];
				indexForward[i] = out->Nodes->degreeOfFreedom[i*N0t+2];
		}

		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, indexForward );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1*2, indexBackward );
	}
	TMPMEMFREE( indexBackward );
	TMPMEMFREE( indexForward );

  /*   set the elements: */
  #pragma omp parallel for private(i0,i1,k,node0) 
  for (i1=0;i1<NE1;i1++) {
    for (i0=0;i0<numElementsLocal;i0++) {
      k=i0+numElementsLocal*i1;
			node0 = (periodicLocal[0] && !i0) ? 2*i1*N0t :  2*i0+2*i1*N0t+periodicLocal[0] ;

      out->Elements->Id[k]=((firstNodeConstruct/2+i0)%NE0) * NE1 + i1;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);
			out->Elements->Dom[k] = ELEMENT_INTERNAL;

      out->Elements->Nodes[INDEX2(0,k,8)]=node0;
      out->Elements->Nodes[INDEX2(1,k,8)]=node0+2;
      out->Elements->Nodes[INDEX2(2,k,8)]=node0+2*N0t+2;
      out->Elements->Nodes[INDEX2(3,k,8)]=node0+2*N0t;
      out->Elements->Nodes[INDEX2(4,k,8)]=node0+1;
      out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0t+2;
      out->Elements->Nodes[INDEX2(6,k,8)]=node0+2*N0t+1;
      out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0t;
    }
		if( boundaryLeft )
			out->Elements->Dom[i1*numElementsLocal] = ELEMENT_BOUNDARY;
		if( boundaryRight )
			out->Elements->Dom[(i1+1)*numElementsLocal-1] = ELEMENT_BOUNDARY;
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0);
  
	/* set the distribution information for the elements */
	Finley_ElementFile_setDomainFlags( out->Elements );

  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=3;
  }
  
	Finley_ElementFile_setDomainFlags( out->FaceElements );
  totalNECount=numElementsLocal*NE1;
  faceNECount=0;
	faceNEBoundary = totalNECount + numElementsInternal*(!periodic[1])*2 + NE1*(!periodic[0])*(domLeft+domRight);
  
	numFaceLeft = domLeft*(!periodic[0])*NE1;
	faceStart = out->FaceElements->elementDistribution->vtxdist[mpi_info->rank];
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x1=0): */
  
     #pragma omp parallel for private(i0,k,kk,facePerm) 
     for (i0=0;i0<numElementsLocal; i0++) {
       k=i0+faceNECount;
			 kk=i0;
			 facePerm = useElementsOnFace ? face0 : face0nc;

       out->FaceElements->Id[k]=faceStart + numFaceLeft + i0*2;
       out->FaceElements->Tag[k]=10;
       out->FaceElements->Color[k]=i0%2;
			 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
				
			 for( j=0; j<NUMNODES; j++ )
				  out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,8)];
     }
		 if( boundaryLeft ){
			 out->FaceElements->Dom[faceNECount]=ELEMENT_BOUNDARY;
			 if( periodicLocal[0] )
			 	 out->FaceElements->Dom[faceNECount+1]=ELEMENT_BOUNDARY;
		 }
		 if( boundaryRight ){
			 out->FaceElements->Dom[faceNECount+numElementsLocal-1]=ELEMENT_BOUNDARY;
			 out->FaceElements->Id[faceNECount+numElementsLocal-1]=out->FaceElements->elementDistribution->vtxdist[mpi_info->rank+1];
		 }
     totalNECount+=numElementsLocal;
     faceNECount+=numElementsLocal;
  
     /* **  elements on boundary 020 (x1=1): */
  
     #pragma omp parallel for private(i0,k,kk,facePerm) 
     for (i0=0;i0<numElementsLocal;i0++) {
       k=i0+faceNECount;
			 kk=i0 + numElementsLocal*(NE1-1);
			 facePerm = useElementsOnFace ? face1 : face1nc;

       out->FaceElements->Id[k]=faceStart + numFaceLeft + i0*2+1;
       out->FaceElements->Tag[k]=20;
       out->FaceElements->Color[k]=i0%2+2;
			 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;

			 for( j=0; j<NUMNODES; j++ )
				  out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,8)];
     }
		 if( boundaryLeft ){
			 out->FaceElements->Dom[faceNECount]=ELEMENT_BOUNDARY;
			 if( periodicLocal[0] )
			 	 out->FaceElements->Dom[faceNECount+1]=ELEMENT_BOUNDARY;
		 }
		 if( boundaryRight ){
			 out->FaceElements->Dom[faceNECount+numElementsLocal-1]=ELEMENT_BOUNDARY;
			 out->FaceElements->Id[faceNECount+numElementsLocal-1]=out->FaceElements->elementDistribution->vtxdist[mpi_info->rank+1]+1;
		 }
     totalNECount+=numElementsLocal;
     faceNECount+=numElementsLocal;
  }
  if (domLeft && !periodicLocal[0]) {
     /* **  elements on boundary 001 (x0=0): */
  
     #pragma omp parallel for private(i1,k,kk,facePerm) 
     for (i1=0;i1<NE1;i1++) {
       k=i1+faceNECount;
			 kk=i1*numElementsLocal;
			 facePerm = useElementsOnFace ? face2 : face2nc;

       out->FaceElements->Id[k]=faceStart + i1;
       out->FaceElements->Tag[k]=1;
       out->FaceElements->Color[k]=(i1%2)+4;
			 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;

			 for( j=0; j<NUMNODES; j++ )
				  out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,8)];
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  
     /* **  elements on boundary 002 (x0=1): */
  } 
  if (domRight && !periodicLocal[1]) {
     #pragma omp parallel for private(i1,k,kk,facePerm) 
     for (i1=0;i1<NE1;i1++) {
       k=i1+faceNECount;
			 kk=(i1+1)*numElementsLocal-1;
			 facePerm = useElementsOnFace ? face3 : face3nc;
   
       out->FaceElements->Id[k]=faceStart + numElementsLocal*2*(!periodic[1])+i1;
       out->FaceElements->Tag[k]=2;
       out->FaceElements->Color[k]=(i1%2)+6;
			 out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
   
			 for( j=0; j<NUMNODES; j++ )
				  out->FaceElements->Nodes[INDEX2(j,k,NUMNODES)]=out->Elements->Nodes[INDEX2(facePerm[j],kk,8)];
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=7;
	out->FaceElements->numElements = faceNECount;

  /* setup distribution info for other elements */
	Finley_ElementFile_setDomainFlags( out->ContactElements );
	Finley_ElementFile_setDomainFlags( out->Points );

	Finley_Mesh_prepareElementDistribution( out );

	/* reorder the degrees of freedom */
	Finley_Mesh_resolveDegreeOfFreedomOrder( out, TRUE );
	
  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds(out);
  if ( !Finley_MPI_noError( mpi_info )) {
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  }

  /* prepare mesh for further calculatuions:*/
  Finley_Mesh_prepare(out) ;

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  if ( !Finley_MPI_noError( mpi_info )) {
    Paso_MPIInfo_dealloc( mpi_info );
    Finley_Mesh_dealloc(out);
    return NULL;
  }
  if (Finley_noError()) {
       if (! Finley_Mesh_isPrepared(out)) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  }
  return out;
}
#endif
