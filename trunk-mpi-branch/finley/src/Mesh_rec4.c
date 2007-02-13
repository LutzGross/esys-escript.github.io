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
static void domain_calculateDimension( index_t rank, dim_t size, dim_t numElementsGlobal, bool_t periodic, dim_t *numNodesLocal, dim_t *numDOFLocal, dim_t *numElementsLocal, dim_t *numElementsInternal, dim_t *firstNode, dim_t *nodesExternal, dim_t *DOFExternal, dim_t *numNodesExternal, bool_t *periodicLocal, dim_t *numDOFInternal )
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
    numDOFLocal[0]--;
		
	numDOFInternal[0] = numDOFLocal[0];
	if( size>1 )
	{	
		if( rank>0 || periodic ){
			numDOFInternal[0] -= 1;
		}	
		if( rank<(size-1) || periodic ){
			numDOFInternal[0] -= 1;
		}	
	}	
  DOFExternal[0] = nodesExternal[0];
  DOFExternal[1] = nodesExternal[1];
}
#endif

#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Rec4_singleCPU(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace,Paso_MPIInfo *mpi_info) 
#else
Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) 
#endif
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
#ifdef PASO_MPI
  out=Finley_Mesh_alloc(name,2,order,mpi_info);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Rec4,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec4Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Rec4Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line2,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line2_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0*N1);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0*NDOF1, 0, 0);
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, NE0*NE1, NE0*NE1);
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, NFaceElements );
#else
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
      node0=i0+i1*N0;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);
#ifdef PASO_MPI
      out->Elements->Dom[k]=ELEMENT_INTERNAL;
#endif

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
#ifdef PASO_MPI
		      out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
   
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
#ifdef PASO_MPI
		      out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
   
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
#ifdef PASO_MPI
		      out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif

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
#ifdef PASO_MPI
		      out->FaceElements->Dom[k]=ELEMENT_INTERNAL;
#endif
   
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
  if (Finley_noError()) {
       if ( ! Finley_Mesh_isPrepared(out) ) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  } else {
      Finley_Mesh_dealloc(out);
  }
  return out;
}
#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) 
{
  dim_t N0,N1,NE0,NE1,i0,i1,N0t, NDOF0t, NE0_local, NDOF0,NDOF1, NFaceElements, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, faceNECount, totalNECount,M0,M1,numDOFInternal;
  index_t NUMNODES,k,firstNode=0, DOFcount=0, node0, node1;
  index_t targetDomain=-1, firstNodeConstruct;
  index_t *forwardDOF=NULL, *backwardDOF=NULL;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE, boundaryLeft=FALSE, boundaryRight=FALSE;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
	index_t faceStart, numFaceLeft, i;

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=NE0+1;
  N1=NE1+1;
  Paso_MPIInfo *mpi_info = NULL;

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError())
        return NULL;

	/* use the serial code to generate the mesh in the 1-CPU case */
	if( mpi_info->size==1 ){
		out = Finley_RectangularMesh_Rec4_singleCPU(numElements,Length,periodic,order,useElementsOnFace,mpi_info);
		return out;
	}

  if( mpi_info->rank==0 )
    domLeft = TRUE;
  if( mpi_info->rank==mpi_info->size-1 )
    domRight = TRUE;
  if( mpi_info->rank>0 && mpi_info->rank<mpi_info->size-1 )
    domInternal = TRUE;

  /* dimensions of the local subdomain */
  domain_calculateDimension( mpi_info->rank, mpi_info->size, NE0, periodic[0], &numNodesLocal, &numDOFLocal, &numElementsLocal, &numElementsInternal, &firstNode, nodesExternal, DOFExternal, &numNodesExternal, periodicLocal, &numDOFInternal );  

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
      NFaceElements+=2*numElementsLocal;
  } else {
      NDOF1=N1-1;
  }
  
	boundaryLeft = !domLeft || periodicLocal[0];
	boundaryRight = !domRight || periodicLocal[1];
	N0t = numNodesLocal + boundaryRight + boundaryLeft;
	NDOF0t = numDOFLocal + boundaryRight + boundaryLeft;
	firstNodeConstruct = firstNode - boundaryLeft;
	firstNodeConstruct = firstNodeConstruct<0 ? N0-2 : firstNodeConstruct;

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
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, NDOF0t*NDOF1, 2*NDOF1, 0 );
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, numElementsLocal*NE1, NE1*(numElementsLocal-boundaryRight*(!periodic[1])) );
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, NFaceElements-2*boundaryRight*(!periodic[1]) );

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /*  set nodes: */
  #pragma omp parallel for private(i0,i1,k)
  for (i1=0;i1<N1;i1++) {
    for ( i0=0;i0<N0t;i0++,k++) {
			k=i0+i1*N0t;
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE((i0+firstNodeConstruct) % N0)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=k;
      out->Nodes->Tag[k]=0;
      out->Nodes->degreeOfFreedom[k]=i0 + (i1%NDOF1)*N0t;
		 	out->Nodes->Dom[k]=NODE_INTERNAL;
    }
  }
	if( boundaryLeft ){
		for( i1=0; i1<N1; i1++ ){
			out->Nodes->Dom[N0t*i1] = NODE_EXTERNAL;
			out->Nodes->Dom[N0t*i1+1] = NODE_BOUNDARY; 
		}
	}
	if( boundaryRight ){
		for( i1=0; i1<N1; i1++ ){
			out->Nodes->Dom[N0t*(i1+1)-1] = NODE_EXTERNAL;
			out->Nodes->Dom[N0t*(i1+1)-2] = NODE_BOUNDARY; 
		}
	}
	if( periodicLocal[0] ){
		for( i1=0; i1<N1; i1++ ){
			out->Nodes->degreeOfFreedom[i1*N0t+2] = out->Nodes->degreeOfFreedom[i1*N0t+1];
			out->Nodes->Dom[N0t*i1+2] = NODE_BOUNDARY; 
		}
	}
		
  /* tag 'em */
  for(i0=0; i0<N0t; i0++ )
  {
    out->Nodes->Tag[i0] += 10;                        // bottom
    out->Nodes->Tag[i0+N0t*(N1-1)] += 20;   // top
  }
  if( domLeft && !periodicLocal[0] )
    for( i1=0; i1<N1; i1++ ) 
      out->Nodes->Tag[i1*N0t] += 1;         // left

  if( domRight && !periodicLocal[0] )
    for( i1=0; i1<N1; i1++ )
      out->Nodes->Tag[(i1+1)*N0t-1] += 2;   // right

	/* form the boudary communication information */
	forwardDOF  = MEMALLOC(NDOF1,index_t);
	backwardDOF = MEMALLOC(NDOF1,index_t);
	if( !(mpi_info->size==2 && periodicLocal[0])){
		if( boundaryLeft  ) {
			targetDomain = mpi_info->rank-1 < 0 ? mpi_info->size-1 : mpi_info->rank-1;
			for( i1=0; i1<NDOF1; i1++ ){
				forwardDOF[i1] = out->Nodes->degreeOfFreedom[i1*N0t+1];
				backwardDOF[i1] = out->Nodes->degreeOfFreedom[i1*N0t];
			} 
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, backwardDOF );
		}
		if( boundaryRight ) {
			targetDomain = mpi_info->rank+1 > mpi_info->size-1 ? 0 : mpi_info->rank+1;
			for( i1=0; i1<NDOF1; i1++ ){
				forwardDOF[i1] = out->Nodes->degreeOfFreedom[(i1+1)*N0t-2];
				backwardDOF[i1] = out->Nodes->degreeOfFreedom[(i1+1)*N0t-1];
			} 
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, backwardDOF );
		}
	} else{
		/* periodic boundary conditions with 2 domains, need to change the order in which domain 0 passes boundary data */
		targetDomain = 1;
		
		for( i1=0; i1<NDOF1; i1++ ){
			forwardDOF[i1] = out->Nodes->degreeOfFreedom[(i1+1)*N0t-2];
			backwardDOF[i1] = out->Nodes->degreeOfFreedom[(i1+1)*N0t-1];
		} 
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, backwardDOF );

		for( i1=0; i1<NDOF1; i1++ ){
			forwardDOF[i1] = out->Nodes->degreeOfFreedom[i1*N0t+1];
			backwardDOF[i1] = out->Nodes->degreeOfFreedom[i1*N0t];
		} 
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, NDOF1, backwardDOF );
	}
	MEMFREE( forwardDOF );
	MEMFREE( backwardDOF );

  /* elements: */
  
 #pragma omp parallel for private(i0,i1,k,node0)  
 for ( i1=0; i1<NE1; i1++) {
    for ( i0=0; i0<numElementsLocal; i0++, k++) {
			k=i1*numElementsLocal + i0;
			node0 = (periodicLocal[0] && !i0) ? i1*N0t :  i1*N0t + i0 + periodicLocal[0];

      out->Elements->Id[k]=((firstNodeConstruct+i0)%NE0) * NE1 + i1;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);
			out->Elements->Dom[k]=ELEMENT_INTERNAL;

      out->Elements->Nodes[INDEX2(0,k,4)]=node0;
      out->Elements->Nodes[INDEX2(1,k,4)]=node0+1;
      out->Elements->Nodes[INDEX2(2,k,4)]=node0+N0t+1;
      out->Elements->Nodes[INDEX2(3,k,4)]=node0+N0t;
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0);
	if( boundaryLeft ) 
		for( i1=0; i1<NE1; i1++ )
			out->Elements->Dom[i1*numElementsLocal]=ELEMENT_BOUNDARY;
	if( boundaryRight ) 
		for( i1=0; i1<NE1; i1++ )
			out->Elements->Dom[(i1+1)*numElementsLocal-1]=ELEMENT_BOUNDARY;
		
	Finley_ElementFile_setDomainFlags( out->Elements );

  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=4;
  } else {
     NUMNODES=2;
  }
  totalNECount=numElementsLocal*NE1;
  faceNECount=0;
  
	numFaceLeft = domLeft*(!periodic[0])*NE1;
	faceStart = out->FaceElements->elementDistribution->vtxdist[mpi_info->rank];
	if( !periodic[1] ){
		
	  /* elements on boundary 010 (x1=0): */
	  #pragma omp parallel for private(i0,k,node0) 
	  for (i0=0;i0<numElementsLocal;i0++) {
	 			k=i0+faceNECount;
				node0 = (periodicLocal[0] && !i0) ? 0 :  i0 + periodicLocal[0];
 
				out->FaceElements->Id[k]=faceStart + numFaceLeft + i0*2;
				out->FaceElements->Tag[k]   = 10;
				out->FaceElements->Color[k] = i0%2;
				out->FaceElements->Dom[k] = ELEMENT_INTERNAL;

				if (useElementsOnFace) {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
						out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0t+1;
						out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0t;
				} else {
						out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
						out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
				}
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

	  totalNECount += numElementsLocal;
	  faceNECount  += numElementsLocal;
	 
	  /* **  elements on boundary 020 (x1=1): */

	  #pragma omp parallel for private(i0,k,node0) 
	  for (i0=0;i0<numElementsLocal;i0++) {
			k=i0+faceNECount;
			node0 = (periodicLocal[0] && !i0) ? (NE1-1)*N0t :  i0 + periodicLocal[0] + (NE1-1)*N0t;

      out->FaceElements->Id[k]=faceStart + numFaceLeft + i0*2+1;
			out->FaceElements->Tag[k]   = 20;
			out->FaceElements->Color[k] = i0%2+2;
			out->FaceElements->Dom[k] = ELEMENT_INTERNAL;

			if (useElementsOnFace) {
				out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0t+1;
				out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0t;
				out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
				out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
			} else {
				out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0t+1;
				out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0t;
			}
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

	  totalNECount += numElementsLocal;
	  faceNECount  += numElementsLocal;
	} 
  /* **  elements on boundary 001 (x0=0): */
  if( domLeft && !periodic[0] )
	{
#pragma omp parallel for private(i1,k,node0) 
		for (i1=0;i1<NE1;i1++) {
			k=i1+faceNECount;
			node0=i1*N0t;

      out->FaceElements->Id[k]=faceStart + i1;
			out->FaceElements->Tag[k]   = 1;
			out->FaceElements->Color[k] = i1%2+4;
			out->FaceElements->Dom[k] = ELEMENT_INTERNAL;

			if (useElementsOnFace) {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0t;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0t+1;
			} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0t;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;

			}
		}
		totalNECount+=NE1;
		faceNECount+=NE1;
		}
		/* **  elements on boundary 002 (x0=1): */
		if( domRight && !periodic[0])
		{
#pragma omp parallel for private(i1,k,node0) 
		for (i1=0;i1<NE1;i1++) {
			k=i1+faceNECount;
			node0=(N0t-2)+i1*N0t;

      out->FaceElements->Id[k]=faceStart + numElementsLocal*2*(!periodic[1])+i1;
			out->FaceElements->Tag[k]   = 2;
			out->FaceElements->Color[k] = i1%2+6;
			out->FaceElements->Dom[k] = ELEMENT_INTERNAL;

			if (useElementsOnFace) {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0t+1;
					out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0t;
					out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
			} else {
					out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
					out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0t+1;
			}
		}
		totalNECount+=NE1;
		faceNECount+=NE1;      
	}
  out->FaceElements->minColor = 0;
  out->FaceElements->maxColor = 7;

	out->FaceElements->numElements = faceNECount;
	Finley_ElementFile_setDomainFlags( out->FaceElements );

  /* setup distribution info for other elements */
	Finley_ElementFile_setDomainFlags( out->ContactElements );
	Finley_ElementFile_setDomainFlags( out->Points );

	Finley_Mesh_prepareElementDistribution( out );

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

	//print_mesh_statistics( out, TRUE );

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif
  if (Finley_noError()) {
       if (! Finley_Mesh_isPrepared(out) ) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  }
  return out; 
}
#endif



