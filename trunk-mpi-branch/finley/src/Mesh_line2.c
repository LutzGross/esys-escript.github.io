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

/* Generates numElements[0] mesh with first order elements (Line2) in
   the interval [0,Length[0]]. order is the desired accuracy of the
   integration scheme. */


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

#endif

#ifdef PASO_MPI
Finley_Mesh* Finley_RectangularMesh_Line2_singleCPU(dim_t* numElements,double* Length,bool_t* periodic, dim_t order,bool_t useElementsOnFace, Paso_MPIInfo *mpi_info) 
#else
Finley_Mesh* Finley_RectangularMesh_Line2(dim_t* numElements,double* Length,bool_t* periodic, dim_t order,bool_t useElementsOnFace) 
#endif
{
  /* Serial/OpenMP version */
  dim_t N0,NE0,i0,NDOF0,NFaceElements;
  index_t NUMNODES,k;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  N0=NE0+1;
 
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements=2;
  } else {
      NDOF0=N0-1;
      NFaceElements=0;
  }
  
  /*  allocate mesh: */
  
  sprintf(name,"Rectangular mesh with %d nodes",N0);
#ifdef PASO_MPI
  out=Finley_Mesh_alloc(name,1,order,mpi_info);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Line2,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Line2Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line2Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Point1_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
  if (! Finley_noError()) {
        Finley_Mesh_dealloc(out);
        return NULL;
  }
  
  /*  allocate tables: */
  Finley_NodeFile_allocTable(out->Nodes,N0);
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, N0, 0, 0);
  Finley_ElementFile_allocTable(out->Elements,NE0);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, NE0, 0);
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, 0 );
#else
  out=Finley_Mesh_alloc(name,1,order);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Line2,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Line2Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Line2Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Point1,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Point1_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
  if (! Finley_noError()) {
        Finley_Mesh_dealloc(out);
        return NULL;
  }
  
  /*  allocate tables: */
  
  Finley_NodeFile_allocTable(out->Nodes,N0);
  Finley_ElementFile_allocTable(out->Elements,NE0);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */
  
  #pragma omp parallel for private(i0,k) 
  for (i0=0;i0<N0;i0++) {
     k=i0;
     out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(i0)/DBLE(N0-1)*Length[0];
     out->Nodes->Id[k]=k;
     out->Nodes->Tag[k]=0;
     out->Nodes->degreeOfFreedom[k]=(i0%NDOF0);
#ifdef PASO_MPI
		 out->Nodes->Dom[k] = NODE_INTERNAL;
#endif
  }
  if (!periodic[0]) {
     out->Nodes->Tag[0]=1;
     out->Nodes->Tag[N0-1]=2;
  }
  
  /*   set the elements: */
  
  #pragma omp parallel for private(i0,k) 
  for (i0=0;i0<NE0;i0++) {
    k=i0;
    out->Elements->Id[k]=k;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=COLOR_MOD(i0);
#ifdef PASO_MPI
		 out->Elements->Dom[k] = ELEMENT_INTERNAL;
#endif

    out->Elements->Nodes[INDEX2(0,k,2)]=i0;
    out->Elements->Nodes[INDEX2(1,k,2)]=i0+1;
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0);
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=2;
  } else {
     NUMNODES=1;
  }
  if (!periodic[0]) {
    out->FaceElements->Id[0]=NE0;
    out->FaceElements->Tag[0]=1;
    out->FaceElements->Color[0]=0;
#ifdef PASO_MPI
		 out->FaceElements->Dom[0] = ELEMENT_INTERNAL;
#endif
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=1;
    } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
    }

    out->FaceElements->Id[1]=NE0+1;
    out->FaceElements->Tag[1]=2;
    out->FaceElements->Color[1]=1;
#ifdef PASO_MPI
		 out->FaceElements->Dom[1] = ELEMENT_INTERNAL;
#endif
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
       out->FaceElements->Nodes[INDEX2(1,1,NUMNODES)]=N0-2;
    } else {
       out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
    }
  }
  out->FaceElements->maxColor=1;
  out->FaceElements->minColor=0;

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
Finley_Mesh* Finley_RectangularMesh_Line2(dim_t* numElements,double* Length,bool_t* periodic, dim_t order,bool_t useElementsOnFace) 
{
/* MPI version */
  dim_t N0, NE0, NE0_local, i0, NDOF0, NFaceElements, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, N0t, NDOF0t;
  index_t NUMNODES,k,firstNode=0, firstNodeConstruct, DOFcount=0, forwardDOF[2], backwardDOF[2];
  index_t targetDomain=-1, node0;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE, boundaryRight, boundaryLeft;
  Finley_Mesh* out=NULL;
  char name[50];
  double time0=Finley_timer();
  Paso_MPIInfo *mpi_info = NULL;

  /* dimensions the global domain */
  NE0=MAX(1,numElements[0]);
  N0=NE0+1;
  NDOF0=N0;
  if( periodic[0] )
    NDOF0--;

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError()) {
        Finley_Mesh_dealloc(out);
        return NULL;
  }
	
	/* use the serial code to generate the mesh in the 1-CPU case */
	if( mpi_info->size==1 ){
		out = Finley_RectangularMesh_Line2_singleCPU(numElements,Length,periodic,order,useElementsOnFace,mpi_info);
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

  /* form face element if the local domain contains a periodic boundary */
  NFaceElements = 0;
  if( !periodic[0] )
  {
    if(domLeft)
      NFaceElements++;
    if(domRight)
      NFaceElements++;
  }  

	boundaryLeft = !domLeft || periodicLocal[0];
	boundaryRight = !domRight || periodicLocal[1];
	N0t = numNodesLocal + boundaryRight + boundaryLeft;
	NDOF0t = numDOFLocal + boundaryRight + boundaryLeft;
	firstNodeConstruct = firstNode - boundaryLeft;
	firstNodeConstruct = firstNodeConstruct<0 ? N0-2 : firstNodeConstruct;

  /*  allocate mesh: */
  sprintf(name,"Rectangular mesh with %d nodes",N0);
  out=Finley_Mesh_alloc(name,1,order,mpi_info);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Line2,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Line2Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line2Face_Contact,out->order,mpi_info);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Point1_Contact,out->order,mpi_info);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order,mpi_info);
  if (! Finley_noError()) {
        Finley_Mesh_dealloc(out);
        return NULL;
  }
  
  /*  allocate tables: */
  
  Finley_NodeFile_allocTable(out->Nodes,numNodesLocal+numNodesExternal);
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, numNodesLocal, numNodesExternal, 0);
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  Finley_ElementDistribution_allocTable( out->Elements->elementDistribution, numElementsLocal, domRight);
  Finley_ElementDistribution_allocTable( out->FaceElements->elementDistribution, NFaceElements, 0 );
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /*  set nodes: */
  #pragma omp parallel for private(i0,k)  
  /* local nodes */
  for (i0=0;i0<N0t;i0++) {
     k=i0;
     out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE((i0+firstNodeConstruct) % N0)/DBLE(N0-1)*Length[0];
     out->Nodes->Id[k]=k;
     out->Nodes->Tag[k]=0;
     out->Nodes->degreeOfFreedom[k]=k;
		 out->Nodes->Dom[k]=NODE_INTERNAL;
  }

  /* setup boundary DOF data */
	if( boundaryLeft ){
		out->Nodes->Dom[0] = NODE_EXTERNAL;
		out->Nodes->Dom[1] = NODE_BOUNDARY;
	}
	else
		out->Nodes->Tag[0] += 1;
	if( boundaryRight ){
		out->Nodes->Dom[N0t-1] = NODE_EXTERNAL;	
		out->Nodes->Dom[N0t-2] = NODE_BOUNDARY;	
	}
	else
		out->Nodes->Tag[N0t-1] += 2;
	if( periodicLocal[0] ){
		out->Nodes->degreeOfFreedom[2] = out->Nodes->degreeOfFreedom[1];
		out->Nodes->Dom[2] = out->Nodes->Dom[1];
	}

	if( !(mpi_info->size==2 && periodicLocal[0])){
		if( boundaryLeft  ) {
			targetDomain = mpi_info->rank-1 < 0 ? mpi_info->size-1 : mpi_info->rank-1;
			forwardDOF[0] = out->Nodes->degreeOfFreedom[1];
			backwardDOF[0] = out->Nodes->degreeOfFreedom[0];
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
		}
		if( boundaryRight ) {
			targetDomain = mpi_info->rank+1 > mpi_info->size-1 ? 0 : mpi_info->rank+1;
			forwardDOF[0] = out->Nodes->degreeOfFreedom[N0t-2];
			backwardDOF[0] = out->Nodes->degreeOfFreedom[N0t-1];
			Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
			Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
		}
  } else{
		/* periodic boundary conditions with 2 domains, need to change the order in which domain 0 passes boundary data */
		targetDomain = 1;

		forwardDOF[0] = out->Nodes->degreeOfFreedom[N0t-2];
		backwardDOF[0] = out->Nodes->degreeOfFreedom[N0t-1];
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
		
		forwardDOF[0] = out->Nodes->degreeOfFreedom[1];
		backwardDOF[0] = out->Nodes->degreeOfFreedom[0];
		Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
		Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
	}
	
  if (! Finley_MPI_noError(mpi_info)) {
    Finley_Mesh_dealloc(out);
    return NULL;
  }

  /*   set the elements: */
  /*   form internal elements */ 
  #pragma omp parallel for private(i0,k) 
  for (i0=0;i0<numElementsLocal;i0++) 
  {
    k=i0;
		node0 = (periodicLocal[0] && !i0) ? 0 :  i0 + periodicLocal[0];

    out->Elements->Id[k]=(i0+firstNodeConstruct) % NE0;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=COLOR_MOD(i0);
		out->Elements->Dom[k]=ELEMENT_INTERNAL;

    out->Elements->Nodes[INDEX2(0,k,2)]=node0;
    out->Elements->Nodes[INDEX2(1,k,2)]=node0+1;
  }
	if( boundaryLeft )
		out->Elements->Dom[0] = ELEMENT_BOUNDARY;
	if( boundaryRight )
		out->Elements->Dom[numElementsLocal-1] = ELEMENT_BOUNDARY;	

  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0);

	Finley_ElementFile_setDomainFlags( out->Elements );
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=2;
  } else {
     NUMNODES=1;
  }
  if ( domLeft && !periodicLocal[0] ) 
  {
    out->FaceElements->Id[0]=NE0;
    out->FaceElements->Tag[0]=1;
    out->FaceElements->Color[0]=0;
    out->FaceElements->Dom[0]=ELEMENT_INTERNAL;
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=1;
    } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
    }
  }
  if( domRight && !periodicLocal[1])
  { 
    out->FaceElements->Id[domLeft]=NE0+1;
    out->FaceElements->Tag[domLeft]=2;
    out->FaceElements->Color[domLeft]=1;
    out->FaceElements->Dom[domLeft]=ELEMENT_INTERNAL;
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=N0t-1;
       out->FaceElements->Nodes[INDEX2(1,domLeft,NUMNODES)]=N0t-2;
    } else {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=N0t-1;
    }
  }
	out->FaceElements->numElements=NFaceElements;
  out->FaceElements->maxColor=1;
  out->FaceElements->minColor=0;
  
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

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif
  if (Finley_noError()) {
       if ( ! Finley_Mesh_isPrepared(out) ) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  }
  return out; 
}
#endif
