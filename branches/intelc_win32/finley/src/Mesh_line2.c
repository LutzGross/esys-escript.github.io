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
  printf( "rank/size = %d/%d\nNodes : %d Local, %d External[%d %d], First = %d\nElements : %d Local\nDOF : %d Local, External [%d %d]\nperiodicLocal [%d %d]\n\n", rank, size, *numNodesLocal, *numNodesExternal, nodesExternal[0], nodesExternal[1], *firstNode, *numElementsLocal, *numDOFLocal, DOFExternal[0], DOFExternal[1], periodicLocal[0], periodicLocal[1] );
}

Finley_Mesh* Finley_RectangularMesh_Line2(dim_t* numElements,double* Length,bool_t* periodic, dim_t order,bool_t useElementsOnFace) 
#ifndef PASO_MPI
{
printf("\nDEBUG OUTPUT: %s(%d)\n: prediodic = %d", __FILE__,__LINE__,periodic[0]);
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
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=1;
    } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
    }

    out->FaceElements->Id[1]=NE0+1;
    out->FaceElements->Tag[1]=2;
    out->FaceElements->Color[1]=1;
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
       out->FaceElements->Nodes[INDEX2(1,1,NUMNODES)]=N0-2;
    } else {
       out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
    }
  }
  out->FaceElements->maxColor=1;
  out->FaceElements->minColor=0;

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
// FIXME: dump out the mesh contents at this point to debug windows
printf("\nDEBUG OUTPUT: %s(%d)\nDumping mesh output:\n", __FILE__,__LINE__);
printf("Name: %s\n",out->Name);
printf("Reference Counter: %d\n",out->reference_counter);
printf("%s %d\n", out->FaceElements->ReferenceElement->Type->Name,out->FaceElements->numElements);
    int NN=out->FaceElements->ReferenceElement->Type->numNodes;
    for (int i=0;i<out->FaceElements->numElements;i++) {
      printf("Id: %d Tag: %d Nodes:",out->FaceElements->Id[i],out->FaceElements->Tag[i]);
      for (int j=0;j<NN;j++) printf(" %d",out->Nodes->Id[out->FaceElements->Nodes[INDEX2(j,i,NN)]]);
      printf("\n");
    }
//END FIXME
  
  return out;
}
#else
{
/* Finley_Mesh* Finley_RectangularMesh_Line2(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) { */
/* MPI version */
  dim_t N0, NE0, NE0_local, i0, NDOF0, NFaceElements, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal;
  index_t NUMNODES,k,firstNode=0, DOFcount=0, forwardDOF[2], backwardDOF[2];
  index_t targetDomain=-1;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE;
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
  Finley_NodeDistibution_allocTable( out->Nodes->degreeOfFreedomDistribution, numNodesLocal, numNodesExternal, 0 );
  Finley_ElementFile_allocTable(out->Elements,numElementsLocal);
  if( NFaceElements )
    Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /*  set nodes: */
  #pragma omp parallel for private(i0,k)  
  /* local nodes */
  for (i0=0;i0<numNodesLocal;i0++) {
     k=i0;
     out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(i0+firstNode)/DBLE(N0-1)*Length[0];
     out->Nodes->Id[k]=k;
     out->Nodes->Tag[k]=0;
     out->Nodes->degreeOfFreedom[k]=i0%(numDOFLocal);
  }

  /* external nodes - do left then right hand side */
  /* the following only applies if more than one domain */
  if( mpi_info->size>1 )
  {
    DOFcount = numNodesLocal;
    k=numNodesLocal;
    if( mpi_info->rank!=0 || periodicLocal[0] )
    {
      /* left hand boundary is periodic - 
         add the nodes/DOF that define the element on the right hand boundary */
      if( periodicLocal[0] )
      {
        /* tag the left hand boundary appropriately */
        out->Nodes->Tag[0]=1;
        k--;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=2;
        out->Nodes->degreeOfFreedom[k]=0;
        DOFcount--;
        k++;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(N0-2)/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount++;
        k++;
      }
      /* left hand boundary with another subdomain, need to add the node/DOF that
         defines the element that spans the boundary */
      else
      {
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(firstNode-1)/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount++;
        k++;
      }
    }
    if( mpi_info->rank!=(mpi_info->size-1) || periodicLocal[1] )
    {
      /* right hand boundary is periodic - add the external reference to the distribution */
      if( periodicLocal[1] )
      {
        out->Nodes->Coordinates[INDEX2(0,k,1)]=Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=2;
        out->Nodes->degreeOfFreedom[k]=k;
        k++;
      }
      /* right hand boundary with another subdomain, need to add the node/DOF that
         defines the element that spans the boundary */
      else
      {
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(firstNode+numNodesLocal-periodicLocal[0])/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount;
        k++;
      }

      /* setup boundary DOF data */
      if( domInternal )
      {
        targetDomain = mpi_info->rank-1;
        forwardDOF[0] = 0;
        backwardDOF[0] = numNodesLocal;
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );

        targetDomain = mpi_info->rank+1;
        forwardDOF[0] = numNodesLocal-1;
        backwardDOF[0] = numNodesLocal+1;
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
      }
      else if( domLeft )
      { 
        if( periodicLocal[0] )
        {
          targetDomain = mpi_info->size-1;
          forwardDOF[0] = 0;
          backwardDOF[0] = numNodesLocal;          
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
        }
        targetDomain = mpi_info->rank+1;
        forwardDOF[0] = numNodesLocal-1-periodicLocal[0];
        backwardDOF[0] = numNodesLocal + periodicLocal[0];
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
      }
      else
      {
        targetDomain = mpi_info->rank-1;
        forwardDOF[0] = 0;
        backwardDOF[0] = numNodesLocal;
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );

        if( periodicLocal[1] )
        {
          targetDomain = 0;
          forwardDOF[0] = numNodesLocal-1;          
          backwardDOF[0] = numNodesLocal+1;
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
        }
      }      
       
    }
    if (! Finley_MPI_noError(mpi_info)) {
      Finley_Mesh_dealloc(out);
      return NULL;
    }

    printf( "\n============NODES=============\n" );
    for( k=0; k<numNodesLocal; k++ )
      printf( "\tI\tId %d\tDOF %d\tcoord [%g]\n", out->Nodes->Id[k], out->Nodes->degreeOfFreedom[k] , out->Nodes->Coordinates[INDEX2(0,k,1)] );
    for( k=numNodesLocal; k<numNodesLocal+numNodesExternal; k++ )
      printf( "\tE\tId %d\tDOF %d\tcoord [%g]\n", out->Nodes->Id[k], out->Nodes->degreeOfFreedom[k] , out->Nodes->Coordinates[INDEX2(0,k,1)] );
  
    for( k=0; k<out->Nodes->degreeOfFreedomDistribution->numNeighbours; k++ )
    {
      if( out->Nodes->degreeOfFreedomDistribution->neighbours[k]>=0 )
      {
        printf( "\t%d boundary DOF { ", out->Nodes->degreeOfFreedomDistribution->edges[k]->numForward ); 
        for( i0=0; i0<out->Nodes->degreeOfFreedomDistribution->edges[k]->numForward; i0++ )
          printf( "%d ", out->Nodes->degreeOfFreedomDistribution->edges[k]->indexForward[i0] );
        printf("} to %d\n", out->Nodes->degreeOfFreedomDistribution->neighbours[k] );
        printf( "\t%d boundary DOF { ", out->Nodes->degreeOfFreedomDistribution->edges[k]->numBackward ); 
        for( i0=0; i0<out->Nodes->degreeOfFreedomDistribution->edges[k]->numBackward; i0++ )
          printf( "%d ", out->Nodes->degreeOfFreedomDistribution->edges[k]->indexBackward[i0] );
        printf("} from %d\n", out->Nodes->degreeOfFreedomDistribution->neighbours[k] );
      }
    }    
  }
  else
  {
    printf( "\n============NODES=============\n" );
    for( k=0; k<numNodesLocal; k++ )
      printf( "\tI\tId %d\tDOF %d\tcoord [%g]\n", out->Nodes->Id[k], out->Nodes->degreeOfFreedom[k] , out->Nodes->Coordinates[INDEX2(0,k,1)] );
  }

  /*   set the elements: */
  /*   form internal elements */ 
  //#pragma omp parallel for private(i0,k) 
  for (i0=0;i0<numElementsInternal;i0++) 
  {
    k=i0;
    out->Elements->Id[k]=k;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=0;

    out->Elements->Nodes[INDEX2(0,k,2)]=i0;
    out->Elements->Nodes[INDEX2(1,k,2)]=i0+1;
  }
 
  /* followed by boundary elements... */
  i0 = numElementsInternal;
  if( mpi_info->size>1 )
  {
    /* left hand boundary */
    if( mpi_info->rank>0 ) /* left hand boundary is an internal boundary */
    {
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,2)]=i0+1;
      out->Elements->Nodes[INDEX2(1,k,2)]=0;
      i0++;
    }
    else if( periodicLocal[0] ) /* left hand boundary is a periodic boundary */
    {
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,2)]=i0+2;
      out->Elements->Nodes[INDEX2(1,k,2)]=i0+1;
      i0++;
    }

    /* right hand boundary */
    if( mpi_info->rank<mpi_info->size-1 ) /* right hand boundary is an internal boundary */
    {
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,2)]=numNodesLocal-1-periodicLocal[0];
      out->Elements->Nodes[INDEX2(1,k,2)]=nodesExternal[0]+numNodesLocal;
      i0++;
    }
    else if( periodicLocal[1] ) /* right hand boundary is a periodic boundary */
    {
      /* no need to do anything */;
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,2)]=numNodesLocal-1;
      out->Elements->Nodes[INDEX2(1,k,2)]=numNodesLocal+1;
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=0;

  out->Elements->elementDistribution->numLocal    = numElementsLocal;
  out->Elements->elementDistribution->numInternal = numElementsInternal;
  out->Elements->elementDistribution->numBoundary = numElementsLocal - numElementsInternal;
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=2;
  } else {
     NUMNODES=1;
  }
  if ( domLeft && !periodicLocal[0] ) 
  {
    out->FaceElements->Id[0]=i0-1;
    out->FaceElements->Tag[0]=1;
    out->FaceElements->Color[0]=0;
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=1;
    } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
    }
    i0++;
  }
  if( domRight && !periodicLocal[1])
  { 
    out->FaceElements->Id[domLeft]=i0;
    out->FaceElements->Tag[domLeft]=2;
    out->FaceElements->Color[domLeft]=1;
    /* TODO */
    /* check that the correct indices have been used */
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=numNodesLocal-2;
       out->FaceElements->Nodes[INDEX2(1,domLeft,NUMNODES)]=numNodesLocal-1;
    } else {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=numNodesLocal-1;
    }
  }
  out->FaceElements->maxColor=0;
  out->FaceElements->minColor=0;
  if( domLeft || domRight && !periodic[0] )
    out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal = domLeft + domRight;

  printf( "\n============ELEMENTS (%d)=============\n", out->Elements->numElements );
  for( k=0; k<out->Elements->elementDistribution->numInternal; k++ )
  {
    printf( "I\tId %d : nodes [%d %d]->DOF [%d %d]\n", out->Elements->Id[k], out->Elements->Nodes[INDEX2(0,k,2)], out->Elements->Nodes[INDEX2(1,k,2)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,2)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,2)]] );
  }
  for( k=out->Elements->elementDistribution->numInternal; k<out->Elements->elementDistribution->numLocal; k++ )
  {
    printf( "E\tId %d : nodes [%d %d]->DOF [%d %d]\n", out->Elements->Id[k], out->Elements->Nodes[INDEX2(0,k,2)], out->Elements->Nodes[INDEX2(1,k,2)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,2)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,2)]] );
  }
  for( k=0; k<out->FaceElements->numElements; k++ )
  {
    if( NUMNODES==2 )
      printf( "F\tId %d : nodes [%d %d]->DOF [%d %d]\n", out->FaceElements->Id[k], out->FaceElements->Nodes[INDEX2(0,k,2)], out->FaceElements->Nodes[INDEX2(1,k,2)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,2)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(1,k,2)]] );
    else
      printf( "F\tId %d : nodes [%d]->DOF [%d]\n", out->FaceElements->Id[k], out->FaceElements->Nodes[INDEX2(0,k,1)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,1)]] );
  }
  /*  face elements done: */
  
  /*   condense the nodes: */
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculatuions:*/

  /* TEMPFIX */
  /* this isn't needed for a mesh generated this way */ 
  //Finley_Mesh_prepare(out);

  if (! Finley_MPI_noError( mpi_info )) {
      Paso_MPIInfo_dealloc( mpi_info );
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  /* free up memory */
  Paso_MPIInfo_dealloc( mpi_info );

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  return out; 
}
#endif




/*
* Revision 1.3  2005/09/01 03:31:36  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-01
*
* Revision 1.2.2.2  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2.2.1  2005/08/24 02:02:18  gross
* timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
*
* Revision 1.2  2005/07/08 04:07:52  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/06/30 01:53:56  gross
* a bug in coloring fixed
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

