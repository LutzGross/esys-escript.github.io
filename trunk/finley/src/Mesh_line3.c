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

/*   Generates numElements[0] mesh with second order elements (Line3) in the interval */
/*   [0,Length[0]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "RectangularMesh.h"



/**************************************************************
  The code for Mesh_line[2/3].c is messy, it was used as a test
  bed for getting the MPI distributed mesh figured out before
  doing the rectangle/brick domains. The code is efficient, it just
  uses too many if/else if/else statements. Unfortunately, calculating
  a distributed mesh is an easy thing to visualise, but a devil to implement!

  Sorry if I don't get around
  to neatening this code up, but rest assured, the more important
  (and complex) rectangle and brick codes are much better organised.
**************************************************************/

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
  /* some debugging printf statements */
  //printf( "rank/size = %d/%d\nNodes : %d Local, %d External[%d %d], First = %d\nElements : %d Local\nDOF : %d Local, External [%d %d], Boundary [%d %d]\nperiodicLocal [%d %d]\n\n", rank, size, *numNodesLocal, *numNodesExternal, nodesExternal[0], nodesExternal[1], *firstNode, *numElementsLocal, *numDOFLocal, DOFExternal[0], DOFExternal[1], DOFBoundary[0], DOFBoundary[1], periodicLocal[0], periodicLocal[1] );
}
#endif


Finley_Mesh* Finley_RectangularMesh_Line3(dim_t* numElements,double* Length,bool_t* periodic,index_t order,bool_t useElementsOnFace) 
#ifndef PASO_MPI
{
  dim_t N0,NE0,i0,NDOF0,NFaceElements,NUMNODES;
  index_t k,node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  N0=2*NE0+1;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements=2;
  } else {
      NDOF0=N0-1;
      NFaceElements=0;
  }
  
  /*  allocate mesh: */
  
  sprintf(name,"Rectangular mesh with %d nodes",N0);
  /* TEMPFIX */

  out=Finley_Mesh_alloc(name,1,order);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Line3,out->order);
  if (useElementsOnFace) {
    out->FaceElements=Finley_ElementFile_alloc(Line3Face,out->order);
    out->ContactElements=Finley_ElementFile_alloc(Line3Face_Contact,out->order);
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
  
  #pragma omp parallel for private(i0,k,node0) 
  for (i0=0;i0<NE0;i0++) {
    k=i0;
    node0=2*i0;

    out->Elements->Id[k]=k;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=COLOR_MOD(i0);

    out->Elements->Nodes[INDEX2(0,k,3)]=node0;
    out->Elements->Nodes[INDEX2(1,k,3)]=node0+2;
    out->Elements->Nodes[INDEX2(2,k,3)]=node0+1;
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0);
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=3;
  } else {
     NUMNODES=1;
  }
  
  if (!periodic[0]) {
     out->FaceElements->Id[0]=NE0;
     out->FaceElements->Tag[0]=1;
     out->FaceElements->Color[0]=0;
     if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=2;
       out->FaceElements->Nodes[INDEX2(2,0,NUMNODES)]=1;
     } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
     }

     out->FaceElements->Id[1]=NE0+1;
     out->FaceElements->Tag[1]=2;
     out->FaceElements->Color[1]=1;
     if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
        out->FaceElements->Nodes[INDEX2(1,1,NUMNODES)]=N0-3;
        out->FaceElements->Nodes[INDEX2(2,1,NUMNODES)]=N0-2;
     } else {
        out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
     }
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=1;

  /*  face elements done: */
  
  /*   condense the nodes: */
  
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculations:*/

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
#else
/* MPI version */
{
  dim_t N0, NE0, NE0_local, i0, NDOF0, NFaceElements, numNodesLocal, numDOFLocal, numElementsLocal, numElementsInternal, nodesExternal[2], DOFExternal[2], numNodesExternal, DOFBoundary[2];
  index_t *numForward=NULL, *numBackward=NULL;
  index_t NUMNODES, node0,k, i,firstNode=0, DOFcount=0, forwardDOF[4], backwardDOF[4];
  index_t targetDomain=-1;
  bool_t periodicLocal[2], domLeft=FALSE, domRight=FALSE, domInternal=FALSE;
  Finley_Mesh* out=NULL;
  char name[50];
  double time0=Finley_timer();
  Paso_MPIInfo *mpi_info = NULL;

  /* dimensions the global domain */
  NE0=MAX(1,numElements[0]);
  N0=2*NE0+1;
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
  domain_calculateDimension( mpi_info->rank, mpi_info->size, NE0, periodic[0], &numNodesLocal, &numDOFLocal, &numElementsLocal, &numElementsInternal, &firstNode, nodesExternal, DOFExternal, &numNodesExternal, periodicLocal, DOFBoundary );  

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

  out->Elements=Finley_ElementFile_alloc(Line3,out->order,mpi_info);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Line3Face,out->order,mpi_info);
     out->ContactElements=Finley_ElementFile_alloc(Line3Face_Contact,out->order,mpi_info);
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
  Finley_NodeDistribution_allocTable( out->Nodes->degreeOfFreedomDistribution, numDOFLocal, DOFExternal[0]+DOFExternal[1], 0 );
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

  /* this is messy, could be done cleaner */
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
        k--;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=0;
        DOFcount--;
        k++;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(N0-2)/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount++;
        k++;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(N0-3)/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount++;
        k++;
      }
      /* left hand boundary with another subdomain, need to add the nodes/DOFs that
         defines the element that spans the boundary */
      else
      {
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(firstNode-1)/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount++;
        k++;
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(firstNode-2)/DBLE(N0-1)*Length[0];
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
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=k;
        k++;
      }
      /* right hand boundary with another subdomain, need to add the nodes/DOFs that
         defines the element that spans the boundary */
      else
      {
        out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(firstNode+numNodesLocal-periodicLocal[0])/DBLE(N0-1)*Length[0];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=DOFcount;
        k++;
      }
    }
    /* setup boundary DOF data */
    if( domInternal )
    {
      targetDomain = mpi_info->rank-1;
      forwardDOF[0] = 0;
      backwardDOF[0] = numNodesLocal+1; backwardDOF[1] = numNodesLocal;
      Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, backwardDOF );

      targetDomain = mpi_info->rank+1;
      forwardDOF[0] = numNodesLocal-2; forwardDOF[1] = numNodesLocal-1;
      backwardDOF[0] = numNodesLocal+2;
      Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, forwardDOF );
      Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
    }
    else if( mpi_info->size>2 || (mpi_info->size==2&&!periodic[0]) )
      if( domLeft )
      { 
        if( periodicLocal[0] )
        {
          targetDomain = mpi_info->size-1;
          forwardDOF[0] = 0;
          backwardDOF[0] = numNodesLocal; backwardDOF[1] = numNodesLocal-1;
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, backwardDOF );
        }
        targetDomain = mpi_info->rank+1;
        forwardDOF[0] = numNodesLocal-2-periodicLocal[0]; 
        forwardDOF[1] = numNodesLocal-1-periodicLocal[0];
        backwardDOF[0] = numNodesLocal + periodicLocal[0];
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
      }
      else
      {
        targetDomain = mpi_info->rank-1;
        forwardDOF[0] = 0;
        backwardDOF[0] = numNodesLocal+1;
        backwardDOF[1] = numNodesLocal;
        Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, forwardDOF );
        Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, backwardDOF );

        if( periodicLocal[1] )
        {
          targetDomain = 0;
          forwardDOF[0] = numNodesLocal-1; forwardDOF[1] = numNodesLocal-1;         
          backwardDOF[0] = numNodesLocal+1;
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 2, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 1, backwardDOF );
        }
      }   
    else if( mpi_info->size==2 && periodic[0]) 
      if( domLeft )
      {
          targetDomain = mpi_info->size-1;
          forwardDOF[0] = 0; 
          forwardDOF[1] = numDOFLocal-2; 
          forwardDOF[2] = numDOFLocal-1;
          backwardDOF[0] = numDOFLocal+1; 
          backwardDOF[1] = numDOFLocal; 
          backwardDOF[2] = numDOFLocal+2;
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 3, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 3, backwardDOF );
      }  
      else
      {
          targetDomain = 0;          
          forwardDOF[0] = numDOFLocal-2; 
          forwardDOF[1] = numDOFLocal-1;
          forwardDOF[2] = 0; 
          backwardDOF[0] = numDOFLocal+2;
          backwardDOF[1] = numDOFLocal+1; 
          backwardDOF[2] = numDOFLocal;           
          Finley_NodeDistribution_addForward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 3, forwardDOF );
          Finley_NodeDistribution_addBackward( out->Nodes->degreeOfFreedomDistribution, targetDomain, 3, backwardDOF );
      } 
       
    
    if (! Finley_MPI_noError(mpi_info)) {
      Finley_Mesh_dealloc(out);
      return NULL;
    }
    out->Nodes->degreeOfFreedomDistribution->numBoundary = DOFBoundary[0] + DOFBoundary[1];
    out->Nodes->degreeOfFreedomDistribution->numInternal = numDOFLocal - out->Nodes->degreeOfFreedomDistribution->numBoundary;
    
    /*printf( "\n============NODES=============\n" );
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
    }*/   
  }

  /*   set the elements: */
  /*   form internal elements */ 
  //#pragma omp parallel for private(i0,k) 
  for (i0=0;i0<numElementsInternal;i0++) 
  {
    k=i0;
    node0=2*i0;
    out->Elements->Id[k]=k;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=0;

    out->Elements->Nodes[INDEX2(0,k,3)]=node0;
    out->Elements->Nodes[INDEX2(1,k,3)]=node0+2;
    out->Elements->Nodes[INDEX2(2,k,3)]=node0+1;
  }
 
  /* followed by boundary elements... */
  i0 = numElementsInternal;
  if( mpi_info->size>1 )
  {
    /* left hand boundary */
    if( mpi_info->rank>0 ) /* left hand boundary is an internal boundary */
    {
      k=i0;
      node0=numNodesLocal;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,3)]=node0+1;
      out->Elements->Nodes[INDEX2(1,k,3)]=0;
      out->Elements->Nodes[INDEX2(2,k,3)]=node0;
      i0++;
    }
    else if( periodicLocal[0] ) /* left hand boundary is a periodic boundary */
    {
      k=i0;
      node0 = numNodesLocal;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,3)]=node0+1;
      out->Elements->Nodes[INDEX2(1,k,3)]=node0-1;
      out->Elements->Nodes[INDEX2(2,k,3)]=node0;
      i0++;
    }

    /* right hand boundary */
    if( mpi_info->rank<mpi_info->size-1 ) /* right hand boundary is an internal boundary */
    {
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,3)]=numNodesLocal-2-periodicLocal[0];
      out->Elements->Nodes[INDEX2(1,k,3)]=nodesExternal[0]+numNodesLocal;
      out->Elements->Nodes[INDEX2(2,k,3)]=numNodesLocal-1-periodicLocal[0];
    }
    else if( periodicLocal[1] ) /* right hand boundary is a periodic boundary */
    {
      /* no need to do anything */;
      k=i0;
      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=0;

      out->Elements->Nodes[INDEX2(0,k,3)]=numNodesLocal-2;
      out->Elements->Nodes[INDEX2(1,k,3)]=numNodesLocal+2;
      out->Elements->Nodes[INDEX2(2,k,3)]=numNodesLocal-1;
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=0;

  out->Elements->elementDistribution->numLocal    = numElementsLocal;
  out->Elements->elementDistribution->numInternal = numElementsInternal;
  out->Elements->elementDistribution->numBoundary = numElementsLocal - numElementsInternal;
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=3;
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
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=2;
       out->FaceElements->Nodes[INDEX2(2,0,NUMNODES)]=1;
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
    if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=numNodesLocal-3;
       out->FaceElements->Nodes[INDEX2(1,domLeft,NUMNODES)]=numNodesLocal-1;
       out->FaceElements->Nodes[INDEX2(2,domLeft,NUMNODES)]=numNodesLocal-2;
    } else {
       out->FaceElements->Nodes[INDEX2(0,domLeft,NUMNODES)]=numNodesLocal-1;
    }
  }
  out->FaceElements->maxColor=0;
  out->FaceElements->minColor=0;
  out->FaceElements->elementDistribution->numBoundary = 0;
  if( domLeft || domRight && !periodic[0] )
    out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal = domLeft + domRight;
  else
    out->FaceElements->elementDistribution->numLocal = out->FaceElements->elementDistribution->numInternal = 0;

  /* setup distribution info for other elements */
  out->ContactElements->elementDistribution->numLocal = out->ContactElements->elementDistribution->numInternal = out->ContactElements->elementDistribution->numInternal = 0;
  out->Points->elementDistribution->numLocal = out->Points->elementDistribution->numInternal = out->Points->elementDistribution->numInternal = 0;

  /*printf( "\n============ELEMENTS (%d)=============\n", out->Elements->numElements );
  for( k=0; k<out->Elements->elementDistribution->numInternal; k++ )
  {
    printf( "I\tId %d : nodes [%d %d %d]->DOF [%d %d %d]\n", out->Elements->Id[k],  out->Elements->Nodes[INDEX2(0,k,3)], out->Elements->Nodes[INDEX2(1,k,3)], out->Elements->Nodes[INDEX2(2,k,3)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,3)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,3)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(2,k,3)]] );
  }
  for( k=out->Elements->elementDistribution->numInternal; k<out->Elements->elementDistribution->numLocal; k++ )
  {
    printf( "B\tId %d : nodes [%d %d %d]->DOF [%d %d %d]\n", out->Elements->Id[k],  out->Elements->Nodes[INDEX2(0,k,3)], out->Elements->Nodes[INDEX2(1,k,3)], out->Elements->Nodes[INDEX2(2,k,3)], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(0,k,3)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(1,k,3)]], out->Nodes->degreeOfFreedom[out->Elements->Nodes[INDEX2(2,k,3)]] );
  }
  printf( "\n" );
  for( k=0; k<out->FaceElements->numElements; k++ )
  {
    if( NUMNODES==3 )
				printf( "F\tId %d : nodes [%d %d %d]->DOF [%d %d %d]\n", out->FaceElements->Id[k],  out->FaceElements->Nodes[INDEX2(0,k,3)], out->FaceElements->Nodes[INDEX2(1,k,3)], out->FaceElements->Nodes[INDEX2(2,k,3)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,3)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(1,k,3)]], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(2,k,3)]] );
    else
      printf( "F\tId %d : nodes [%d]->DOF [%d]\n", out->FaceElements->Id[k], out->FaceElements->Nodes[INDEX2(0,k,1)], out->Nodes->degreeOfFreedom[out->FaceElements->Nodes[INDEX2(0,k,1)]] );
  }*/
  /*  face elements done: */
  
  /*   condense the nodes: */

  Finley_Mesh_resolveNodeIds(out);

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
  Finley_Mesh_prepare(out);

  if( !Finley_MPI_noError(mpi_info) )
  {
    if( Finley_noError() )
      Finley_setError( PASO_MPI_ERROR, "Error on another MPI process" );
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
* Revision 1.2  2005/07/08 04:07:53  jgs
* Merge of development branch back to main trunk on 2005-07-08
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

