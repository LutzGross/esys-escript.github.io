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
/*   prepare nodes does */

/*   - creates a dense labeling of  degressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedDegressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedTo */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

void Finley_Mesh_prepareNodes(Finley_Mesh* in) {
  dim_t n,len;
  index_t id,max_id,min_id,*maskReducedDOF=NULL,*maskDOF=NULL,*reducedNodesMask=NULL,*index=NULL;
#ifdef PASO_MPI
  index_t *indexLocal=NULL, *maskReducedDOFLocation=NULL;
  index_t thisDom = in->MPIInfo->rank, i;
  dim_t bufferLength=0, numReducedLocal, numReducedInternal, numReducedBoundary, numReducedExternal, numForward, numBackward, totalDOF;

#endif

/* TODO */
/* this is ok for the automatically generated rectangular meshes, but for arbitrary meshes it
   may be neccesary to make sure that
      -  the [internal][boundary][external] ordering is respected
      -  the distribution data in Nodes->degreeOfFreedomDistribution reflects any changes in the
         DOF ordering. */

  max_id=Finley_Util_getMaxInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  min_id=Finley_Util_getMinInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  len=max_id-min_id+1;

  reducedNodesMask=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  maskDOF=TMPMEMALLOC(len,index_t);
  maskReducedDOF=TMPMEMALLOC(len,index_t);
  index=TMPMEMALLOC(MAX(in->Nodes->numNodes,len),index_t);
  if  (! (Finley_checkPtr(maskDOF) || Finley_checkPtr(maskReducedDOF) 
                                        || Finley_checkPtr(reducedNodesMask) || Finley_checkPtr(index) ) ) {
      /* initialize everything */
      #pragma omp parallel
      {
         #pragma omp for private(n) schedule(static)
         for (n=0;n<in->Nodes->numNodes;n++) {
              in->Nodes->toReduced[n]=-1;
              in->Nodes->reducedDegreeOfFreedom[n]=-1;
              reducedNodesMask[n]=-1;
         }
         #pragma omp for private(n) schedule(static)
         for (n=0;n<len;n++) {
              maskDOF[n]=-1;
              maskReducedDOF[n]=-1;
         }
      }
      /* mark all nodes used by reduced elements */
      Finley_Mesh_markNodes(reducedNodesMask,0,in,TRUE);
      
      /* mark used degrees of freedom */
      /* OMP */
      for (n=0;n<in->Nodes->numNodes;n++) {
              id=in->Nodes->degreeOfFreedom[n]-min_id;
              maskDOF[id]=1;
              if (reducedNodesMask[n]>=0) maskReducedDOF[id]=1;
      }
      /* get a list of all nodes used in the reduced mesh and convert into in->Nodes->toReduced: */
      in->Nodes->reducedNumNodes=Finley_Util_packMask(in->Nodes->numNodes,reducedNodesMask,index);

      #pragma omp parallel for private(n) schedule(static)
      for (n=0;n<in->Nodes->reducedNumNodes;n++) 
        in->Nodes->toReduced[index[n]]=n;

      /* get a list of the DOFs in the reduced mesh and convert it into reducedDegreeOfFreedom */
      in->Nodes->reducedNumDegreesOfFreedom=Finley_Util_packMask(len,maskReducedDOF,index);
      #pragma omp parallel for private(n) schedule(static)
      for (n=0;n<in->Nodes->reducedNumDegreesOfFreedom;n++) 
        maskReducedDOF[index[n]]=n;

      /* get a list of the DOFs and convert it into degreeOfFreedom */
      in->Nodes->numDegreesOfFreedom=Finley_Util_packMask(len,maskDOF,index);

      #pragma omp parallel 
      {
          #pragma omp for private(n) schedule(static)
          for (n=0;n<in->Nodes->numDegreesOfFreedom;n++) 
            maskDOF[index[n]]=n;
          #pragma omp for private(n,id) schedule(static)
          for (n=0;n<in->Nodes->numNodes;n++) {
                id=in->Nodes->degreeOfFreedom[n]-min_id;
                in->Nodes->degreeOfFreedom[n]=maskDOF[id];
                in->Nodes->reducedDegreeOfFreedom[n]=maskReducedDOF[id];
          }
      }
   }

#ifdef PASO_MPI

  /*********************************************************** 
    update the distribution data 
   ***********************************************************/

  totalDOF = in->Nodes->degreeOfFreedomDistribution->numLocal + in->Nodes->degreeOfFreedomDistribution->numExternal;

  /* update the forward and backward indices for each neighbour */
  for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numNeighbours; n++ ){
    for( i=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward; i++ )
      in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i] = maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]];
    for( i=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward; i++ )
      in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i] = maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]];
  }

  /* update the global indices for the external DOF on this subdomain */
  Finley_NodeDistribution_calculateIndexExternal( in->Nodes->degreeOfFreedomDistribution, in->Nodes->CommBuffer );

  /*********************************************************** 
    compile distribution data for the reduced degrees of freedom 
   ***********************************************************/
  
  /* determnine the number of internal, boundary and external DOF in the reduced distribution */
  totalDOF = in->Nodes->degreeOfFreedomDistribution->numLocal + in->Nodes->degreeOfFreedomDistribution->numExternal;

  maskReducedDOFLocation = MEMALLOC( in->Nodes->numDegreesOfFreedom, index_t );
  for( n=0; n<in->Nodes->numDegreesOfFreedom; n++ )
    maskReducedDOFLocation[n] = -1;
  Finley_Mesh_markOrderedDegreesOfFreedomLocation( maskReducedDOFLocation, 0, in, TRUE);

  numReducedInternal = numReducedBoundary = numReducedLocal = numReducedExternal = 0;
  for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numLocal + in->Nodes->degreeOfFreedomDistribution->numExternal; n++ ) {
    switch( maskReducedDOFLocation[n] ) {
      case 1 :
        numReducedInternal++;
        break;
      case 2 :
        numReducedBoundary++;
        break;
      case 3 :        
        numReducedExternal++;
        break;
      }
  }
  numReducedLocal = numReducedInternal + numReducedBoundary;
  MEMFREE( maskReducedDOFLocation );

  if( Finley_MPI_noError( in->MPIInfo ) )  
  {
    Finley_NodeDistribution_allocTable( in->Nodes->reducedDegreeOfFreedomDistribution, numReducedLocal, numReducedExternal, 0 );
    in->Nodes->reducedDegreeOfFreedomDistribution->numInternal = numReducedInternal;  
    in->Nodes->reducedDegreeOfFreedomDistribution->numBoundary = numReducedBoundary;

    /* determine the forward and backward distributions for each neighbour */
    bufferLength = 0;  
    for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numNeighbours; n++ ){
      if( in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward>bufferLength )
        bufferLength = in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward;
      else if( in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward>bufferLength )    
        bufferLength = in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward;
    }
    indexLocal = TMPMEMALLOC( bufferLength, index_t );

    for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numNeighbours; n++ ){
      numForward = numBackward = 0;
      for( i=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward; i++  )
        if( maskReducedDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]]>=0 )
          indexLocal[numForward++] = maskReducedDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]];      
      Finley_NodeDistribution_addForward(  in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->degreeOfFreedomDistribution->neighbours[n], numForward,  indexLocal  );
      for( i=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward; i++  )
        if( maskReducedDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]]>=0 )
          indexLocal[numBackward++] = maskReducedDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]]; 
      Finley_NodeDistribution_addBackward( in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->degreeOfFreedomDistribution->neighbours[n], numBackward, indexLocal  );
    }

    /* update the global indices for the reduced external DOF on this subdomain */

    Finley_NodeDistribution_calculateIndexExternal( in->Nodes->degreeOfFreedomDistribution, in->Nodes->CommBuffer );
    Finley_NodeDistribution_formCommBuffer( in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->reducedCommBuffer );
    Finley_NodeDistribution_calculateIndexExternal( in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->reducedCommBuffer );
  }

  TMPMEMFREE( indexLocal );
#endif
  TMPMEMFREE(reducedNodesMask);
  TMPMEMFREE(maskDOF);
  TMPMEMFREE(maskReducedDOF);
  TMPMEMFREE(index);
}

/*
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:53  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.3  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/


