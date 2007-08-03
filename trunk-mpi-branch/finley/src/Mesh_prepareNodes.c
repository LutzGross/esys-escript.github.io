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
#if 0
  dim_t n,len;
  index_t id,max_id,min_id,*maskReducedDOF=NULL,*maskDOF=NULL,*reducedNodesMask=NULL,*index=NULL;
#ifdef PASO_MPI
  index_t *maskDOFLocation=NULL, *indexLocal=NULL, *maskReducedDOFLocation=NULL, *mask=NULL;
  index_t thisDom = in->MPIInfo->rank, i, iI, iB, iE;
  dim_t bufferLength=0, numInternal, numBoundary, numExternal, numLocal, numReducedLocal, numReducedInternal, numReducedBoundary, numReducedExternal, numForward, numBackward, totalDOF;
	Finley_NodeDistribution *distribution=NULL;

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
      for (n=0;n<in->Nodes->reducedNumDegreesOfFreedom;n++) {
           maskReducedDOF[index[n]]=n;  
      }

      /* get a list of the DOFs and convert it into degreeOfFreedom */
      in->Nodes->numDegreesOfFreedom=Finley_Util_packMask(len,maskDOF,index);
      MEMFREE(in->Nodes->degreeOfFreedomId);
      MEMFREE(in->Nodes->reducedDegreeOfFreedomId);
      in->Nodes->degreeOfFreedomId=MEMALLOC(in->Nodes->numDegreesOfFreedom,index_t);
      in->Nodes->reducedDegreeOfFreedomId=MEMALLOC(in->Nodes->reducedNumDegreesOfFreedom,index_t);
      if (! ( Finley_checkPtr(in->Nodes->degreeOfFreedomId) || Finley_checkPtr(in->Nodes->reducedDegreeOfFreedomId) ) ) {
             #pragma omp parallel 
             {
                 #pragma omp for private(n) schedule(static)
                 for (n=0;n<in->Nodes->numDegreesOfFreedom;n++) {
                   maskDOF[index[n]]=n;
                 }
                 #pragma omp for private(n,id) schedule(static)
                 for (n=0;n<in->Nodes->numNodes;n++) {
                       id=in->Nodes->degreeOfFreedom[n]-min_id;
                       in->Nodes->degreeOfFreedom[n]=maskDOF[id];
                       in->Nodes->reducedDegreeOfFreedom[n]=maskReducedDOF[id];
                       if (maskReducedDOF[id]>-1) in->Nodes->reducedDegreeOfFreedomId[maskReducedDOF[id]]=in->Nodes->Id[n];
                       if (maskDOF[id]>-1) in->Nodes->degreeOfFreedomId[maskDOF[id]]=in->Nodes->Id[n];
                 }
             }
      }
   }
#ifdef PASO_MPI
  /*********************************************************** 
    update the distribution data 
   ***********************************************************/
	
	in->Nodes->degreeOfFreedomDistribution->numInternal = 0;
	in->Nodes->degreeOfFreedomDistribution->numBoundary = 0;
	in->Nodes->degreeOfFreedomDistribution->numExternal = 0;
	mask = MEMALLOC(in->Nodes->numDegreesOfFreedom,index_t);
	for( i=0; i<in->Nodes->numNodes; i++ )
		mask[in->Nodes->degreeOfFreedom[i]] = in->Nodes->Dom[i]; 
	for( i=0; i<in->Nodes->numDegreesOfFreedom; i++ )
		switch( mask[i] ){
			case NODE_INTERNAL :
				in->Nodes->degreeOfFreedomDistribution->numInternal++;
				break; 
			case NODE_BOUNDARY :
				in->Nodes->degreeOfFreedomDistribution->numBoundary++;
				break; 
			case NODE_EXTERNAL :
				in->Nodes->degreeOfFreedomDistribution->numExternal++;
				break; 
		}
		
	in->Nodes->degreeOfFreedomDistribution->numLocal = in->Nodes->degreeOfFreedomDistribution->numInternal + in->Nodes->degreeOfFreedomDistribution->numBoundary;
	
	/* reform the vtxdist */
	distribution = in->Nodes->degreeOfFreedomDistribution;
	
	MPI_Allgather( &distribution->numLocal, 1, MPI_INT, distribution->vtxdist+1, 1, MPI_INT, in->MPIInfo->comm );

	distribution->vtxdist[0] = 0;
	for( i=2; i<in->MPIInfo->size+1; i++ )
		distribution->vtxdist[i] += distribution->vtxdist[i-1];
		
	MPI_Allreduce( &distribution->numLocal, &distribution->numGlobal, 1, MPI_INT, MPI_SUM, in->MPIInfo->comm );
	
  /* update the forward and backward indices for each neighbour */
  for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numNeighbours; n++ ){
    for( i=0, iI=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward; i++ ){
			if( maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]]>=0 )
      	in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[iI++] = maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]];
		}
		in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward = iI;
    for( i=0, iI=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward; i++ ){
			if(maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]]>=0)	 
				in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[iI++] = maskDOF[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]];
		}
		in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward = iI;
  }

  Finley_NodeDistribution_formCommBuffer( in->Nodes->degreeOfFreedomDistribution, in->Nodes->CommBuffer );
  if ( !Finley_MPI_noError( in->MPIInfo )) {
    goto clean;
  }

  /* update the global indices for the external DOF on this subdomain */
  Finley_NodeDistribution_calculateIndexExternal( in->Nodes->degreeOfFreedomDistribution, in->Nodes->CommBuffer );
  if( !Finley_MPI_noError( in->MPIInfo ) ) {
    goto clean;
  }

  /*********************************************************** 
    compile distribution data for the reduced degrees of freedom 
   ***********************************************************/
  
  /* determnine the number of internal, boundary and external DOF in the reduced distribution */
  totalDOF = in->Nodes->degreeOfFreedomDistribution->numLocal + in->Nodes->degreeOfFreedomDistribution->numExternal;

  numReducedInternal = numReducedBoundary = numReducedLocal = numReducedExternal = 0;
	n=0;
  for( n=0; n<len; n++ ) 
		if( maskReducedDOF[n]>=0 )
			switch( mask[maskDOF[n]] ){
				case NODE_INTERNAL :
					numReducedInternal++;
					break;
				case NODE_BOUNDARY :
					numReducedBoundary++;
					break;
				case NODE_EXTERNAL :
					numReducedExternal++;
					break;
			}

  numReducedLocal = numReducedInternal + numReducedBoundary;

  Finley_NodeDistribution_allocTable( in->Nodes->reducedDegreeOfFreedomDistribution, numReducedLocal, numReducedExternal, 0 );
  if( Finley_MPI_noError( in->MPIInfo ) )  
  {
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

		/* find the new mapping from full to reduced DOF */
		for( i=0; i<in->Nodes->numNodes; i++ )
			maskReducedDOF[in->Nodes->degreeOfFreedom[i]] = in->Nodes->reducedDegreeOfFreedom[i];

		/* use the mapping to map the full distribution to the reduced distribution */	
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
  	if( !Finley_MPI_noError( in->MPIInfo ) ) {
			TMPMEMFREE( indexLocal );
			goto clean;
		}
    Finley_NodeDistribution_formCommBuffer( in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->reducedCommBuffer );
  	if( !Finley_MPI_noError( in->MPIInfo ) ) {
  		TMPMEMFREE( indexLocal );
			goto clean;
		}
    Finley_NodeDistribution_calculateIndexExternal( in->Nodes->reducedDegreeOfFreedomDistribution, in->Nodes->reducedCommBuffer );
		TMPMEMFREE( indexLocal );
  }
	MEMFREE( mask );
	Finley_MPI_noError( in->MPIInfo );
clean:
#endif
  TMPMEMFREE(reducedNodesMask);
  TMPMEMFREE(maskDOF);
  TMPMEMFREE(maskReducedDOF);
  TMPMEMFREE(index);
  if (Finley_noError()) in->Nodes->isPrepared=TRUE;
#endif
}

