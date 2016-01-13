
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*    assembles the mass matrix in lumped form                */

/*    The coefficient D has to be defined on the integration points or not present. */

/*    lumpedMat has to be initialized before the routine is called. */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/************************************************************************************/

void Finley_Assemble_LumpedSystem(Finley_NodeFile* nodes,Finley_ElementFile* elements, escriptDataC* lumpedMat, escriptDataC* D, const bool_t useHRZ) 
{

  bool_t reducedIntegrationOrder=FALSE, expandedD;
  char error_msg[LenErrorMsg_MAX];
  Finley_Assemble_Parameters p;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK], k, e, len_EM_lumpedMat, q, s, isub;
  type_t funcspace;
  index_t color,*row_index=NULL;
  __const double *D_p=NULL;
  double *S=NULL, *EM_lumpedMat=NULL, *Vol=NULL, *lumpedMat_p=NULL;
  register double rtmp;
  register double m_t=0., diagS=0., rtmp2=0.;

  Finley_resetError();

  if (nodes==NULL || elements==NULL) return;
  if (isEmpty(lumpedMat) || isEmpty(D)) return;
  if (isEmpty(lumpedMat) && !isEmpty(D)) { 
        Finley_setError(TYPE_ERROR,"Finley_Assemble_LumpedSystem: coefficients are non-zero but no lumped matrix is given.");
        return;
  }
  funcspace=getFunctionSpaceType(D);
  /* check if all function spaces are the same */
  if (funcspace==FINLEY_ELEMENTS) {
       reducedIntegrationOrder=FALSE;
  } else if (funcspace==FINLEY_FACE_ELEMENTS)  {
       reducedIntegrationOrder=FALSE;
  } else if (funcspace==FINLEY_REDUCED_ELEMENTS) {
       reducedIntegrationOrder=TRUE;
  } else if (funcspace==FINLEY_REDUCED_FACE_ELEMENTS)  {
       reducedIntegrationOrder=TRUE;
    } else if (funcspace==FINLEY_POINTS)  {
       reducedIntegrationOrder=TRUE;
  } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_LumpedSystem: assemblage failed because of illegal function space.");
  }
  if (! Finley_noError()) return;

  /* set all parameters in p*/
  Finley_Assemble_getAssembleParameters(nodes,elements,NULL,lumpedMat, reducedIntegrationOrder, &p);
  if (! Finley_noError()) return;
 
  /* check if all function spaces are the same */
  if (! numSamplesEqual(D,p.numQuadTotal,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_LumpedSystem: sample points of coefficient D don't match (%d,%d)",p.numQuadSub,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  /*  check the dimensions: */
  if (p.numEqu==1) {
    if (!isEmpty(D)) {
       if (!isDataPointShapeEqual(D,0,dimensions)) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_LumpedSystem: coefficient D, rank 0 expected.");
       }

    }
  } else {
    if (!isEmpty(D)) {
      dimensions[0]=p.numEqu;
      if (!isDataPointShapeEqual(D,1,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_LumpedSystem: coefficient D, expected shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
  }
  if (Finley_noError()) {
    requireWrite(lumpedMat);
    lumpedMat_p=getSampleDataRW(lumpedMat,0);
    if (funcspace==FINLEY_POINTS) {
              #pragma omp parallel private(color, D_p)
	      {
		    for (color=elements->minColor;color<=elements->maxColor;color++) {
		      /*  open loop over all elements: */
		      #pragma omp for private(e) schedule(static)
		      for(e=0;e<elements->numElements;e++){
			  if (elements->Color[e]==color) {
			    D_p=getSampleDataRO(D, e);
			    if (NULL!=D_p)  Finley_Util_AddScatter(1,
                                                        &(p.row_DOF[elements->Nodes[INDEX2(0,e,p.NN)]]),
                                                        p.numEqu,
                                                        D_p,
                                                        lumpedMat_p, 
                                                        p.row_DOF_UpperBound);
			  } /* end color check */
		      } /* end element loop */
		  } /* end color loop */
	    } /* end parallel region */
    } else {  
	  
	  len_EM_lumpedMat=p.row_numShapesTotal*p.numEqu;
	  
	  expandedD=isExpanded(D);
	  S=p.row_jac->BasisFunctions->S;

	  #pragma omp parallel private(color, EM_lumpedMat, row_index, Vol, D_p, s, q, k, rtmp, diagS, m_t, isub, rtmp2)
	  {
	    EM_lumpedMat=THREAD_MEMALLOC(len_EM_lumpedMat,double);
	    row_index=THREAD_MEMALLOC(p.row_numShapesTotal,index_t);
	    if ( !Finley_checkPtr(EM_lumpedMat) && !Finley_checkPtr(row_index) ) {
		if (p.numEqu == 1) { /* single equation */
		  if (expandedD) {	/* with expanded D */	     
		      for (color=elements->minColor;color<=elements->maxColor;color++) {
			  /*  open loop over all elements: */
			  #pragma omp for private(e) schedule(static)
			  for(e=0;e<elements->numElements;e++){              
			      if (elements->Color[e]==color) {
				  for (isub=0; isub<p.numSub; isub++) {
				    Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
				    D_p=getSampleDataRO(D,e);                          
				    if (useHRZ) {

					m_t=0; /* mass of the element: m_t */
					#pragma ivdep
					for (q=0;q<p.numQuadSub;q++) m_t+=Vol[q]*D_p[INDEX2(q, isub,p.numQuadSub) ];
				
					diagS=0; /* diagonal sum: S */
					for (s=0;s<p.row_numShapes;s++) {
					    rtmp=0;
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++) {
					      rtmp2=S[INDEX2(s,q,p.row_numShapes)];
					      rtmp+=Vol[q]*D_p[INDEX2(q, isub,p.numQuadSub)]*rtmp2*rtmp2;
					    }
					    EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
					    diagS+=rtmp;
					}
					/* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					rtmp=m_t/diagS;
					#pragma ivdep
					for (s=0;s<p.row_numShapes;s++) EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=rtmp;
				
				    } else { /* row-sum lumping */
					for (s=0;s<p.row_numShapes;s++) {
					    rtmp=0;
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_p[INDEX2(q, isub,p.numQuadSub)];
					    EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
					}

				    }
				    for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
				    Finley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
				  } /* end of isub loop */ 
			    } /* end color check */	
			  } /* end element loop */
			} /* end color loop */


		  } else  {	/* with constant D */	

		      for (color=elements->minColor;color<=elements->maxColor;color++) {
			  /*  open loop over all elements: */
			  #pragma omp for private(e) schedule(static)
			  for(e=0;e<elements->numElements;e++){
			      if (elements->Color[e]==color) {
				  for (isub=0; isub<p.numSub; isub++) {
				    Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);			     
				    D_p=getSampleDataRO(D,e);                          
				    if (useHRZ) { /* HRZ lumping */
					m_t=0; /* mass of the element: m_t */
					#pragma ivdep
					for (q=0;q<p.numQuadSub;q++) m_t+=Vol[q];
					diagS=0; /* diagonal sum: S */
					for (s=0;s<p.row_numShapes;s++) {
					    rtmp=0;
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++){
					      rtmp2=S[INDEX2(s,q,p.row_numShapes)];
						rtmp+=Vol[q]*rtmp2*rtmp2;
					    }
					    EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
					    diagS+=rtmp;
					}
					/* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					rtmp=m_t/diagS*D_p[0];
					#pragma ivdep
					for (s=0;s<p.row_numShapes;s++) EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=rtmp;
				    } else { /* row-sum lumping */
					for (s=0;s<p.row_numShapes;s++) {
					    rtmp=0;
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
					    EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp*D_p[0];
					}
				    }
				    for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
				    Finley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
							      } /* end of isub loop */ 
			    } /* end color check */ 
			  } /* end element loop */
			} /* end color loop */
		      
		  }
		} else { /* system of  equation */
		  if (expandedD) { /* with expanded D */	
		      for (color=elements->minColor;color<=elements->maxColor;color++) {
			  /*  open loop over all elements: */
			  #pragma omp for private(e) schedule(static)
			  for(e=0;e<elements->numElements;e++){
			    if (elements->Color[e]==color) {
				for (isub=0; isub<p.numSub; isub++) {     
				    Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
				    D_p=getSampleDataRO(D,e);  
				  
				    if (useHRZ) { /* HRZ lumping */
					for (k=0;k<p.numEqu;k++) {
					    m_t=0; /* mass of the element: m_t */
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++) m_t+=Vol[q]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)];
							
					    diagS=0; /* diagonal sum: S */
					    for (s=0;s<p.row_numShapes;s++) {
						rtmp=0;
						#pragma ivdep
						for (q=0;q<p.numQuadSub;q++) {
						  rtmp2=S[INDEX2(s,q,p.row_numShapes)];
						  rtmp+=Vol[q]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)]*rtmp2*rtmp2;
						}
						EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
						diagS+=rtmp;
					    }
					    /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					    rtmp=m_t/diagS;
					    #pragma ivdep
					    for (s=0;s<p.row_numShapes;s++) EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=rtmp;
					  }				  
				    } else { /* row-sum lumping */
					for (s=0;s<p.row_numShapes;s++) {
					    for (k=0;k<p.numEqu;k++) {
					      rtmp=0.;
					      #pragma ivdep
					      for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)];
					      EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
					    }
					}
				    }
				    for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
				    Finley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
				} /* end of isub loop */
			    } /* end color check */
			  } /* end element loop */
		      } /* end color loop */
		  } else { /* with constant D */
		      for (color=elements->minColor;color<=elements->maxColor;color++) {
			  /*  open loop over all elements: */
			  #pragma omp for private(e) schedule(static)
			  for(e=0;e<elements->numElements;e++){
			    if (elements->Color[e]==color) {
				for (isub=0; isub<p.numSub; isub++) {     
				    Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
				    D_p=getSampleDataRO(D,e);
				
				    if (useHRZ) { /* HRZ lumping */
					    m_t=0; /* mass of the element: m_t */
					    #pragma ivdep
					    for (q=0;q<p.numQuadSub;q++) m_t+=Vol[q]; 
					    diagS=0; /* diagonal sum: S */
					    for (s=0;s<p.row_numShapes;s++) {
						rtmp=0;
						#pragma ivdep
						for (q=0;q<p.numQuadSub;q++) {
						  rtmp2=S[INDEX2(s,q,p.row_numShapes)];
						  rtmp+=Vol[q]*rtmp2*rtmp2;
						}
						#pragma ivdep
						for (k=0;k<p.numEqu;k++) EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
						diagS+=rtmp;
					    }
					    
					    /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					    rtmp=m_t/diagS;
					    for (s=0;s<p.row_numShapes;s++) {
						  #pragma ivdep
						  for (k=0;k<p.numEqu;k++) EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=rtmp*D_p[k];
					    }
				    } else { /* row-sum lumping */
				      for (s=0;s<p.row_numShapes;s++) {
					  for (k=0;k<p.numEqu;k++) {
					      rtmp=0.;
					      #pragma ivdep
					      for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
					      EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp*D_p[k];
					  }
				      }
				    }
				    for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
				    Finley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
				} /* end of isub loop */
			    } /* end color check */
			  } /* end element loop */
		      } /* end color loop */
		  }
		}
	    } /* end of pointer check */
	    THREAD_MEMFREE(EM_lumpedMat);
	    THREAD_MEMFREE(row_index);
	  } /* end parallel region */
    }
  }
}
