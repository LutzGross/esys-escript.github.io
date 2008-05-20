
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*    assembles the mass matrix in lumped form                */

/*    The coefficient D has to be defined on the integration points or not present. */

/*    lumpedMat has to be initialized before the routine is called. */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/* Disabled until the tests pass */
/* #define NEW_LUMPING */ /* */

/**************************************************************/

void Finley_Assemble_LumpedSystem(Finley_NodeFile* nodes,Finley_ElementFile* elements, escriptDataC* lumpedMat, escriptDataC* D) 
{

  bool_t reducedIntegrationOrder=FALSE, expandedD;
  char error_msg[LenErrorMsg_MAX];
  Assemble_Parameters p;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK], k, e, len_EM_lumpedMat, q, s;
  type_t funcspace;
  index_t color,*row_index=NULL;
  double *S=NULL, *EM_lumpedMat=NULL, *Vol=NULL, *D_p=NULL, *lumpedMat_p=NULL;
  register double m_t, diagS;
  register double rtmp;
  size_t len_EM_lumpedMat_size;

  Finley_resetError();

  if (nodes==NULL || elements==NULL) return;
  if (isEmpty(lumpedMat) || isEmpty(D)) return;
  if (isEmpty(lumpedMat) && !isEmpty(D)) { 
        Finley_setError(TYPE_ERROR,"Assemble_LumpedSystem: coefficients are non-zero but no lumped matrix is given.");
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
  } else {
       Finley_setError(TYPE_ERROR,"Assemble_LumpedSystem: assemblage failed because of illegal function space.");
  }
  if (! Finley_noError()) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,NULL,lumpedMat, reducedIntegrationOrder, &p);
  if (! Finley_noError()) return;

  /* check if all function spaces are the same */

  if (! numSamplesEqual(D,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Assemble_LumpedSystem: sample points of coefficient D don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  /*  check the dimensions: */
  
  if (p.numEqu==1) {
    if (!isEmpty(D)) {
       if (!isDataPointShapeEqual(D,0,dimensions)) {
          Finley_setError(TYPE_ERROR,"Assemble_LumpedSystem: coefficient D, rank 0 expected.");
       }

    }
  } else {
    if (!isEmpty(D)) {
      dimensions[0]=p.numEqu;
      if (!isDataPointShapeEqual(D,1,dimensions)) {
          sprintf(error_msg,"Assemble_LumpedSystem: coefficient D, expected shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
  }

  if (Finley_noError()) {
    lumpedMat_p=getSampleData(lumpedMat,0);
    len_EM_lumpedMat=p.row_NN*p.numEqu;
    len_EM_lumpedMat_size=len_EM_lumpedMat*sizeof(double);
    expandedD=isExpanded(D);
    S=p.row_jac->ReferenceElement->S;

    #pragma omp parallel private(color, EM_lumpedMat, row_index, Vol, D_p, s, q, k, rtmp, m_t, diagS)
    {
       EM_lumpedMat=THREAD_MEMALLOC(len_EM_lumpedMat,double);
       row_index=THREAD_MEMALLOC(p.row_NN,index_t);
       if ( !Finley_checkPtr(EM_lumpedMat) && !Finley_checkPtr(row_index) ) {
          if (p.numEqu == 1) {
             if (expandedD) {
                 for (color=elements->minColor;color<=elements->maxColor;color++) {
                    /*  open loop over all elements: */
                    #pragma omp for private(e) schedule(static)
                    for(e=0;e<elements->numElements;e++){
                       if (elements->Color[e]==color) {
                          Vol=&(p.row_jac->volume[INDEX2(0,e,p.numQuad)]);
                          memset(EM_lumpedMat,0,len_EM_lumpedMat_size);
                          D_p=getSampleData(D,e);
                          #ifdef NEW_LUMPING /* HRZ lumping */
                          /*
                           *           Number of PDEs: 1
                           *  D_p varies over element: True
                           */
                          m_t=0; /* mass of the element: m_t */
                          for (q=0;q<p.numQuad;q++) {
                              m_t+=Vol[q]*D_p[q];
                          }                           
                          diagS=0; /* diagonal sum: S */
                          for (s=0;s<p.row_NS;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) {
                                  rtmp+=Vol[q]*D_p[q]*S[INDEX2(s,q,p.row_NS)]*S[INDEX2(s,q,p.row_NS)];
                              }
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]+=rtmp;
                              diagS+=EM_lumpedMat[INDEX2(0,s,p.numEqu)];
                          }
                          /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
                          for (s=0;s<p.row_NS;s++) {
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=m_t/diagS;
                          }
                          #else /* row-sum lumping */
                          for (s=0;s<p.row_NS;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)]*D_p[q];
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]+=rtmp;
                          }
                          #endif
                          for (q=0;q<p.row_NN;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                          Finley_Util_AddScatter(p.row_NN,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
                       } /* end color check */
                    } /* end element loop */
                  } /* end color loop */
             } else  {
                 for (color=elements->minColor;color<=elements->maxColor;color++) {
                 /*  open loop over all elements: */
                 #pragma omp for private(e) schedule(static)
                 for(e=0;e<elements->numElements;e++){
                    if (elements->Color[e]==color) {
                          Vol=&(p.row_jac->volume[INDEX2(0,e,p.numQuad)]);
                          memset(EM_lumpedMat,0,len_EM_lumpedMat_size);
                          D_p=getSampleData(D,e);
                          #ifdef NEW_LUMPING /* HRZ lumping */
                          /*
                           *           Number of PDEs: 1
                           *  D_p varies over element: False
                           */
                          m_t=0; /* mass of the element: m_t */
                          for (q=0;q<p.numQuad;q++) {
                              m_t+=Vol[q]*D_p[0];
                          }                           
                          diagS=0; /* diagonal sum: S */
                          for (s=0;s<p.row_NS;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) {
                                  rtmp+=Vol[q]*D_p[0]*S[INDEX2(s,q,p.row_NS)]*S[INDEX2(s,q,p.row_NS)];
                              }
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]+=rtmp;
                              diagS+=EM_lumpedMat[INDEX2(0,s,p.numEqu)];
                          }
                          /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
                          for (s=0;s<p.row_NS;s++) {
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=m_t/diagS;
                          }
                          #else /* row-sum lumping */
                          for (s=0;s<p.row_NS;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)];
                              EM_lumpedMat[INDEX2(0,s,p.numEqu)]+=rtmp*D_p[0];
                          }
                          #endif
                          for (q=0;q<p.row_NN;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                          Finley_Util_AddScatter(p.row_NN,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
                      } /* end color check */
                    } /* end element loop */
                 } /* end color loop */
             }
          } else {
             if (expandedD) {
                 for (color=elements->minColor;color<=elements->maxColor;color++) {
                    /*  open loop over all elements: */
                    #pragma omp for private(e) schedule(static)
                    for(e=0;e<elements->numElements;e++){
                       if (elements->Color[e]==color) {
                          Vol=&(p.row_jac->volume[INDEX2(0,e,p.numQuad)]);
                          memset(EM_lumpedMat,0,len_EM_lumpedMat_size);
                          D_p=getSampleData(D,e);
                          #ifdef NEW_LUMPING /* HRZ lumping */
                          /*
                           *           Number of PDEs: Multiple
                           *  D_p varies over element: True
                           */
                          for (k=0;k<p.numEqu;k++) {
                              m_t=0; /* mass of the element: m_t */
                              for (q=0;q<p.numQuad;q++) {
                                  m_t+=Vol[q]*D_p[INDEX2(k,q,p.numEqu)];
                              }                           
                              diagS=0; /* diagonal sum: S */
                              for (s=0;s<p.row_NS;s++) {
                                  rtmp=0;
                                  for (q=0;q<p.numQuad;q++) {
                                      rtmp+=Vol[q]*D_p[INDEX2(k,q,p.numEqu)]*S[INDEX2(s,q,p.row_NS)]*S[INDEX2(s,q,p.row_NS)];
                                  }
                                  EM_lumpedMat[INDEX2(k,s,p.numEqu)]+=rtmp;
                                  diagS+=EM_lumpedMat[INDEX2(k,s,p.numEqu)];
                              }
                              /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
                              for (s=0;s<p.row_NS;s++) {
                                  EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=m_t/diagS;
                              }
                          }
                          #else /* row-sum lumping */
                          for (s=0;s<p.row_NS;s++) {
                              for (k=0;k<p.numEqu;k++) {
                                  rtmp=0.;
                                  for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)]*D_p[INDEX2(k,q,p.numEqu)];
                                  EM_lumpedMat[INDEX2(k,s,p.numEqu)]+=rtmp;
                              }
                          }
                          #endif
                          for (q=0;q<p.row_NN;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                          Finley_Util_AddScatter(p.row_NN,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
                       } /* end color check */
                    } /* end element loop */
                } /* end color loop */
             } else {
                 /*  open loop over all elements: */
                 for (color=elements->minColor;color<=elements->maxColor;color++) {
                    #pragma omp for private(e) schedule(static)
                    for(e=0;e<elements->numElements;e++){
                       if (elements->Color[e]==color) {
                          Vol=&(p.row_jac->volume[INDEX2(0,e,p.numQuad)]);
                          memset(EM_lumpedMat,0,len_EM_lumpedMat_size);
                          D_p=getSampleData(D,e);
                          #ifdef NEW_LUMPING /* HRZ lumping */
                          /*
                           *           Number of PDEs: Multiple
                           *  D_p varies over element: False
                           */
                          for (k=0;k<p.numEqu;k++) {
                              m_t=0; /* mass of the element: m_t */
                              for (q=0;q<p.numQuad;q++) {
                                  m_t+=Vol[q]*D_p[k];
                              }                           
                              diagS=0; /* diagonal sum: S */
                              for (s=0;s<p.row_NS;s++) {
                                  rtmp=0;
                                  for (q=0;q<p.numQuad;q++) {
                                      rtmp+=Vol[q]*D_p[k]*S[INDEX2(s,q,p.row_NS)]*S[INDEX2(s,q,p.row_NS)];
                                  }
                                  EM_lumpedMat[INDEX2(k,s,p.numEqu)]+=rtmp;
                                  diagS+=EM_lumpedMat[INDEX2(k,s,p.numEqu)];
                              }
                              /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
                              for (s=0;s<p.row_NS;s++) {
                                  EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=m_t/diagS;
                              }
                          }
                          #else /* row-sum lumping */
                          for (s=0;s<p.row_NS;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)];
                              for (k=0;k<p.numEqu;k++) EM_lumpedMat[INDEX2(k,s,p.numEqu)]+=rtmp*D_p[k];
                          }
                          #endif
                          for (q=0;q<p.row_NN;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                          Finley_Util_AddScatter(p.row_NN,row_index,p.numEqu,EM_lumpedMat,lumpedMat_p, p.row_DOF_UpperBound);
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
