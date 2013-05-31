
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


/****************************************************************************

  Assembles the mass matrix in lumped form

  The coefficient D has to be defined on the integration points.
  lumpedMat has to be initialized before the routine is called.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <sstream>
#include <vector>

void Finley_Assemble_LumpedSystem(NodeFile* nodes,
                                  ElementFile* elements,
                                  escriptDataC* lumpedMat, escriptDataC* D,
                                  const bool_t useHRZ)
{
    Finley_resetError();

    if (!nodes || !elements || isEmpty(lumpedMat) || isEmpty(D))
        return;

    const type_t funcspace=getFunctionSpaceType(D);
    bool_t reducedIntegrationOrder;
    // check function space of D
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
        Finley_setError(TYPE_ERROR, "Finley_Assemble_LumpedSystem: assemblage failed because of illegal function space.");
        return;
    }

    // set all parameters in p
    Finley_Assemble_Parameters p;
    Finley_Assemble_getAssembleParameters(nodes, elements, NULL, lumpedMat,
                                          reducedIntegrationOrder, &p);
    if (!Finley_noError())
        return;

    // check if all function spaces are the same
    if (!numSamplesEqual(D, p.numQuadTotal, elements->numElements) ) {
        std::stringstream ss;
        ss << "Finley_Assemble_LumpedSystem: sample points of coefficient D "
            "don't match (" << p.numQuadSub << "," << elements->numElements
            << ").";
        std::string errorMsg = ss.str();
        Finley_setError(TYPE_ERROR, errorMsg.c_str());
        return;
    }

    // check the dimensions:
    if (p.numEqu==1) {
        const int dimensions = 0; //dummy
        if (!isDataPointShapeEqual(D, 0, &dimensions)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_LumpedSystem: coefficient D, rank 0 expected.");
            return;
        }
    } else {
        const int dimensions = p.numEqu;
        if (!isDataPointShapeEqual(D, 1, &dimensions)) {
            std::stringstream ss;
            ss << "Finley_Assemble_LumpedSystem: coefficient D does not have "
                "expected shape (" << p.numEqu << ",).";
            std::string errorMsg = ss.str();
            Finley_setError(TYPE_ERROR, errorMsg.c_str());
            return;
        }
    }

    requireWrite(lumpedMat);
    double *lumpedMat_p=getSampleDataRW(lumpedMat, 0);
    if (funcspace==FINLEY_POINTS) {
#pragma omp parallel
        {
            for (int color=elements->minColor; color<=elements->maxColor; color++) {
                // loop over all elements:
#pragma omp for
                for (int e=0; e<elements->numElements; e++) {
                    if (elements->Color[e]==color) {
                        const double *D_p=getSampleDataRO(D, e);
                        Finley_Util_AddScatter(1,
                                &(p.row_DOF[elements->Nodes[INDEX2(0,e,p.NN)]]),
                                p.numEqu, D_p, lumpedMat_p,
                                p.row_DOF_UpperBound);
                    } // end color check
                } // end element loop
            } // end color loop
        } // end parallel region
    } else { // function space not points
        bool expandedD=isExpanded(D);
        const double *S=p.row_jac->BasisFunctions->S;

#pragma omp parallel
        {
            std::vector<double> EM_lumpedMat(p.row_numShapesTotal*p.numEqu);
            std::vector<int> row_index(p.row_numShapesTotal);
            if (p.numEqu == 1) { // single equation
                if (expandedD) { // with expanded D
                    for (int color=elements->minColor; color<=elements->maxColor; color++) {
                        // loop over all elements:
#pragma omp for
                        for (int e=0; e<elements->numElements; e++) {
                            if (elements->Color[e]==color) {
                                for (int isub=0; isub<p.numSub; isub++) {
                                    const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
                                    const double *D_p=getSampleDataRO(D,e);
                                    if (useHRZ) {
                                        double m_t=0; // mass of the element
                                        double diagS=0; // diagonal sum
                                        double rtmp;
                                        #pragma ivdep
                                        for (int q=0; q<p.numQuadSub; q++)
                                            m_t+=Vol[q]*D_p[INDEX2(q, isub,p.numQuadSub) ];

                                        for (int s=0; s<p.row_numShapes; s++) {
                                            rtmp=0.;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++) {
                                                const double Sq=S[INDEX2(s,q,p.row_numShapes)];
                                                rtmp+=Vol[q]*D_p[INDEX2(q, isub,p.numQuadSub)]*Sq*Sq;
                                            }
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
                                            diagS+=rtmp;
                                        }
                                        // rescale diagonals by m_t/diagS to
                                        // ensure consistent mass over element
                                        rtmp=m_t/diagS;
                                        #pragma ivdep
                                        for (int s=0; s<p.row_numShapes; s++)
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=rtmp;
                                    } else { // row-sum lumping
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            double rtmp=0.;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++)
                                                rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_p[INDEX2(q, isub,p.numQuadSub)];
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
                                        }
                                    }
                                    for (int q=0; q<p.row_numShapesTotal; q++)
                                        row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
                                    Finley_Util_AddScatter(p.row_numShapesTotal,
                                            &row_index[0], p.numEqu,
                                            &EM_lumpedMat[0], lumpedMat_p,
                                            p.row_DOF_UpperBound);
                                } // end of isub loop
                            } // end color check
                        } // end element loop
                    } // end color loop

                } else  { // with constant D
                    for (int color=elements->minColor; color<=elements->maxColor; color++) {
                        // loop over all elements:
#pragma omp for
                        for (int e=0; e<elements->numElements; e++) {
                            if (elements->Color[e]==color) {
                                for (int isub=0; isub<p.numSub; isub++) {
                                    const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
                                    const double *D_p=getSampleDataRO(D,e);
                                    if (useHRZ) { // HRZ lumping
                                        double m_t=0; // mass of the element
                                        double diagS=0; // diagonal sum
                                        double rtmp;
                                        #pragma ivdep
                                        for (int q=0; q<p.numQuadSub; q++)
                                            m_t+=Vol[q];
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            rtmp=0.;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++) {
                                                const double Sq=S[INDEX2(s,q,p.row_numShapes)];
                                                rtmp+=Vol[q]*Sq*Sq;
                                            }
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp;
                                            diagS+=rtmp;
                                        }
                                        // rescale diagonals by m_t/diagS to
                                        // ensure consistent mass over element
                                        rtmp=m_t/diagS*D_p[0];
                                        #pragma ivdep
                                        for (int s=0; s<p.row_numShapes; s++)
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]*=rtmp;
                                    } else { // row-sum lumping
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            double rtmp=0.;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++)
                                                rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                            EM_lumpedMat[INDEX2(0,s,p.numEqu)]=rtmp*D_p[0];
                                        }
                                    }
                                    for (int q=0; q<p.row_numShapesTotal; q++)
                                        row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
                                    Finley_Util_AddScatter(p.row_numShapesTotal,
                                            &row_index[0], p.numEqu,
                                            &EM_lumpedMat[0], lumpedMat_p,
                                            p.row_DOF_UpperBound);
                                } // end of isub loop
                            } // end color check
                        } // end element loop
                    } // end color loop
                }

            } else { // system of equations
                if (expandedD) { // with expanded D
                    for (int color=elements->minColor; color<=elements->maxColor; color++) {
                        // loop over all elements:
#pragma omp for
                        for (int e=0; e<elements->numElements; e++) {
                            if (elements->Color[e]==color) {
                                for (int isub=0; isub<p.numSub; isub++) {
                                    const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                                    const double *D_p=getSampleDataRO(D,e);

                                    if (useHRZ) { // HRZ lumping
                                        for (int k=0; k<p.numEqu; k++) {
                                            double m_t=0.; // mass of element
                                            double diagS=0; // diagonal sum
                                            double rtmp;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++)
                                                m_t+=Vol[q]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)];

                                            for (int s=0; s<p.row_numShapes; s++) {
                                                rtmp=0;
                                                #pragma ivdep
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    const double Sq=S[INDEX2(s,q,p.row_numShapes)];
                                                    rtmp+=Vol[q]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)]*Sq*Sq;
                                                }
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
                                                diagS+=rtmp;
                                            }
                                            // rescale diagonals by m_t/diagS
                                            // to ensure consistent mass over
                                            // element
                                            rtmp=m_t/diagS;
                                            #pragma ivdep
                                            for (int s=0; s<p.row_numShapes; s++)
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=rtmp;
                                        }
                                    } else { // row-sum lumping
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            for (int k=0; k<p.numEqu; k++) {
                                                double rtmp=0.;
                                                #pragma ivdep
                                                for (int q=0; q<p.numQuadSub; q++)
                                                    rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_p[INDEX3(k,q,isub,p.numEqu,p.numQuadSub)];
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
                                            }
                                        }
                                    }
                                    for (int q=0; q<p.row_numShapesTotal; q++)
                                        row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
                                    Finley_Util_AddScatter(p.row_numShapesTotal,
                                            &row_index[0], p.numEqu,
                                            &EM_lumpedMat[0], lumpedMat_p,
                                            p.row_DOF_UpperBound);
                                } // end of isub loop
                            } // end color check
                        } // end element loop
                    } // end color loop

                } else { // with constant D
                    for (int color=elements->minColor; color<=elements->maxColor; color++) {
                        // loop over all elements:
#pragma omp for
                        for (int e=0; e<elements->numElements; e++) {
                            if (elements->Color[e]==color) {
                                for (int isub=0; isub<p.numSub; isub++) {
                                    const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
                                    const double *D_p=getSampleDataRO(D,e);

                                    if (useHRZ) { // HRZ lumping
                                        double m_t=0.; // mass of the element
                                        double diagS=0; // diagonal sum
                                        double rtmp;
                                        #pragma ivdep
                                        for (int q=0; q<p.numQuadSub; q++)
                                            m_t+=Vol[q];
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            rtmp=0.;
                                            #pragma ivdep
                                            for (int q=0; q<p.numQuadSub; q++) {
                                                const double Sq=S[INDEX2(s,q,p.row_numShapes)];
                                                rtmp+=Vol[q]*Sq*Sq;
                                            }
                                            #pragma ivdep
                                            for (int k=0; k<p.numEqu; k++)
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp;
                                            diagS+=rtmp;
                                        }

                                        // rescale diagonals by m_t/diagS to
                                        // ensure consistent mass over element
                                        rtmp=m_t/diagS;
                                        for (int s=0; s<p.row_numShapes; s++)
                                            #pragma ivdep
                                            for (int k=0; k<p.numEqu; k++)
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]*=rtmp*D_p[k];
                                    } else { // row-sum lumping
                                        for (int s=0; s<p.row_numShapes; s++) {
                                            for (int k=0; k<p.numEqu; k++) {
                                                double rtmp=0.;
                                                #pragma ivdep
                                                for (int q=0; q<p.numQuadSub; q++)
                                                    rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                                EM_lumpedMat[INDEX2(k,s,p.numEqu)]=rtmp*D_p[k];
                                            }
                                        }
                                    }
                                    for (int q=0; q<p.row_numShapesTotal; q++)
                                        row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
                                    Finley_Util_AddScatter(p.row_numShapesTotal,
                                            &row_index[0], p.numEqu,
                                            &EM_lumpedMat[0], lumpedMat_p,
                                            p.row_DOF_UpperBound);
                                } // end of isub loop
                            } // end color check
                        } // end element loop
                    } // end color loop
                }
            } // number of equations
        } // end parallel region
    } // function space
}

