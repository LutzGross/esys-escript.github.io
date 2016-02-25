
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************

  Assembles the mass matrix in lumped form.

  The coefficient D has to be defined on the integration points or not present.
  lumpedMat has to be initialized before the routine is called.

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

namespace dudley {

void Assemble_LumpedSystem(Dudley_NodeFile* nodes, Dudley_ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ)
{
    Dudley_resetError();

    if (!nodes || !elements || lumpedMat.isEmpty() || D.isEmpty())
        return;

    const int funcspace = D.getFunctionSpace().getTypeCode();
    bool reducedIntegrationOrder;
    // check function space of D
    if (funcspace == DUDLEY_ELEMENTS) {
        reducedIntegrationOrder = false;
    } else if (funcspace == DUDLEY_FACE_ELEMENTS) {
        reducedIntegrationOrder = false;
    } else if (funcspace == DUDLEY_REDUCED_ELEMENTS) {
        reducedIntegrationOrder = true;
    } else if (funcspace == DUDLEY_REDUCED_FACE_ELEMENTS) {
        reducedIntegrationOrder = true;
    } else {
        Dudley_setError(TYPE_ERROR, "Assemble_LumpedSystem: assemblage failed because of illegal function space.");
        return;
    }

    // initialize parameters
    Assemble_Parameters p;
    Assemble_getAssembleParameters(nodes, elements, escript::ASM_ptr(),
                                      lumpedMat, reducedIntegrationOrder, &p);
    if (!Dudley_noError())
        return;

    // check if all function spaces are the same
    if (!D.numSamplesEqual(p.numQuad, elements->numElements)) {
        std::stringstream ss;
        ss << "Assemble_LumpedSystem: sample points of coefficient D "
              "don't match (" << p.numQuad << ","
           << elements->numElements << ")";
        const std::string msg(ss.str());
        Dudley_setError(TYPE_ERROR, msg.c_str());
        return;
    }

    // check the dimensions
    if (p.numEqu == 1) {
        const escript::DataTypes::ShapeType dimensions; //dummy
        if (D.getDataPointShape() != dimensions) {
            Dudley_setError(TYPE_ERROR, "Assemble_LumpedSystem: coefficient D, rank 0 expected.");
            return;
        }
    } else {
        const escript::DataTypes::ShapeType dimensions(1, p.numEqu);
        if (D.getDataPointShape() != dimensions) {
            std::stringstream ss;
            ss << "Assemble_LumpedSystem: coefficient D, expected "
                  "shape (" << p.numEqu << ",)";
            const std::string msg(ss.str());
            Dudley_setError(TYPE_ERROR, msg.c_str());
        }
    }

    lumpedMat.requireWrite();
    double* lumpedMat_p = lumpedMat.getSampleDataRW(0);
    
    if (funcspace==DUDLEY_POINTS) {
#pragma omp parallel
        {
            for (int color=elements->minColor; color<=elements->maxColor; color++) {
                // loop over all elements
#pragma omp for
                for (index_t e=0; e<elements->numElements; e++) {
                    if (elements->Color[e]==color) {
                        const double* D_p = D.getSampleDataRO(e);
                        Dudley_Util_AddScatter(1,
                                      &p.row_DOF[elements->Nodes[INDEX2(0,e,p.NN)]],
                                      p.numEqu, D_p, lumpedMat_p, 
                                      p.row_DOF_UpperBound);
                    } // end color check
                } // end element loop
            } // end color loop
        } // end parallel region
    } else {
        bool expandedD = D.actsExpanded();
        const double *S = NULL;
        if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &S))
        {
            Dudley_setError(TYPE_ERROR, "Assemble_LumpedSystem: Unable to locate shape function.");
        }
#pragma omp parallel
        {
            std::vector<double> EM_lumpedMat(p.numShapes * p.numEqu);
            std::vector<index_t> row_index(p.numShapes);

            if (p.numEqu == 1) { // single equation
                if (expandedD) { // with expanded D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e);
                                if (useHRZ) {
                                    double m_t = 0; // mass of the element
                                    for (int q = 0; q < p.numQuad; q++)
                                        m_t += vol * D_p[INDEX2(q, 0, p.numQuad)];
                                    double diagS = 0;  // diagonal sum
                                    double rtmp;
                                    for (int s = 0; s < p.numShapes; s++) {
                                        rtmp = 0.;
                                        for (int q = 0; q < p.numQuad; q++)
                                            rtmp +=
                                                vol * D_p[INDEX2(q, 0, p.numQuad)] * S[INDEX2(s, q, p.numShapes)] *
                                                S[INDEX2(s, q, p.numShapes)];
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
                                        diagS += rtmp;
                                    }
                                    // rescale diagonals by m_t/diagS to ensure
                                    // consistent mass over element
                                    rtmp = m_t / diagS;
                                    for (int s = 0; s < p.numShapes; s++)
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] *= rtmp;
                                } else { // row-sum lumping
                                    for (int s = 0; s < p.numShapes; s++) {
                                        double rtmp = 0.;
                                        for (int q = 0; q < p.numQuad; q++)
                                            rtmp += vol * S[INDEX2(s, q, p.numShapes)] * D_p[INDEX2(q, 0, p.numQuad)];
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
                                    }
                                }
                                for (int q = 0; q < p.numShapes; q++)
                                    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                Dudley_Util_AddScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.row_DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                } else { // with constant D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e);
                                if (useHRZ) { // HRZ lumping
                                    // mass of the element
                                    const double m_t = vol*p.numQuad;
                                    double diagS = 0; // diagonal sum
                                    double rtmp;
                                    for (int s = 0; s < p.numShapes; s++) {
                                        rtmp = 0.;
                                        for (int q = 0; q < p.numQuad; q++) {
                                            rtmp += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
                                        }
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
                                        diagS += rtmp;
                                    }
                                    // rescale diagonals by m_t/diagS to ensure
                                    // consistent mass over element
                                    rtmp = m_t / diagS * D_p[0];
                                    for (int s = 0; s < p.numShapes; s++)
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] *= rtmp;
                                } else { // row-sum lumping
                                    for (int s = 0; s < p.numShapes; s++) {
                                        double rtmp = 0.;
                                        for (int q = 0; q < p.numQuad; q++)
                                            rtmp += vol * S[INDEX2(s, q, p.numShapes)];
                                        EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp * D_p[0];
                                    }
                                }
                                for (int q = 0; q < p.numShapes; q++)
                                    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                Dudley_Util_AddScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.row_DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                }

            } else { // system of equations
                if (expandedD) { // with expanded D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e);

                                if (useHRZ) { // HRZ lumping
                                    for (int k = 0; k < p.numEqu; k++) {
                                        double m_t = 0; // mass of the element
                                        for (int q = 0; q < p.numQuad; q++)
                                            m_t += vol * D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)];

                                        double diagS = 0; // diagonal sum
                                        double rtmp;
                                        for (int s = 0; s < p.numShapes; s++) {
                                            rtmp = 0.;
                                            for (int q = 0; q < p.numQuad; q++)
                                                rtmp +=
                                                    vol * D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)] *
                                                    S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
                                            diagS += rtmp;
                                        }
                                        // rescale diagonals by m_t/diagS to
                                        // ensure consistent mass over element
                                        rtmp = m_t / diagS;
                                        for (int s = 0; s < p.numShapes; s++)
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] *= rtmp;
                                    }
                                } else { // row-sum lumping
                                    for (int s = 0; s < p.numShapes; s++) {
                                        for (int k = 0; k < p.numEqu; k++) {
                                            double rtmp = 0.;
                                            for (int q = 0; q < p.numQuad; q++)
                                                rtmp +=
                                                    vol * S[INDEX2(s, q, p.numShapes)] *
                                                    D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)];
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
                                        }
                                    }
                                }
                                for (int q = 0; q < p.numShapes; q++)
                                    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                Dudley_Util_AddScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.row_DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                } else { // with constant D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e);

                                if (useHRZ) { // HRZ lumping
                                    double m_t = vol * p.numQuad; // mass of the element
                                    double diagS = 0; // diagonal sum
                                    double rtmp;
                                    for (int s = 0; s < p.numShapes; s++) {
                                        rtmp = 0.;
                                        for (int q = 0; q < p.numQuad; q++)
                                            rtmp += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
                                        for (int k = 0; k < p.numEqu; k++)
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
                                        diagS += rtmp;
                                    }
                                    // rescale diagonals by m_t/diagS to ensure
                                    // consistent mass over element
                                    rtmp = m_t / diagS;
                                    for (int s = 0; s < p.numShapes; s++) {
                                        for (int k = 0; k < p.numEqu; k++)
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] *= rtmp * D_p[k];
                                    }
                                } else { // row-sum lumping
                                    for (int s = 0; s < p.numShapes; s++) {
                                        for (int k = 0; k < p.numEqu; k++) {
                                            double rtmp = 0.;
                                            for (int q = 0; q < p.numQuad; q++)
                                                rtmp += vol * S[INDEX2(s, q, p.numShapes)];
                                            EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp * D_p[k];
                                        }
                                    }
                                }
                                for (int q = 0; q < p.numShapes; q++)
                                    row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                Dudley_Util_AddScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.row_DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                }
            }
        } // end parallel region
    }
}

} // namespace dudley

