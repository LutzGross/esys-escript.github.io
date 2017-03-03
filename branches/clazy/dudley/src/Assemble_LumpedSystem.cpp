
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_LumpedSystem(const NodeFile* nodes, const ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ)
{
    if (!nodes || !elements || lumpedMat.isEmpty() || D.isEmpty())
        return;

    if (D.isComplex() || lumpedMat.isComplex())
    {
        throw escript::ValueError("Assemble_LumpedSystem: complex arguments not supported.");
    }
    escript::DataTypes::real_t wantreal=0;    
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
        throw escript::ValueError("Assemble_LumpedSystem: assemblage failed because of illegal function space.");
    }

    // initialize parameters
    AssembleParameters p(nodes, elements, escript::ASM_ptr(),
                         lumpedMat, reducedIntegrationOrder);

    // check if all function spaces are the same
    if (!D.numSamplesEqual(p.numQuad, elements->numElements)) {
        std::stringstream ss;
        ss << "Assemble_LumpedSystem: sample points of coefficient D "
              "don't match (" << p.numQuad << ","
           << elements->numElements << ")";
        throw escript::ValueError(ss.str());
    }

    // check the dimensions
    if (p.numEqu == 1) {
        const escript::DataTypes::ShapeType dimensions; //dummy
        if (D.getDataPointShape() != dimensions) {
            throw escript::ValueError("Assemble_LumpedSystem: coefficient D, rank 0 expected.");
        }
    } else {
        const escript::DataTypes::ShapeType dimensions(1, p.numEqu);
        if (D.getDataPointShape() != dimensions) {
            std::stringstream ss;
            ss << "Assemble_LumpedSystem: coefficient D, expected "
                  "shape (" << p.numEqu << ",)";
            throw escript::ValueError(ss.str());
        }
    }

    lumpedMat.requireWrite();
    double* lumpedMat_p = lumpedMat.getExpandedVectorReference(wantreal).data();

    if (funcspace==DUDLEY_POINTS) {
#pragma omp parallel
        {
            for (int color=elements->minColor; color<=elements->maxColor; color++) {
                // loop over all elements
#pragma omp for
                for (index_t e=0; e<elements->numElements; e++) {
                    if (elements->Color[e]==color) {
                        const double* D_p = D.getSampleDataRO(e, wantreal);
                        util::addScatter(1,
                                      &p.DOF[elements->Nodes[INDEX2(0,e,p.NN)]],
                                      p.numEqu, D_p, lumpedMat_p,
                                      p.DOF_UpperBound);
                    } // end color check
                } // end element loop
            } // end color loop
        } // end parallel region
    } else {
        bool expandedD = D.actsExpanded();
        const double *S = NULL;
        if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &S)) {
            throw DudleyException("Assemble_LumpedSystem: Unable to locate shape function.");
        }
#pragma omp parallel
        {
            std::vector<double> EM_lumpedMat(p.numShapes * p.numEqu);
            IndexVector row_index(p.numShapes);

            if (p.numEqu == 1) { // single equation
                if (expandedD) { // with expanded D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.jac->absD[e] * p.jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e, wantreal);
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
                                    row_index[q] = p.DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                util::addScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                } else { // with constant D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.jac->absD[e] * p.jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e, wantreal);
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
                                    row_index[q] = p.DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                util::addScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.DOF_UpperBound);
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
                                const double vol = p.jac->absD[e] * p.jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e, wantreal);

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
                                    row_index[q] = p.DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                util::addScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                } else { // with constant D
                    for (int color = elements->minColor; color <= elements->maxColor; color++) {
                        // loop over all elements
#pragma omp for
                        for (index_t e = 0; e < elements->numElements; e++) {
                            if (elements->Color[e] == color) {
                                const double vol = p.jac->absD[e] * p.jac->quadweight;
                                const double* D_p = D.getSampleDataRO(e, wantreal);

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
                                    row_index[q] = p.DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
                                util::addScatter(p.numShapes, &row_index[0],
                                       p.numEqu, &EM_lumpedMat[0], lumpedMat_p,
                                       p.DOF_UpperBound);
                            } // end color check
                        } // end element loop
                    } // end color loop
                }
            }
        } // end parallel region
    }
}

} // namespace dudley

