
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

/****************************************************************************

  Assembles a single PDE into the stiffness matrix S and right hand side F

      -(A_{i,j} u_,j)_i-(B_{i} u)_i+C_{j} u_,j-D u_m  and -(X_,i)_i + Y

  in a 2D domain. The shape functions for test and solution must be identical
  and row_NS == row_NN.

  Shape of the coefficients:

      A = 2 x 2
      B = 2
      C = 2
      D = scalar
      X = 2
      Y = scalar

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_PDE_Single_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y)
{
    const int DIM = 2;
    bool expandedA = A.actsExpanded();
    bool expandedB = B.actsExpanded();
    bool expandedC = C.actsExpanded();
    bool expandedD = D.actsExpanded();
    bool expandedX = X.actsExpanded();
    bool expandedY = Y.actsExpanded();
    double* F_p = NULL;
    if (!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0);
    }
    const double* S = p.shapeFns;
    const int len_EM_S = p.numShapes * p.numShapes;
    const int len_EM_F = p.numShapes;

#pragma omp parallel
    {
        std::vector<double> EM_S(len_EM_S);
        std::vector<double> EM_F(len_EM_F);
        std::vector<index_t> row_index(len_EM_F);

        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e] == color) {
                    const double vol = p.jac->absD[e] * p.jac->quadweight;
                    const double* DSDX = &p.jac->DSDX[INDEX5(0, 0, 0, 0, e, p.numShapes, DIM, p.numQuad, 1)];
                    std::fill(EM_S.begin(), EM_S.end(), 0);
                    std::fill(EM_F.begin(), EM_F.end(), 0);
                    bool add_EM_F = false;
                    bool add_EM_S = false;
                    /////////////////
                    //  process A  //
                    /////////////////
                    if (!A.isEmpty()) {
                        const double* A_p = A.getSampleDataRO(e);
                        add_EM_S = true;
                        if (expandedA) {
                            const double* A_q = &A_p[INDEX4(0, 0, 0, 0, DIM, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f = 0.;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        f += vol *
                                            (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                             A_q[INDEX3(0, 0, q, DIM, DIM)] *
                                             DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                             DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                             A_q[INDEX3(0, 1, q, DIM, DIM)] *
                                             DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
                                             DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                             A_q[INDEX3(1, 0, q, DIM, DIM)] *
                                             DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                             DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                             A_q[INDEX3(1, 1, q, DIM, DIM)] *
                                             DSDX[INDEX3(r, 1, q, p.numShapes, DIM)]);
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f00 = 0;
                                    double f01 = 0;
                                    double f10 = 0;
                                    double f11 = 0;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const double f0 = vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                        const double f1 = vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                        f00 += f0 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f01 += f0 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                        f10 += f1 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f11 += f1 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                        f00 * A_p[INDEX2(0, 0, DIM)] + f01 * A_p[INDEX2(0, 1, DIM)] +
                                        f10 * A_p[INDEX2(1, 0, DIM)] + f11 * A_p[INDEX2(1, 1, DIM)];
                                }
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const double* B_p = B.getSampleDataRO(e);
                        add_EM_S = true;
                        if (expandedB) {
                            const double* B_q = &B_p[INDEX3(0, 0, 0, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f = 0.;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        f +=
                                            vol * S[INDEX2(r, q, p.numShapes)] *
                                            (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                             B_q[INDEX2(0, q, DIM)] +
                                             DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] * B_q[INDEX2(1, q, DIM)]);
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f0 = 0;
                                    double f1 = 0;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const double f = vol * S[INDEX2(r, q, p.numShapes)];
                                        f0 += f * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                        f1 += f * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                        f0 * B_p[0] + f1 * B_p[1];
                                }
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty())
                    {
                        const double* C_p = C.getSampleDataRO(e);
                        add_EM_S = true;
                        if (expandedC) {
                            const double* C_q = &C_p[INDEX3(0, 0, 0, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f = 0;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        f += vol * S[INDEX2(s, q, p.numShapes)]*
                                            (C_q[INDEX2(0, q, DIM)] *
                                            DSDX[INDEX3(r, 0, q, p.numShapes, DIM)]
                                            + C_q[INDEX2(1, q, DIM)] *
                                            DSDX[INDEX3(r, 1, q, p.numShapes, DIM)]);
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f0 = 0;
                                    double f1 = 0;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const double f = vol * S[INDEX2(s, q, p.numShapes)];
                                        f0 += f * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f1 += f * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                    }
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                        f0 * C_p[0] + f1 * C_p[1];
                                }
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty())
                    {
                        const double* D_p = D.getSampleDataRO(e);
                        add_EM_S = true;
                        if (expandedD) {
                            const double* D_q = &D_p[INDEX2(0, 0, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f = 0;
                                    for (int q = 0; q < p.numQuad; q++)
                                        f +=
                                            vol * S[INDEX2(s, q, p.numShapes)] * D_q[q] *
                                            S[INDEX2(r, q, p.numShapes)];
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    double f = 0;
                                    for (int q = 0; q < p.numQuad; q++)
                                        f += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(r, q, p.numShapes)];
                                    EM_S[INDEX4(0, 0, s, r, p.numEqu, p.numEqu, p.numShapes)] += f * D_p[0];
                                }
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const double* X_p = X.getSampleDataRO(e);
                        add_EM_F = true;
                        if (expandedX) {
                            const double* X_q = &X_p[INDEX3(0, 0, 0, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                double f = 0.;
                                for (int q = 0; q < p.numQuad; q++) {
                                    f += vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                               X_q[INDEX2(0, q, DIM)] +
                                               DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] * X_q[INDEX2(1, q, DIM)]);
                                }
                                EM_F[INDEX2(0, s, p.numEqu)] += f;
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                double f0 = 0.;
                                double f1 = 0.;
                                for (int q = 0; q < p.numQuad; q++) {
                                    f0 += vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                    f1 += vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                }
                                EM_F[INDEX2(0, s, p.numEqu)] += f0*X_p[0] + f1*X_p[1];
                            }
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const double* Y_p = Y.getSampleDataRO(e);
                        add_EM_F = true;
                        if (expandedY) {
                            const double* Y_q = &Y_p[INDEX2(0, 0, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                double f = 0;
                                for (int q = 0; q < p.numQuad; q++)
                                    f += vol * S[INDEX2(s, q, p.numShapes)] * Y_q[q];
                                EM_F[INDEX2(0, s, p.numEqu)] += f;
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                double f = 0;
                                for (int q = 0; q < p.numQuad; q++)
                                    f += vol * S[INDEX2(s, q, p.numShapes)];
                                EM_F[INDEX2(0, s, p.numEqu)] += f * Y_p[0];
                            }
                        }
                    }
                    // add the element matrices onto the matrix and
                    // right hand side
                    for (int q = 0; q < p.numShapes; q++)
                        row_index[q] = p.DOF[p.elements->Nodes[INDEX2(q, e, p.NN)]];
                    if (add_EM_F)
                        util::addScatter(p.numShapes, &row_index[0],
                                    p.numEqu, &EM_F[0], F_p, p.DOF_UpperBound);
                    if (add_EM_S)
                        Assemble_addToSystemMatrix(p.S, row_index, p.numEqu,
                                                   EM_S);
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel region
}

} // namespace dudley

