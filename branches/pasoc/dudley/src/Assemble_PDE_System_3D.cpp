
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

  Assembles the system of numEqu PDEs into the stiffness matrix S and right
  hand side F

      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m
  and
      -(X_{k,i})_i + Y_k

  u has p.numEqu components in a 3D domain. The shape functions for test and
  solution must be identical and row_NS == row_NN

  Shape of the coefficients:

      A = p.numEqu x 3 x p.numEqu x 3
      B = 3 x p.numEqu x p.numEqu
      C = p.numEqu x 3 x p.numEqu
      D = p.numEqu x p.numEqu
      X = p.numEqu x 3
      Y = p.numEqu

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

template<typename Scalar>
void Assemble_PDE_System_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y)
{
    const int DIM = 3;
    bool expandedA = A.actsExpanded();
    bool expandedB = B.actsExpanded();
    bool expandedC = C.actsExpanded();
    bool expandedD = D.actsExpanded();
    bool expandedX = X.actsExpanded();
    bool expandedY = Y.actsExpanded();
    const Scalar zero = static_cast<Scalar>(0);
    Scalar* F_p = NULL;
    if (!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0, zero);
    }
    const double* S = p.shapeFns;
    const size_t len_EM_S = p.numShapes * p.numShapes * p.numEqu * p.numEqu;
    const size_t len_EM_F = p.numShapes * p.numEqu;

#pragma omp parallel
    {
        std::vector<Scalar> EM_S(len_EM_S);
        std::vector<Scalar> EM_F(len_EM_F);
        std::vector<index_t> row_index(p.numShapes);

        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e] == color) {
                    const double vol = p.jac->absD[e] * p.jac->quadweight;
                    const double* DSDX = &p.jac->DSDX[INDEX5(0, 0, 0, 0, e, p.numShapes, DIM, p.numQuad, 1)];
                    std::fill(EM_S.begin(), EM_S.end(), zero);
                    std::fill(EM_F.begin(), EM_F.end(), zero);
                    bool add_EM_F = false;
                    bool add_EM_S = false;

                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const Scalar* A_p = A.getSampleDataRO(e, zero);
                        add_EM_S = true;
                        if (expandedA) {
                            const Scalar* A_q = &A_p[INDEX6(0, 0, 0, 0, 0, 0, p.numEqu, DIM, p.numEqu, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            Scalar f = zero;
                                            for (int q = 0; q < p.numQuad; q++) {
                                                f +=
                                                    vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 0, m, 0, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 0, m, 1, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 0, m, 2, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 1, m, 0, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 1, m, 1, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 1, m, 2, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 2, m, 0, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 2, m, 1, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
                                                           DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
                                                           A_q[INDEX5(k, 2, m, 2, q, p.numEqu, DIM, p.numEqu, DIM)]
                                                           * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)]);

                                            }
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    Scalar f00 = zero;
                                    Scalar f01 = zero;
                                    Scalar f02 = zero;
                                    Scalar f10 = zero;
                                    Scalar f11 = zero;
                                    Scalar f12 = zero;
                                    Scalar f20 = zero;
                                    Scalar f21 = zero;
                                    Scalar f22 = zero;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const Scalar f0 = vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                        f00 += f0 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f01 += f0 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                        f02 += f0 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];

                                        const Scalar f1 = vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                        f10 += f1 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f11 += f1 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                        f12 += f1 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];

                                        const Scalar f2 = vol * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
                                        f20 += f2 * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f21 += f2 * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                        f22 += f2 * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];
                                    }
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                                f00 * A_p[INDEX4(k, 0, m, 0, p.numEqu, DIM, p.numEqu)] +
                                                f01 * A_p[INDEX4(k, 0, m, 1, p.numEqu, DIM, p.numEqu)] +
                                                f02 * A_p[INDEX4(k, 0, m, 2, p.numEqu, DIM, p.numEqu)] +
                                                f10 * A_p[INDEX4(k, 1, m, 0, p.numEqu, DIM, p.numEqu)] +
                                                f11 * A_p[INDEX4(k, 1, m, 1, p.numEqu, DIM, p.numEqu)] +
                                                f12 * A_p[INDEX4(k, 1, m, 2, p.numEqu, DIM, p.numEqu)] +
                                                f20 * A_p[INDEX4(k, 2, m, 0, p.numEqu, DIM, p.numEqu)] +
                                                f21 * A_p[INDEX4(k, 2, m, 1, p.numEqu, DIM, p.numEqu)] +
                                                f22 * A_p[INDEX4(k, 2, m, 2, p.numEqu, DIM, p.numEqu)];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const Scalar* B_p = B.getSampleDataRO(e, zero);
                        add_EM_S = true;
                        if (expandedB) {
                            const Scalar* B_q = &B_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, DIM, p.numEqu, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            Scalar f = zero;
                                            for (int q = 0; q < p.numQuad; q++) {
                                                f +=
                                                    vol * S[INDEX2(r, q, p.numShapes)] *
                                                    (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                                     B_q[INDEX4(k, 0, m, q, p.numEqu, DIM, p.numEqu)] +
                                                     DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                                     B_q[INDEX4(k, 1, m, q, p.numEqu, DIM, p.numEqu)] +
                                                     DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
                                                     B_q[INDEX4(k, 2, m, q, p.numEqu, DIM, p.numEqu)]);
                                            }
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                        }
                                    }
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    Scalar f0 = zero;
                                    Scalar f1 = zero;
                                    Scalar f2 = zero;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const Scalar f = vol * S[INDEX2(r, q, p.numShapes)];
                                        f0 += f * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                        f1 += f * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                        f2 += f * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
                                    }
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                                f0 * B_p[INDEX3(k, 0, m, p.numEqu, DIM)] +
                                                f1 * B_p[INDEX3(k, 1, m, p.numEqu, DIM)] +
                                                f2 * B_p[INDEX3(k, 2, m, p.numEqu, DIM)];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        const Scalar* C_p = C.getSampleDataRO(e, zero);
                        add_EM_S = true;
                        if (expandedC) {
                            const Scalar* C_q = &C_p[INDEX5(0, 0, 0, 0, 0, p.numEqu, p.numEqu, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            Scalar f = zero;
                                            for (int q = 0; q < p.numQuad; q++) {
                                                f +=
                                                    vol * S[INDEX2(s, q, p.numShapes)] *
                                                    (C_q[INDEX4(k, m, 0, q, p.numEqu, p.numEqu, DIM)] *
                                                     DSDX[INDEX3(r, 0, q, p.numShapes, DIM)] +
                                                     C_q[INDEX4(k, m, 1, q, p.numEqu, p.numEqu, DIM)] *
                                                     DSDX[INDEX3(r, 1, q, p.numShapes, DIM)] +
                                                     C_q[INDEX4(k, m, 2, q, p.numEqu, p.numEqu, DIM)] *
                                                     DSDX[INDEX3(r, 2, q, p.numShapes, DIM)]);
                                            }
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                        }
                                    }
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    Scalar f0 = zero;
                                    Scalar f1 = zero;
                                    Scalar f2 = zero;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        const Scalar f = vol * S[INDEX2(s, q, p.numShapes)];
                                        f0 += f * DSDX[INDEX3(r, 0, q, p.numShapes, DIM)];
                                        f1 += f * DSDX[INDEX3(r, 1, q, p.numShapes, DIM)];
                                        f2 += f * DSDX[INDEX3(r, 2, q, p.numShapes, DIM)];
                                    }
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++)
                                        {
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                                f0 * C_p[INDEX3(k, m, 0, p.numEqu, p.numEqu)] +
                                                f1 * C_p[INDEX3(k, m, 1, p.numEqu, p.numEqu)] +
                                                f2 * C_p[INDEX3(k, m, 2, p.numEqu, p.numEqu)];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        const Scalar* D_p = D.getSampleDataRO(e, zero);
                        add_EM_S = true;
                        if (expandedD) {
                            const Scalar* D_q = &D_p[INDEX4(0, 0, 0, 0, p.numEqu, p.numEqu, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            Scalar f = zero;
                                            for (int q = 0; q < p.numQuad; q++) {
                                                f +=
                                                    vol * S[INDEX2(s, q, p.numShapes)] *
                                                    D_q[INDEX3(k, m, q, p.numEqu, p.numEqu)] *
                                                    S[INDEX2(r, q, p.numShapes)];
                                            }
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] += f;
                                        }
                                    }
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int r = 0; r < p.numShapes; r++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuad; q++)
                                        f += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(r, q, p.numShapes)];
                                    for (int k = 0; k < p.numEqu; k++) {
                                        for (int m = 0; m < p.numEqu; m++) {
                                            EM_S[INDEX4(k, m, s, r, p.numEqu, p.numEqu, p.numShapes)] +=
                                                f * D_p[INDEX2(k, m, p.numEqu)];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const Scalar* X_p = X.getSampleDataRO(e, zero);
                        add_EM_F = true;
                        if (expandedX) {
                            const Scalar* X_q = &X_p[INDEX4(0, 0, 0, 0, p.numEqu, DIM, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int k = 0; k < p.numEqu; k++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuad; q++) {
                                        f +=
                                            vol * (DSDX[INDEX3(s, 0, q, p.numShapes, DIM)] *
                                                   X_q[INDEX3(k, 0, q, p.numEqu, DIM)] +
                                                   DSDX[INDEX3(s, 1, q, p.numShapes, DIM)] *
                                                   X_q[INDEX3(k, 1, q, p.numEqu, DIM)] +
                                                   DSDX[INDEX3(s, 2, q, p.numShapes, DIM)] *
                                                   X_q[INDEX3(k, 2, q, p.numEqu, DIM)]);
                                    }
                                    EM_F[INDEX2(k, s, p.numEqu)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                Scalar f0 = zero;
                                Scalar f1 = zero;
                                Scalar f2 = zero;
                                for (int q = 0; q < p.numQuad; q++) {
                                    f0 += vol * DSDX[INDEX3(s, 0, q, p.numShapes, DIM)];
                                    f1 += vol * DSDX[INDEX3(s, 1, q, p.numShapes, DIM)];
                                    f2 += vol * DSDX[INDEX3(s, 2, q, p.numShapes, DIM)];
                                }
                                for (int k = 0; k < p.numEqu; k++) {
                                    EM_F[INDEX2(k, s, p.numEqu)] += f0 * X_p[INDEX2(k, 0, p.numEqu)]
                                        + f1 * X_p[INDEX2(k, 1, p.numEqu)] + f2 * X_p[INDEX2(k, 2, p.numEqu)];
                                }
                            }
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const Scalar* Y_p = Y.getSampleDataRO(e, zero);
                        add_EM_F = true;
                        if (expandedY) {
                            const Scalar* Y_q = &Y_p[INDEX3(0, 0, 0, p.numEqu, p.numQuad)];
                            for (int s = 0; s < p.numShapes; s++) {
                                for (int k = 0; k < p.numEqu; k++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuad; q++)
                                        f += vol * S[INDEX2(s, q, p.numShapes)] * Y_q[INDEX2(k, q, p.numEqu)];
                                    EM_F[INDEX2(k, s, p.numEqu)] += f;
                                }
                            }
                        } else {
                            for (int s = 0; s < p.numShapes; s++) {
                                Scalar f = zero;
                                for (int q = 0; q < p.numQuad; q++)
                                    f += vol * S[INDEX2(s, q, p.numShapes)];
                                for (int k = 0; k < p.numEqu; k++)
                                    EM_F[INDEX2(k, s, p.numEqu)] += f * Y_p[k];
                            }
                        }
                    }

                    // add the element matrices onto the matrix and right
                    // hand side
                    for (int q = 0; q < p.numShapes; q++)
                        row_index[q] = p.DOF[p.elements->Nodes[INDEX2(q, e, p.NN)]];

                    if (add_EM_F)
                        util::addScatter(p.numShapes, &row_index[0], p.numEqu,
                                         &EM_F[0], F_p, p.DOF_UpperBound);
                    if (add_EM_S)
                        Assemble_addToSystemMatrix(p.S, row_index, p.numEqu,
                                                   EM_S);
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel region
}

// instantiate our two supported versions
template void Assemble_PDE_System_3D<escript::DataTypes::real_t>(
                            const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);
template void Assemble_PDE_System_3D<escript::DataTypes::cplx_t>(
                            const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

} // namespace dudley

