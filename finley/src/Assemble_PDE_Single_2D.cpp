
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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

namespace finley {

template<typename Scalar>
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
    const Scalar zero = static_cast<Scalar>(0);
    Scalar* F_p = NULL;
    if (!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0, zero);
    }
    const std::vector<double>& S(p.row_jac->BasisFunctions->S);
    const int len_EM_S = p.row_numShapesTotal*p.col_numShapesTotal;
    const int len_EM_F = p.row_numShapesTotal;

#pragma omp parallel
    {
        std::vector<Scalar> EM_S(len_EM_S);
        std::vector<Scalar> EM_F(len_EM_F);
        IndexVector row_index(len_EM_F);

        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e] == color) {
                    for (int isub = 0; isub < p.numSub; isub++) {
                        const double* Vol = &(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                        const double* DSDX = &(p.row_jac->DSDX[INDEX5(0,0,0,isub,e, p.row_numShapesTotal,DIM,p.numQuadSub,p.numSub)]);
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
                                const Scalar* A_q = &A_p[INDEX4(0,0,0,isub,DIM,DIM,p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                + DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant A
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f00 = zero;
                                        Scalar f01 = zero;
                                        Scalar f10 = zero;
                                        Scalar f11 = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            const Scalar f0 = Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                            const Scalar f1 = Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                            f00 += f0*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f01 += f0*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                            f10 += f1*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f11 += f1*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)] +=
                                              f00 * A_p[INDEX2(0,0,DIM)]
                                            + f01 * A_p[INDEX2(0,1,DIM)]
                                            + f10 * A_p[INDEX2(1,0,DIM)]
                                            + f11 * A_p[INDEX2(1,1,DIM)];
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
                                const Scalar* B_q = &B_p[INDEX3(0,0,isub,DIM,p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*S[INDEX2(r,q,p.row_numShapes)]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(0,q,DIM)]
                                                    + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(1,q,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant B
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f0 = zero;
                                        Scalar f1 = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            const Scalar f = Vol[q]*S[INDEX2(r,q,p.row_numShapes)];
                                            f0 += f * DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                            f1 += f * DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f0*B_p[0]+f1*B_p[1];
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
                                const Scalar* C_q = &C_p[INDEX3(0,0,isub,DIM,p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*(C_q[INDEX2(0,q,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                    + C_q[INDEX2(1,q,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant C
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f0 = zero;
                                        Scalar f1 = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            const Scalar f = Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                            f0 += f * DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f1 += f * DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f0*C_p[0]+f1*C_p[1];
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
                                const Scalar* D_q = &D_p[INDEX2(0, isub, p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[q]*S[INDEX2(r,q,p.row_numShapes)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant D
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f*D_p[0];
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
                                const Scalar* X_q = &X_p[INDEX3(0,0,isub,DIM,p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuadSub; q++) {
                                        f += Vol[q]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(0,q,DIM)]
                                                 + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(1,q,DIM)]);
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)] += f;
                                }
                            } else { // constant X
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    Scalar f0 = zero;
                                    Scalar f1 = zero;
                                    for (int q = 0; q < p.numQuadSub; q++) {
                                        f0 += Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                        f1 += Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)] += f0*X_p[0] + f1*X_p[1];
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
                                const Scalar* Y_q = &Y_p[INDEX2(0,isub,p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuadSub; q++) {
                                        f += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[q];
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)] += f;
                                }
                            } else { // constant Y
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuadSub; q++) {
                                        f += Vol[q] * S[INDEX2(s,q,p.row_numShapes)];
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)] += f * Y_p[0];
                                }
                            }
                        }
                        // add the element matrices onto the matrix and
                        // right hand side
                        for (int q = 0; q < p.row_numShapesTotal; q++) {
                            row_index[q] = p.row_DOF[p.elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
                        }
                        if (add_EM_F) {
                            util::addScatter(p.row_numShapesTotal,
                                    &row_index[0], p.numEqu, &EM_F[0], F_p,
                                    p.row_DOF_UpperBound);
                        }
                        if (add_EM_S)
                            Assemble_addToSystemMatrix(p.S,
                                    p.row_numShapesTotal, &row_index[0],
                                    p.numEqu, p.col_numShapesTotal,
                                    &row_index[0], p.numComp, &EM_S[0]);
                    } // end of isub
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel region
}

// instantiate our two supported versions
template void Assemble_PDE_Single_2D<escript::DataTypes::real_t>(
                            const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);
template void Assemble_PDE_Single_2D<escript::DataTypes::cplx_t>(
                            const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

} // namespace finley

