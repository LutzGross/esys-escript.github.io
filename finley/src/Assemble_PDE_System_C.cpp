
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

/****************************************************************************

  Assembles the system of numEqu PDEs into the stiffness matrix S and right
  hand side F

      D_{k,m} u_m  and Y_k

  u has p.numComp components in a 3D domain. The shape functions for test and
  solution must be identical and 2*row_NS == row_NN (contact elements).

  Shape of the coefficients:

      D = p.numEqu x p.numComp
      Y = p.numEqu

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace finley {

template<typename Scalar>
void Assemble_PDE_System_C(const AssembleParameters& p, const escript::Data& D,
                           const escript::Data& Y)
{
    bool expandedD = D.actsExpanded();
    bool expandedY = Y.actsExpanded();
    const Scalar zero = static_cast<Scalar>(0);
    Scalar* F_p = NULL;
    if (!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0, zero);
    }
    const std::vector<double>& S(p.row_jac->BasisFunctions->S);

#pragma omp parallel
    {
        std::vector<Scalar> EM_S(p.row_numShapesTotal*p.col_numShapesTotal*p.numEqu*p.numComp);
        std::vector<Scalar> EM_F(p.row_numShapesTotal*p.numEqu);
        IndexVector row_index(p.row_numShapesTotal);

        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e] == color) {
                    for (int isub = 0; isub < p.numSub; isub++) {
                        const double* Vol = &(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                        bool add_EM_F = false;
                        bool add_EM_S = false;
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            const Scalar* D_p = D.getSampleDataRO(e, zero);
                            add_EM_S = true;
                            if (expandedD) {
                                const Scalar* D_q = &D_p[INDEX4(0,0,0,isub,p.numEqu,p.numComp, p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        for (int k = 0; k < p.numEqu; k++) {
                                            for (int m = 0; m < p.numComp; m++) {
                                                Scalar val = zero;
                                                for (int q = 0; q < p.numQuadSub; q++) {
                                                    val += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[INDEX3(k,m,q,p.numEqu,p.numComp)]*S[INDEX2(r,q,p.row_numShapes)];
                                                }
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]= val;
                                                EM_S[INDEX4(k,m,s,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]=-val;
                                                EM_S[INDEX4(k,m,s+p.row_numShapes,r,p.numEqu,p.numComp,p.row_numShapesTotal)]=-val;
                                                EM_S[INDEX4(k,m,s+p.row_numShapes,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]= val;
                                            }
                                        }
                                    }
                                }
                            } else { // constant D
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int r = 0; r < p.col_numShapes; r++) {
                                        Scalar f = zero;
                                        for (int q = 0; q < p.numQuadSub; q++) {
                                            f += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        }
                                        for (int k = 0; k < p.numEqu; k++) {
                                            for (int m = 0; m < p.numComp; m++) {
                                                const Scalar fD = f * D_p[INDEX2(k,m,p.numEqu)];
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]= fD;
                                                EM_S[INDEX4(k,m,s,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]=-fD;
                                                EM_S[INDEX4(k,m,s+p.row_numShapes,r,p.numEqu,p.numComp,p.row_numShapesTotal)]=-fD;
                                                EM_S[INDEX4(k,m,s+p.row_numShapes,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]= fD;
                                            }
                                        }
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
                                const Scalar* Y_q = &Y_p[INDEX3(0, 0, isub, p.numEqu, p.numQuadSub)];
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    for (int k = 0; k < p.numEqu; k++) {
                                        Scalar val = zero;
                                        for (int q = 0; q < p.numQuadSub; q++)
                                            val += Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[INDEX2(k,q,p.numEqu)];
                                        EM_F[INDEX2(k,s,p.numEqu)] = -val;
                                        EM_F[INDEX2(k,s+p.row_numShapes,p.numEqu)] = val;
                                    }
                                }
                            } else { // constant Y
                                for (int s = 0; s < p.row_numShapes; s++) {
                                    Scalar f = zero;
                                    for (int q = 0; q < p.numQuadSub; q++) {
                                        f += Vol[q] * S[INDEX2(s,q,p.row_numShapes)];
                                    }
                                    for (int k = 0; k < p.numEqu; k++) {
                                        EM_F[INDEX2(k,s,p.numEqu)] = -f * Y_p[k];
                                        EM_F[INDEX2(k,s+p.row_numShapes,p.numEqu)] = f * Y_p[k];
                                    }
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
template void Assemble_PDE_System_C<escript::DataTypes::real_t>(
                            const AssembleParameters& p,
                            const escript::Data& D, const escript::Data& Y);
template void Assemble_PDE_System_C<escript::DataTypes::cplx_t>(
                            const AssembleParameters& p,
                            const escript::Data& D, const escript::Data& Y);

} // namespace finley

