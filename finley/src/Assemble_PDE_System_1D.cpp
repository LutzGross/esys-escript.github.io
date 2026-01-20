
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
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

      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m
  and
      -(X_{k,i})_i + Y_k

  u has p.numComp components in a 1D domain. The shape functions for test and
  solution must be identical and row_NS == row_NN.

  Shape of the coefficients:

      A = p.numEqu x 1 x p.numComp x 1
      B = 1 x numEqu x p.numComp
      C = p.numEqu x 1 x p.numComp
      D = p.numEqu x p.numComp
      X = p.numEqu x 1
      Y = p.numEqu

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace finley {

void Assemble_PDE_System_1D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y)
{
    const int DIM = 1;
    bool expandedA = A.actsExpanded();
    bool expandedB = B.actsExpanded();
    bool expandedC = C.actsExpanded();
    bool expandedD = D.actsExpanded();
    bool expandedX = X.actsExpanded();
    bool expandedY = Y.actsExpanded();
    double *F_p = NULL;
    if(!p.F.isEmpty()) {
        p.F.requireWrite();
        F_p = p.F.getSampleDataRW(0);
    }
    const std::vector<double>& S(p.row_jac->BasisFunctions->S);
    const size_t len_EM_S = p.row_numShapesTotal*p.col_numShapesTotal*p.numEqu*p.numComp;
    const size_t len_EM_F = p.row_numShapesTotal*p.numEqu;

#pragma omp parallel
    {
        for (index_t color = p.elements->minColor; color <= p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (index_t e = 0; e < p.elements->numElements; e++) {
                if (p.elements->Color[e]==color) {
                    for (int isub=0; isub<p.numSub; isub++) {
                        const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                        const double *DSDX=&(p.row_jac->DSDX[INDEX5(0,0,0,isub,e,p.row_numShapesTotal,DIM,p.numQuadSub,p.numSub)]);
                        std::vector<double> EM_S(len_EM_S);
                        std::vector<double> EM_F(len_EM_F);
                        bool add_EM_F=false;
                        bool add_EM_S=false;
                        ///////////////
                        // process A //
                        ///////////////
                        if (!A.isEmpty()) {
                            const double *A_p=A.getSampleDataRO(e);
                            add_EM_S=true;
                            if (expandedA) {
                                const double *A_q=&(A_p[INDEX6(0,0,0,0,0,isub,p.numEqu,DIM,p.numComp,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                double f=0.;
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*
                                                        A_q[INDEX5(k,0,m,0,q,p.numEqu,DIM,p.numComp,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                                }
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                            }
                                        }
                                    }
                                }
                            } else { // constant A
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)] +=
                                                    f*A_p[INDEX4(k,0,m,0,p.numEqu,DIM,p.numComp)];
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
                            const double *B_p=B.getSampleDataRO(e);
                            add_EM_S=true;
                            if (expandedB) {
                                const double *B_q=&(B_p[INDEX5(0,0,0,0,isub,p.numEqu,DIM,p.numComp,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                double f=0.;
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)] *
                                                        B_q[INDEX4(k,0,m,q,p.numEqu,DIM,p.numComp)]*S[INDEX2(r,q,p.row_numShapes)];
                                                }
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                            }
                                        }
                                    }
                                }
                            } else { // constant B
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*S[INDEX2(r,q,p.row_numShapes)];
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)] +=
                                                    f*B_p[INDEX3(k,0,m,p.numEqu,DIM)];
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
                            const double *C_p=C.getSampleDataRO(e);
                            add_EM_S=true;
                            if (expandedC) {
                                const double *C_q=&(C_p[INDEX5(0,0,0,0,isub,p.numEqu,p.numComp,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                double f=0.;
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)] *
                                                        C_q[INDEX4(k,m,0,q,p.numEqu,p.numComp,DIM)] *
                                                        DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                                }
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                            }
                                        }
                                    }
                                }
                            } else { // constant C
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)] +=
                                                    f*C_p[INDEX3(k,m,0,p.numEqu,p.numComp)];
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
                            const double *D_p=D.getSampleDataRO(e);
                            add_EM_S=true;
                            if (expandedD) {
                                const double *D_q=&(D_p[INDEX4(0,0,0,isub,p.numEqu,p.numComp,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                double f=0.;
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[INDEX3(k,m,q,p.numEqu,p.numComp)]*S[INDEX2(r,q,p.row_numShapes)];
                                                }
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                            }
                                        }
                                    }
                                }
                            } else { // constant D
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f*D_p[INDEX2(k,m,p.numEqu)];
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
                            const double *X_p=X.getSampleDataRO(e);
                            add_EM_F=true;
                            if (expandedX) {
                                const double *X_q=&(X_p[INDEX4(0,0,0,isub,p.numEqu,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int k=0; k<p.numEqu; k++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*X_q[INDEX3(k,0,q,p.numEqu,DIM)];
                                        EM_F[INDEX2(k,s,p.numEqu)]+=f;
                                    }
                                }
                            } else { // constant X
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0.;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                    for (int k=0; k<p.numEqu; k++)
                                        EM_F[INDEX2(k,s,p.numEqu)]+=f*X_p[INDEX2(k,0,p.numEqu)];
                                }
                            }
                        }
                        ///////////////
                        // process Y //
                        ///////////////
                        if (!Y.isEmpty()) {
                            const double *Y_p=Y.getSampleDataRO(e);
                            add_EM_F=true;
                            if (expandedY) {
                                const double *Y_q=&(Y_p[INDEX3(0,0,isub,p.numEqu,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int k=0; k<p.numEqu; k++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[INDEX2(k,q,p.numEqu)];
                                        EM_F[INDEX2(k,s,p.numEqu)]+=f;
                                    }
                                }
                            } else { // constant Y
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0.;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                    for (int k=0; k<p.numEqu; k++)
                                        EM_F[INDEX2(k,s,p.numEqu)]+=f*Y_p[k];
                                }
                            }
                        }
                        // add the element matrices onto the matrix and
                        // right hand side
                        std::vector<index_t> row_index(p.row_numShapesTotal);
                        for (int q=0; q<p.row_numShapesTotal; q++)
                            row_index[q]=p.row_DOF[p.elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];

                        if (add_EM_F)
                            util::addScatter(p.row_numShapesTotal,
                                    &row_index[0], p.numEqu, &EM_F[0], F_p,
                                    p.row_DOF_UpperBound);
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

} // namespace finley

