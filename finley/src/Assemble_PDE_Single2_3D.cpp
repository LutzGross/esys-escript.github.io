
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

  Assembles a single PDE into the stiffness matrix S and right hand side F

      -(A_{i,j} u_,j)_i-(B_{i} u)_i+C_{j} u_,j-D u_m  and -(X_,i)_i + Y

  in a 3D domain. The shape functions for test and solution must be identical
  and row_NS == row_NN.

  Shape of the coefficients:

      A = 3 x 3
      B = 3
      C = 3
      D = scalar
      X = 3
      Y = scalar

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_PDE_Single2_3D(AssembleParameters p, ElementFile* elements,
                             Paso_SystemMatrix* Mat, escript::Data& F,
                             escript::Data& A, escript::Data& B,
                             escript::Data& C, escript::Data& D,
                             escript::Data& X, escript::Data& Y)
{
    const int DIM = 3;
    bool expandedA=A.actsExpanded();
    bool expandedB=B.actsExpanded();
    bool expandedC=C.actsExpanded();
    bool expandedD=D.actsExpanded();
    bool expandedX=X.actsExpanded();
    bool expandedY=Y.actsExpanded();
    F.requireWrite();
    double *F_p=F.getSampleDataRW(0);
    const double *S=p.row_jac->BasisFunctions->S;
    const int len_EM_S=p.row_numShapesTotal*p.col_numShapesTotal;
    const int len_EM_F=p.row_numShapesTotal;

#pragma omp parallel
    {
        for (int color=elements->minColor; color<=elements->maxColor; color++) {
            // loop over all elements:
#pragma omp for
            for (int e=0; e<elements->numElements; e++) {
                if (elements->Color[e]==color) {
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
                                const double *A_q=&(A_p[INDEX4(0,0,0,isub, DIM,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0.;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            f+=Vol[q]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant A
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f00=0;
                                        double f01=0;
                                        double f02=0;
                                        double f10=0;
                                        double f11=0;
                                        double f12=0;
                                        double f20=0;
                                        double f21=0;
                                        double f22=0;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            const double f0=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                            f00+=f0*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f01+=f0*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                            f02+=f0*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];

                                            const double f1=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                            f10+=f1*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f11+=f1*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                            f12+=f1*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];

                                            const double f2=Vol[q]*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                                            f20+=f2*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f21+=f2*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                            f22+=f2*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f00*A_p[INDEX2(0,0,DIM)]
                                            + f01*A_p[INDEX2(0,1,DIM)]
                                            + f02*A_p[INDEX2(0,2,DIM)]
                                            + f10*A_p[INDEX2(1,0,DIM)]
                                            + f11*A_p[INDEX2(1,1,DIM)]
                                            + f12*A_p[INDEX2(1,2,DIM)]
                                            + f20*A_p[INDEX2(2,0,DIM)]
                                            + f21*A_p[INDEX2(2,1,DIM)]
                                            + f22*A_p[INDEX2(2,2,DIM)];
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
                                const double *B_q=&(B_p[INDEX3(0,0,isub,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            f+=Vol[q]*S[INDEX2(r,q,p.row_numShapes)]*
                                                 (DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(0,q,DIM)]
                                                 +DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(1,q,DIM)]
                                                 +DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(2,q,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant B
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f0=0;
                                        double f1=0;
                                        double f2=0;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            const double f=Vol[q]*S[INDEX2(r,q,p.row_numShapes)];
                                            f0+=f*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                            f1+=f*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                            f2+=f*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f0*B_p[0]+f1*B_p[1]+f2*B_p[2];
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
                                const double *C_q=&(C_p[INDEX3(0,0,isub,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*
                                                (C_q[INDEX2(0,q,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                +C_q[INDEX2(1,q,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                +C_q[INDEX2(2,q,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]);
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant C
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f0=0;
                                        double f1=0;
                                        double f2=0;
                                        for (int q=0; q<p.numQuadSub; q++) {
                                            const double f=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                            f0+=f*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                            f1+=f*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                            f2+=f*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
                                        }
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f0*C_p[0]+f1*C_p[1]+f2*C_p[2];
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
                                const double *D_q=&(D_p[INDEX2(0,isub,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[q]*S[INDEX2(r,q,p.row_numShapes)];
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f;
                                    }
                                }
                            } else { // constant D
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=f*D_p[0];
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
                                const double *X_q=&(X_p[INDEX3(0,0,isub,DIM,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0;
                                    for (int q=0; q<p.numQuadSub; q++) {
                                        f+=Vol[q]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(0,q,DIM)]
                                                 + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(1,q,DIM)]
                                                 + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(2,q,DIM)]);
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)]+=f;
                                }
                            } else { // constant X
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f0=0;
                                    double f1=0;
                                    double f2=0;
                                    for (int q=0; q<p.numQuadSub; q++) {
                                        f0+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                        f1+=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                        f2+=Vol[q]*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                                    }
                                    EM_F[INDEX2(0,s,p.numEqu)]+=f0*X_p[0]+f1*X_p[1]+f2*X_p[2];
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
                                const double *Y_q=&(Y_p[INDEX2(0,isub, p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[q];
                                    EM_F[INDEX2(0,s,p.numEqu)]+=f;
                                }
                            } else { // constant Y
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                    EM_F[INDEX2(0,s,p.numEqu)]+=f*Y_p[0];
                                }
                            }
                        }
                        // add the element matrices onto the matrix and
                        // right hand side
                        std::vector<int> row_index(p.row_numShapesTotal);
                        for (int q=0; q<p.row_numShapesTotal; q++)
                            row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];

                        if (add_EM_F)
                            util::addScatter(p.row_numShapesTotal,
                                    &row_index[0], p.numEqu, &EM_F[0], F_p,
                                    p.row_DOF_UpperBound);
                        if (add_EM_S)
                            Assemble_addToSystemMatrix(Mat,
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

