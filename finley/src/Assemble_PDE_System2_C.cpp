
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

  Assembles the system of numEqu PDEs into the stiffness matrix S and right
  hand side F

      D_{k,m} u_m  and Y_k

  u has p.numComp components in a 3D domain. The shape functions for test and
  solution must be identical and 2*row_NS == row_NN (contact elements).

  Shape of the coefficients:

      D = p.numEqu x p.numComp
      Y = p.numEqu

*****************************************************************************/

/*  Author: Lutz Gross, l.gross@uq.edu.au */

#include "Assemble.h"
#include "Util.h"

#include <vector>

void Finley_Assemble_PDE_System2_C(Finley_Assemble_Parameters p,
                                   Finley_ElementFile* elements,
                                   Paso_SystemMatrix* Mat, escriptDataC* F,
                                   escriptDataC* D, escriptDataC* Y)
{
    bool_t expandedD=isExpanded(D);
    bool_t expandedY=isExpanded(Y);
    requireWrite(F);
    double *F_p=getSampleDataRW(F,0);
    const double *S=p.row_jac->BasisFunctions->S;

#pragma omp parallel
    {
        std::vector<double> EM_S(p.row_numShapesTotal*p.col_numShapesTotal*p.numEqu*p.numComp);
        std::vector<double> EM_F(p.row_numShapesTotal*p.numEqu);

        for (int color=elements->minColor; color<=elements->maxColor; color++) {
            // loop over all elements:
#pragma omp for
            for (int e=0; e<elements->numElements; e++) {
                if (elements->Color[e]==color) {
                    const double *D_p=getSampleDataRO(D,e);
                    const double *Y_p=getSampleDataRO(Y,e);

                    for (int isub=0; isub<p.numSub; isub++) {
                        const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                        bool add_EM_F=false;
                        bool add_EM_S=false;
                        ///////////////
                        // process D //
                        ///////////////
                        if (NULL!=D_p) {
                            add_EM_S=true;
                            if (expandedD) {
                                const double *D_q=&(D_p[INDEX4(0,0,0,isub,p.numEqu,p.numComp, p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                double val=0;
                                                for (int q=0; q<p.numQuadSub; q++) {
                                                    val+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[INDEX3(k,m,q,p.numEqu,p.numComp)]*S[INDEX2(r,q,p.row_numShapes)];
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
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        for (int k=0; k<p.numEqu; k++) {
                                            for (int m=0; m<p.numComp; m++) {
                                                const double fD=f*D_p[INDEX2(k,m,p.numEqu)];
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
                        if (NULL!=Y_p) {
                            add_EM_F=true;
                            if (expandedY) {
                                const double *Y_q=&(Y_p[INDEX3(0,0,isub, p.numEqu,p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int k=0; k<p.numEqu; k++) {
                                        double val=0;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            val+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[INDEX2(k,q,p.numEqu)];
                                        EM_F[INDEX2(k,s,p.numEqu)]=-val;
                                        EM_F[INDEX2(k,s+p.row_numShapes,p.numEqu)]= val;
                                    }
                                }
                            } else { // constant Y
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                    for (int k=0; k<p.numEqu; k++) {
                                        EM_F[INDEX2(k,s,p.numEqu)]=-f*Y_p[k];
                                        EM_F[INDEX2(k,s+p.row_numShapes,p.numEqu)]=f*Y_p[k];
                                    }
                                }
                            }
                        }
                        // add the element matrices onto the matrix and
                        // right hand side
                        std::vector<int> row_index(p.row_numShapesTotal);
                        for (int q=0; q<p.row_numShapesTotal; q++)
                            row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];
               
                        if (add_EM_F)
                            Finley_Util_AddScatter(p.row_numShapesTotal,
                                    &row_index[0], p.numEqu, &EM_F[0], F_p,
                                    p.row_DOF_UpperBound);
                        if (add_EM_S)
                            Finley_Assemble_addToSystemMatrix(Mat,
                                    p.row_numShapesTotal, &row_index[0],
                                    p.numEqu, p.col_numShapesTotal,
                                    &row_index[0], p.numComp, &EM_S[0]);
                    } // end of isub
                } // end color check
            } // end element loop
        } // end color loop
    } // end parallel region
}

