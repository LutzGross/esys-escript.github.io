
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

      D u_m and Y

  in a 3D domain. The shape functions for test and solution must be identical
  and 2*row_NS == row_NN (contact elements).

  Shape of the coefficients:

      D = scalar
      Y = scalar

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_PDE_Single_C(const AssembleParameters& p, escript::Data& D,
                           escript::Data& Y)
{
    bool expandedD=D.actsExpanded();
    bool expandedY=Y.actsExpanded();
    p.F.requireWrite();
    double *F_p=p.F.getSampleDataRW(0);
    const double *S=p.row_jac->BasisFunctions->S;

#pragma omp parallel
    {
        std::vector<double> EM_S(p.row_numShapesTotal*p.col_numShapesTotal);
        std::vector<double> EM_F(p.row_numShapesTotal);

        for (int color=p.elements->minColor; color<=p.elements->maxColor; color++) {
            // loop over all elements
#pragma omp for
            for (int e=0; e<p.elements->numElements; e++) {
                if (p.elements->Color[e]==color) {
                    for (int isub=0; isub<p.numSub; isub++) {
                        const double *Vol=&(p.row_jac->volume[INDEX3(0,isub,e,p.numQuadSub,p.numSub)]);
                        bool add_EM_F=false;
                        bool add_EM_S=false;
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            const double *D_p=D.getSampleDataRO(e);
                            add_EM_S=true;
                            if (expandedD) {
                                const double *D_q=&(D_p[INDEX2(0,isub, p.numQuadSub)]);
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double val=0.;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            val+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[q]*S[INDEX2(r,q,p.row_numShapes)];
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]= val;
                                        EM_S[INDEX4(0,0,s,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]=-val;
                                        EM_S[INDEX4(0,0,s+p.row_numShapes,r,p.numEqu,p.numComp,p.row_numShapesTotal)]=-val;
                                        EM_S[INDEX4(0,0,s+p.row_numShapes,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]= val;
                                    }
                                }
                            } else { // constant D
                                for (int s=0; s<p.row_numShapes; s++) {
                                    for (int r=0; r<p.col_numShapes; r++) {
                                        double f=0;
                                        for (int q=0; q<p.numQuadSub; q++)
                                            f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                        const double fD=f*D_p[0];
                                        EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]= fD;
                                        EM_S[INDEX4(0,0,s,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]=-fD;
                                        EM_S[INDEX4(0,0,s+p.row_numShapes,r,p.numEqu,p.numComp,p.row_numShapesTotal)]=-fD;
                                        EM_S[INDEX4(0,0,s+p.row_numShapes,r+p.col_numShapes,p.numEqu,p.numComp,p.row_numShapesTotal)]= fD;
                                    }
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
                                    double val=0.;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        val+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[q];
                                    EM_F[INDEX2(0,s,p.numEqu)]=-val;
                                    EM_F[INDEX2(0,s+p.row_numShapes,p.numEqu)]= val;
                                }
                            } else { // constant Y
                                for (int s=0; s<p.row_numShapes; s++) {
                                    double f=0.;
                                    for (int q=0; q<p.numQuadSub; q++)
                                        f+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                    EM_F[INDEX2(0,s,p.numEqu)] = -f*Y_p[0];
                                    EM_F[INDEX2(0,s+p.row_numShapes,p.numEqu)] = f*Y_p[0];
                                }
                            }
                        }
                        // add the element matrices onto the matrix and
                        // right hand side
                        std::vector<int> row_index(p.row_numShapesTotal);
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

